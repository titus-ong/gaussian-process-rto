classdef GPCreator
    properties
        system                  % Either HYSYSFile or system object
        objective               % System objective function
        linear_con_A            % System linear inequality constraints LHS (Ax = b)
        linear_con_b            % System linear inequality constraints RHS
        lineq_con_A             % System linear equality constraints LHS
        lineq_con_b             % System linear equality constraints RHS
        lb                      % System lower bounds
        ub                      % System upper bounds
        nonlin_con              % System nonlinear constraints
        options                 % System options for GP
        output_fields           % System field names for outputs
        model                   % GP model
        
        training_input          % Training input data
        mean_input              % Mean of input
        std_input               % Standard deviation of input
        norm_input              % Normalised input
        
        training_output         % Training output data
        mean_output             % Mean of output
        std_output              % Standard deviation of output
        norm_output             % Normalised output
        
        centre                  % Centre matrix
        delta                   % Delta matrix
        fval_min                % GP objective function value matrix
        opt_min                 % Optimal point matrix
        rho                     % Rho matrix
    end
    properties (Constant)
        eta_low = 0.1           % Rho constant
        eta_high = 0.9          % Rho constant
        delta_reduction = 0.5   % Reduction in delta when Rho < eta_low
        delta_expansion = 1.5   % Expansion in delta when Rho > eta_high
    end
    methods
        function obj = GPCreator(system, training_input, training_output)
            obj.system = system;
            obj.objective = system.objective;
            obj.linear_con_A = system.linear_con_A;
            obj.linear_con_b = system.linear_con_b;
            obj.lineq_con_A = system.lineq_con_A;
            obj.lineq_con_b = system.lineq_con_b;
            obj.lb = system.lb;
            obj.ub = system.ub;
            obj.nonlin_con = system.nonlin_con;
            obj.options = system.options;
            obj.output_fields = system.output_fields;
            
            obj.training_input = training_input;
            obj.training_output = training_output;
            obj.delta = system.delta_mat;
            obj.centre = obj.training_input(1, :);
            obj.fval_min(1) = obj.training_output(1, :);
            obj.opt_min(1, :) = obj.training_input(1, :);
            obj.rho = zeros();
            
            obj.update_model();
        end
        
        function obj = update_model(obj)
            % Normalise data
            obj.mean_input = mean(obj.training_input);
            obj.std_input = std(obj.training_input);
            obj.norm_input = normalize(obj.training_input);
            obj.mean_output = mean(obj.training_output);
            obj.std_output = std(obj.training_output);
            obj.norm_output = normalize(obj.training_output);
            
            % Initialise GP model
            obj.model = fitrgp(obj.norm_input, obj.norm_output);
        end
        
        function optimise(obj, iter)
            % Pre-allocation
            columns = length(obj.output_fields);
            rows = size(obj.centre, 1);
            obj.centre = [obj.centre; zeros(iter, columns)];
            obj.delta = [obj.delta; zeros(iter, columns)];
            obj.fval_min = [obj.fval_min; zeros(iter, 1)];
            obj.opt_min = [obj.opt_min; zeros(iter, columns)];
            obj.rho = [obj.rho; zeros(iter, 1)];
            
            % Initialise
            [~, output_true_last] = obj.system.get_output(obj.centre(rows, :));
            pointer = Pointer(obj.centre(rows, :), obj.delta(rows, :));
            
            for i = rows+1:rows+iter
                % Update nonlinear constraints with current iteration data
                nonlincon = @(x) obj.nonlin_con( ...
                    x, obj.centre(i-1, :), obj.delta(i-1, :), obj.model, obj.mean_input, ...
                    obj.std_input, obj.mean_output, obj.std_output ...
                );
                
                % Optimise from various starting points
                starting_pts = pointer.random_sampling();
                dim = size(starting_pts);
                opt_points = zeros(dim(1), dim(2));
                fvals = zeros(dim(1), 1);
                for each = 1:dim(1)
                    [opt, fval] = fmincon( ...
                        obj.objective, starting_pts(each, :), obj.linear_con_A, ...
                        obj.linear_con_b, obj.lineq_con_A, obj.lineq_con_b, ...
                        obj.lb, obj.ub, nonlincon, obj.options);
                    opt_points(each, :) = opt;
                    fvals(each) = fval;
                end
                
                % Get lowest (local) optima out of the starting points
                [obj.fval_min(i), idx] = min(fvals);
                obj.opt_min(i, :) = opt_points(idx, :);
                
                % Rho calculation
                [~, output_true_curr] = obj.system.get_output(obj.opt_min(i, :));
                output_pred_last = predict( ...
                    obj.model, (obj.centre(i-1, :) - obj.mean_input) ./ obj.std_input);
                output_pred_curr = predict( ...
                    obj.model, (obj.opt_min(i, :) - obj.mean_input) ./ obj.std_input);
                obj.rho(i) = ( ...
                    (obj.objective(output_true_curr) - obj.objective(output_true_last)) / ...
                    (obj.objective(output_pred_curr) - obj.objective(output_pred_last)) ...
                );
            
                switch obj.rho(i)
                    case obj.rho(i) < obj.eta_low
                        obj.centre(i, :) = obj.centre(i - 1, :);
                        obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_reduction;
                    case obj.rho(i) < obj.eta_high
                        obj.centre(i, :) = obj.opt_min(i, :);
                        obj.delta(i, :) = obj.delta(i - 1, :);
                    otherwise
                        obj.centre(i, :) = obj.opt_min(i, :);
                        obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_expansion;
                end
                    
                pointer.update(obj.centre(i, :), obj.delta(i, :));
                obj.training_input = [obj.training_input; obj.opt_min(i, :)];
                obj.training_output = [obj.training_output; output_true_curr];
                obj.update_model();
            end
        end
        
        function plot(obj)
            
        end
    end
end