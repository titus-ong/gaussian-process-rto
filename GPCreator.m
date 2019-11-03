classdef GPCreator  < handle
    properties
        system                  % Either HYSYSFile or system object
        objective_model         % Objective function for model
        objective_true          % Objective function for system
        obj_fn                  % Objective function with model inputs
        linear_con_A            % System linear inequality constraints LHS (Ax = b)
        linear_con_b            % System linear inequality constraints RHS
        lineq_con_A             % System linear equality constraints LHS
        lineq_con_b             % System linear equality constraints RHS
        lb                      % System lower bounds
        ub                      % System upper bounds
        nonlin_con              % System nonlinear constraints
        options                 % System options for GP
        system_violation        % Test for true violation of system constraints
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
            obj.objective_model = system.objective_model;
            obj.objective_true = system.objective_true;
            obj.linear_con_A = system.linear_con_A;
            obj.linear_con_b = system.linear_con_b;
            obj.lineq_con_A = system.lineq_con_A;
            obj.lineq_con_b = system.lineq_con_b;
            obj.lb = system.lb;
            obj.ub = system.ub;
            obj.nonlin_con = system.nonlin_con;
            obj.options = system.options;
            obj.system_violation = system.system_violation;
            obj.output_fields = system.output_fields;
            
            obj.training_input = training_input;
            obj.training_output = training_output;
            obj.delta = system.delta_mat;
            obj.centre = obj.training_input(1, :);
            obj.fval_min(1) = 0;
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
            for i = 1:size(obj.training_output, 2)
                obj.model.(obj.output_fields{i}) = fitrgp(obj.training_input, obj.training_output(:, i));
            end
            
            % Update objective function with current iteration data
            obj.obj_fn = @(x) obj.objective_model( ...
                x, obj.model, obj.mean_input, obj.std_input, ...
                obj.mean_output, obj.std_output ...
            );
        end
        
        function optimise(obj, iter)
            % Pre-allocation
            cols = size(obj.centre, 2);
            rows = size(obj.centre, 1);
            obj.centre = [obj.centre; zeros(iter, cols)];
            obj.delta = [obj.delta; zeros(iter, cols)];
            obj.fval_min = [obj.fval_min; zeros(iter, 1)];
            obj.opt_min = [obj.opt_min; zeros(iter, cols)];
            obj.rho = [obj.rho; zeros(iter, 1)];
            
            % Initialise
            pointer = Pointer(obj.centre(rows, :), obj.delta(rows, :));
            
            for i = rows+1:rows+iter
                % Update nonlinear constraints with current iteration data
                nonlincon = @(x) obj.nonlin_con( ...
                    x, obj.centre(i-1, :), obj.delta(i-1, :), obj.model, obj.mean_input, ...
                    obj.std_input, obj.mean_output, obj.std_output ...
                );
                
                % Optimise from various starting points
                starting_pts = pointer.random_sampling(10);
                dim = size(starting_pts);
                opt_points = zeros(dim(1), dim(2));
                fvals = zeros(dim(1), 1);
                for each = 1:dim(1)
                    [opt, fval] = fmincon( ...
                        obj.obj_fn, starting_pts(each, :), obj.linear_con_A, ...
                        obj.linear_con_b, obj.lineq_con_A, obj.lineq_con_b, ...
                        obj.lb, obj.ub, nonlincon, obj.options);
                    opt_points(each, :) = opt;
                    fvals(each) = fval;
                end
                
                % Get lowest (local) optima out of the starting points
                [obj.fval_min(i), idx] = min(fvals);
                obj.opt_min(i, :) = opt_points(idx, :);
                
                % True values
                % Use objective method for true because we don't want to have model inputs
                true_curr = obj.objective_true(obj.opt_min(i, :));
                true_last = obj.objective_true(obj.centre(i-1, :));
                
                % Train new GP
                [true_output, ~] = obj.system.get_output(obj.opt_min(i, :));
                obj.training_input = [obj.training_input; obj.opt_min(i, :)];
                obj.training_output = [obj.training_output; true_output];
                obj.update_model();
                
                % Predicted values
                predicted_curr = obj.obj_fn(obj.opt_min(i, :));
                predicted_last = obj.obj_fn(obj.centre(i-1, :));
                
                % Rho calculation
                obj.rho(i) = ( ...
                    (true_curr - true_last) / (predicted_curr - predicted_last) ...
                );
            
                if obj.system_violation(true_output)
                    % Current point violates system constraints
                    obj.centre(i, :) = obj.centre(i - 1, :);
                    obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_reduction;
                    obj.rho(i) = NaN;
                elseif obj.rho(i) < obj.eta_low
                    obj.centre(i, :) = obj.centre(i - 1, :);
                    obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_reduction;
                elseif obj.rho(i) < obj.eta_high
                    obj.centre(i, :) = obj.opt_min(i, :);
                    obj.delta(i, :) = obj.delta(i - 1, :);
                elseif obj.rho(i) >= obj.eta_high
                    obj.centre(i, :) = obj.opt_min(i, :);
                    obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_expansion;
                else  % Rho calculated is NaN
                    obj.centre(i, :) = obj.centre(i - 1, :);
                    obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_reduction;
                end
                    
                % Update values
                pointer.update(obj.centre(i, :), obj.delta(i, :));
            end
        end
        
        function plot(obj)
            % Plot centre moving against objective function
            % Plot delta around centres
            if size(obj.centre, 2)==2
                obj.plot2d_indiv();
            end
        end
        
        function plot2d_indiv(obj)
            f = figure;
            hold on;
            
            idx = size(obj.centre, 1);
            points = zeros(idx, 2);
            deltas = zeros(idx, 2);
            syms x y a b h k
            for i = 1:size(obj.centre, 1)
                points(i, :) = obj.centre(i, :);
                a = obj.delta(i, 1);
                b = obj.delta(i, 2);
                h = points(i, 1);
                k = points(i, 2);
                ellipse = (((x-h)^2)/(a^2))+(((y-k)^2)/(b^2))==1;
                fimplicit(ellipse, [obj.lb(1) obj.ub(1) obj.lb(2) obj.ub(2)]);
                hold on;
            end
            plot(points(:, 1), points(:, 2), '-*');
%             x1 = linspace(obj.lb(2), obj.ub(2), 10);
%             y1 = linspace(obj.lb(1), obj.ub(1), 10);
%             fvals = zeros(10, 10);
%             co2 = zeros(10, 10);
%             for i = 1:10
%                 for j = 1:10
%                     fvals(j, i) = obj.obj_fn([y1(i), x1(j)]);
%                     co2(j, i) = predict(obj.model.clean_gas_co2, [y1(i), x1(j)]);
%                 end
%             end
%             figure;
%             hold on;
%             surf(x1, y1, fvals);
%             surf(x1, y1, co2);
        end
    end
end