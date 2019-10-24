classdef GPCreator
    properties
        system          % Either HYSYSFile or system object
        delta           % System delta
        objective       % System objective function
        linear_con_A    % System linear inequality constraints LHS (Ax = b)
        linear_con_b    % System linear inequality constraints RHS
        lineq_con_A     % System linear equality constraints LHS
        lineq_con_b     % System linear equality constraints RHS
        lb              % System lower bounds
        ub              % System upper bounds
        nonlin_con      % System nonlinear constraints
        options         % System options for GP
        model           % GP model
        
        init_input      % Initial training input data
        mean_input      % Mean of input
        std_input       % Standard deviation of input
        norm_input      % Normalised input
        
        init_output     % Initial training output data
        mean_output     % Mean of output
        std_output      % Standard deviation of output
        norm_output     % Normalised output
    end
    properties (Constant)
        eta_low = 0.1   % Rho constant
        eta_high = 0.9  % Rho constant
    end
    methods
        function obj = GPCreator(system, training_input, training_output)
            obj.system = system;
            obj.delta = system.delta_mat;
            obj.objective = system.objective;
            obj.linear_con_A = system.linear_con_A;
            obj.linear_con_b = system.linear_con_b;
            obj.lineq_con_A = system.lineq_con_A;
            obj.lineq_con_b = system.lineq_con_b;
            obj.lb = system.lb;
            obj.ub = system.ub;
            obj.nonlin_con = system.nonlin_con;
            obj.options = system.options;
            
            obj.init_input = training_input;
            obj.init_output = training_output;
            
            % Normalise data
            obj.mean_input = mean(obj.init_input);
            obj.std_input = std(obj.init_input);
            obj.norm_input = normalize(obj.init_input);
            obj.mean_output = mean(obj.init_output);
            obj.std_output = std(obj.init_output);
            obj.norm_output = normalize(obj.init_output);
            
            % Initialise GP model
            obj.model = fitrgp(obj.norm_input, obj.norm_output);
        end
        
        function optimise(obj, iter)
            centre = obj.init_input(1, :);
            delta_ = obj.delta;
            pointer = Pointer(centre, delta_);
            
            for i = 1:iter
                % Update nonlinear constraints with current iteration data
                nonlincon = @(x) obj.nonlin_con( ...
                    x, centre, delta_, obj.model, obj.mean_input, ...
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
                    
                
                
                pointer.update(centre, delta_);
            end
        end
    end
end