classdef GPCreator  < matlab.mixin.Copyable
    properties
        system                  % Either HYSYSFile or system object
        linear_con_A            % System linear inequality constraints LHS (Ax = b)
        linear_con_b            % System linear inequality constraints RHS
        lineq_con_A             % System linear equality constraints LHS
        lineq_con_b             % System linear equality constraints RHS
        con_ineq                % System nonlinear inequality constraints
        con_eq                  % System nonlinear equality constraints
        lb                      % System lower bounds
        ub                      % System upper bounds
        options                 % System options for GP
        model                   % GP model
        
        training_input          % Training input data
        training_output         % Training output data
        values_adj              % Mean, std and normalised values
        training_starter        % Starter training input points
        
        centre                  % Centre matrix
        delta                   % Delta matrix
        fval_min                % GP objective function value matrix
        opt_min                 % Optimal point matrix
        rho                     % Rho matrix
        
        forget                  % Forgetting factor boolean
        
        min_TR                  % Minimum trust region as percentage of original delta
        max_TR                  % Maximum trust region as percentage of original delta
        eta_low                 % Rho constant
        eta_high                % Rho constant
        delta_reduction         % Reduction in delta when Rho < eta_low
        delta_expansion         % Expansion in delta when Rho > eta_high
        forgetting_factor       % Allowance for inaccuracies in GP due to outdated data
    end
    methods
        function obj = GPCreator(system, training_input, training_output)
            obj.system = system;
            obj.linear_con_A = system.linear_con_A;
            obj.linear_con_b = system.linear_con_b;
            obj.lineq_con_A = system.lineq_con_A;
            obj.lineq_con_b = system.lineq_con_b;
            obj.lb = system.lb;
            obj.ub = system.ub;
            obj.options = system.options;
            obj.con_ineq = system.constraints_ineq;
            obj.con_eq = system.constraints_eq;
            
            obj.min_TR = system.min_TR;
            obj.max_TR = system.max_TR;
            obj.eta_low = system.eta_low;
            obj.eta_high = system.eta_high;
            obj.delta_reduction = system.delta_reduction;
            obj.delta_expansion = system.delta_expansion;
            obj.forgetting_factor = system.forgetting_factor;
            
            obj.training_input = training_input;
            obj.training_starter = training_input;
            obj.training_output = training_output;
            obj.delta = system.delta_mat;
            obj.centre = obj.training_input(1, :);
            obj.fval_min(1) = 0;
            obj.opt_min(1, :) = obj.training_input(1, :);
            obj.rho = zeros();
            
            obj.model = struct();
            obj.update_model();
            
            if isprop(system, "forget")
                obj.forget = system.forget;
            else
                obj.forget = false;
            end
        end
        
        function obj = update_model(obj)  
            % Initialise or update model
            idx = size(obj.model, 2) + 1;
            if isempty(fieldnames(obj.model))  % model is empty
                idx = 1;  % So model fills up from first row
            end
            obj.values_adj(idx).input.mean = mean(obj.training_input);
            obj.values_adj(idx).input.std = std(obj.training_input);
            obj.values_adj(idx).input.norm = normalize(obj.training_input);
            
            % Objective model
            obj.values_adj(idx).objective.mean = mean(obj.training_output.objective);
            obj.values_adj(idx).objective.std = std(obj.training_output.objective);
            obj.values_adj(idx).objective.norm = normalize(obj.training_output.objective);
            obj.model(idx).objective = fitrgp( ...
                obj.values_adj(idx).input.norm, obj.values_adj(idx).objective.norm ...
                );
            
            % Inequality constraints models
            for i = 1:size(obj.training_output.con_ineq, 2)
                obj.values_adj(idx).(obj.con_ineq{i}).mean = mean(obj.training_output.con_ineq(:, i));
                obj.values_adj(idx).(obj.con_ineq{i}).std = std(obj.training_output.con_ineq(:, i));
                obj.values_adj(idx).(obj.con_ineq{i}).norm = normalize(obj.training_output.con_ineq(:, i));
                obj.model(idx).(obj.con_ineq{i}) = fitrgp( ...
                    obj.values_adj(idx).input.norm, obj.values_adj(idx).(obj.con_ineq{i}).norm ...
                    );
            end
            
            % Equality constraints models
            for i = 1:size(obj.training_output.con_eq, 2)
                obj.values_adj(idx).(obj.con_eq{i}).mean = mean(obj.training_output.con_eq(:, i));
                obj.values_adj(idx).(obj.con_eq{i}).std = std(obj.training_output.con_eq(:, i));
                obj.values_adj(idx).(obj.con_eq{i}).norm = normalize(obj.training_output.con_eq(:, i));
                obj.model(idx).(obj.con_eq{i}) = fitrgp( ...
                    obj.values_adj(idx).input.norm, obj.values_adj(idx).(obj.con_eq{i}).norm ...
                    );
            end
        end
        
        function forget_outdated(obj, true_objective)
            % Check if old training inputs are irrelevant and remove them
            % objective_long: with all training inputs
            % objective_short: less the first (oldest) training input
            inaccurate = true;
            while size(obj.training_input, 1) > 10 && inaccurate
                % Initialise short GP model
                mean_val = mean(obj.training_output.objective(2:end));
                std_val = std(obj.training_output.objective(2:end));
                norm_val = normalize(obj.training_output.objective(2:end));
                mean_input = mean(obj.training_input(2:end, :));
                std_input = std(obj.training_input(2:end, :));
                norm_input = normalize(obj.training_input(2:end, :));
                model_short = fitrgp(norm_input, norm_val);
            
                % Get predicted objective
                objective_short = predict( ...
                    model_short, ((obj.opt_min(end, :)) - mean_input) ./ std_input ...
                    ) * std_val + mean_val;
                objective_long = predict( ...
                    obj.model(end).objective, ( ...
                        (obj.opt_min(end, :)) - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                        ) * obj.values_adj(end).objective.std + obj.values_adj(end).objective.mean;
            
                ratio = abs((objective_long - true_objective) / (objective_short - true_objective));
                
                if ratio > obj.forgetting_factor
                    obj.training_input = obj.training_input(2:end, :);
                    obj.training_output.objective = obj.training_output.objective(2:end, :);
                    obj.training_output.con_eq = obj.training_output.con_eq(2:end, :);
                    obj.training_output.con_ineq = obj.training_output.con_ineq(2:end, :);
                    obj.update_model();
                else
                    inaccurate = false;
                end
            end
        end        
        
        function objective = obj_fn(obj, x)
            objective = predict( ...
                obj.model(end).objective, ...
                    (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                    ) * obj.values_adj(end).objective.std + obj.values_adj(end).objective.mean;
        end
        
        function bool = should_excite(obj, idx)
            tol = 1e-2;
            grad_prev = obj.centre(idx - 1, :) - obj.centre(idx - 2, :);
            unit_prev = grad_prev ./ norm(grad_prev);
            grad_curr = obj.centre(idx, :) - obj.centre(idx - 1, :);
            unit_curr = grad_curr ./ norm(grad_curr);
            x_pdt = cross(unit_prev, unit_curr);
            vec_len = norm(x_pdt);
            if abs(vec_len) < tol
                bool = true;
            else
                bool = false;
            end
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
            true_last = obj.system.get_output(obj.centre(rows, :));
            
            for i = rows+1:rows+iter
                % Update parameters for constraints with current iteration data
                par = obj.system.create_par(obj, i-1);
                func_obj = @(x, ~) obj.obj_fn(x);
                nonlin_con = @(x, par) obj.system.nonlin_con(x, par);
                
%                 % Optimise from various starting points
%                 starting_pts = pointer.random_sampling(1);
%                 dim = size(starting_pts);
%                 opt_points = zeros(dim(1), dim(2));
%                 fvals = zeros(dim(1), 1);
%                 for each = 1:dim(1)
%                     [opt, fval] = fmincon( ...
%                         func_obj, starting_pts(each, :), obj.linear_con_A, ...
%                         obj.linear_con_b, obj.lineq_con_A, obj.lineq_con_b, ...
%                         obj.lb, obj.ub, nonlin_con, obj.options, par ...
%                         );
%                     opt_points(each, :) = opt;
%                     fvals(each) = fval;
%                 end
%                 
%                 % Get lowest (local) optima out of the starting points
%                 [obj.fval_min(i), idx] = min(fvals);
%                 obj.opt_min(i, :) = opt_points(idx, :);

                % Optimise using fmincon
                [opt, fval] = fmincon( ...
                    func_obj, obj.centre(i-1, :), obj.linear_con_A, ...
                    obj.linear_con_b, obj.lineq_con_A, obj.lineq_con_b, ...
                    obj.lb, obj.ub, nonlin_con, obj.options, par ...
                    );
                obj.fval_min(i) = fval;
                obj.opt_min(i, :) = opt;
                
                % Move to new point and get true values
                [true_curr, ineq, eq] = obj.system.get_output(obj.opt_min(i, :));
                
                % Train new GP
                obj.training_input = [obj.training_input; obj.opt_min(i, :)];
                obj.training_output.objective = [obj.training_output.objective; true_curr];
                obj.training_output.con_ineq = [obj.training_output.con_ineq; ineq];
                obj.training_output.con_eq = [obj.training_output.con_eq; eq];
                obj.update_model();
                
                % Account for outdated data
                if obj.forget
                    obj.forget_outdated(true_curr);
                end
                
                % Predicted values
                predicted_curr = obj.obj_fn(obj.opt_min(i, :));
                predicted_last = obj.obj_fn(obj.centre(i-1, :));
                
                % Rho calculation
                obj.rho(i) = ( ...
                    (true_curr - true_last) / (predicted_curr - predicted_last) ...
                );
            
                if obj.system.system_violated(ineq, eq)
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
                    obj.rho(i) = Inf;
                end
                
                % Check minimum and maximum trust region size
                ratio = obj.delta(i, 1) / obj.delta(rows, 1);
                if ratio > obj.max_TR || ratio < obj.min_TR
                    obj.delta(i, :) = obj.delta(i - 1, :);
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
            obj.system.op_region_script();
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
            scatter(obj.training_starter(:,1), obj.training_starter(:,2));
        end
    end
end