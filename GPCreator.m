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
        dir_vec                 % Directional vector of current gradient
        excited                 % Logical array of whether point is an excited point
        
        forget                  % Forgetting factor boolean
        excite                  % Excite boolean
        
        min_TR                  % Minimum trust region as percentage of original delta
        max_TR                  % Maximum trust region as percentage of original delta
        eta_low                 % Rho constant
        eta_high                % Rho constant
        delta_reduction         % Reduction in delta when Rho < eta_low
        delta_expansion         % Expansion in delta when Rho > eta_high
        forgetting_factor       % Allowance for inaccuracies in GP due to outdated data
        constraint_tol          % Tolerance when system constraint is violated
        excite_tol              % Tolerance of points being aligned for excitation
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
            obj.constraint_tol = system.constraint_tol;
            obj.excite_tol = system.excite_tol;
            
            obj.training_input = training_input;
            obj.training_starter = training_input;
            obj.training_output = training_output;
            obj.delta = system.delta_mat;
            obj.centre = obj.training_input(1, :);
            obj.fval_min(1) = 0;
            obj.opt_min(1, :) = obj.training_input(1, :);
            obj.rho = zeros();
            obj.excited = zeros();
            
            obj.model = struct();
            obj.update_model();
            
            if isprop(system, "forget")
                obj.forget = system.forget;
            else
                obj.forget = false;
            end
            
            obj.excite = true;
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
        
        function bool = system_violated(obj, con_ineq, con_eq)
            % Test if output violates system constraints
            % Return true if violated, false if not
            bool = false;
            
            % Check inequality constraints
            for i = 1:length(con_ineq)
                if con_ineq(i) > obj.constraint_tol
                    bool = true;
                end
            end
            
            % Check equality constraints
            for i = 1:length(con_eq)
                if abs(con_eq(i)) > obj.constraint_tol
                    bool = true;
                end
            end
        end 
        
        function objective = obj_fn(obj, x)
            objective = predict( ...
                obj.model(end).objective, ...
                    (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                    ) * obj.values_adj(end).objective.std + obj.values_adj(end).objective.mean;
        end
        
        function [con_ineq, con_eq] = model_con(obj, x)
            con_ineq = zeros(size(obj.training_output.con_ineq, 2), 1);
            for i = 1:size(con_ineq, 1)
                con_ineq(i) = predict( ...
                    obj.model(end).(obj.con_ineq{i}), ...
                        (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                        ) * obj.values_adj(end).(obj.con_ineq{i}).std + obj.values_adj(end).(obj.con_ineq{i}).mean;
            end
            
            con_eq = zeros(size(obj.training_output.con_eq, 2), 1);
            for i = 1:size(con_eq, 1)
                con_eq(i) = predict( ...
                    obj.model(end).(obj.con_eq{i}), ...
                        (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                        ) * obj.values_adj(end).(obj.con_eq{i}).std + obj.values_adj(end).(obj.con_eq{i}).mean;
            end
        end
        
        function bool = should_excite(obj, idx)
            % Check if past two gradients have been close to parallel
            % Excite if parallel within tolerance
            if idx < 3
                bool = false;
            elseif ~obj.excite
                bool = false;
            else
                grad_prev = obj.centre(idx - 1, :) - obj.centre(idx - 2, :);
                unit_prev = grad_prev ./ norm(grad_prev);
                grad_curr = obj.centre(idx, :) - obj.centre(idx - 1, :);
                unit_curr = grad_curr ./ norm(grad_curr);
                if length(unit_prev) < 3
                    vec_len = unit_prev(1)*unit_curr(2) - unit_prev(2)*unit_curr(1);
                else
                    x_pdt = cross(unit_prev, unit_curr);
                    vec_len = norm(x_pdt);
                end
                if abs(vec_len) < obj.excite_tol
                    bool = true;
                elseif sum(isnan(unit_prev)) && sum(isnan(unit_curr))
                    bool = true;
                else
                    bool = false;
                end
                obj.dir_vec = unit_curr;
            end
        end
        
        function excited = get_excited_points(obj, idx, direction_vec)
            % Get excited points given previous direction vector and
            % current centre
            if sum(isnan(direction_vec))
                excited = obj.stationary(idx);
            else
                excited = obj.straight_line(idx, direction_vec);
            end
        end
        
        function excited = stationary(obj, idx)
            % Find a random point within the trust region from the centre
            pointer = Pointer(obj.centre(idx, :), obj.delta(idx, :), obj.lb, obj.ub);
            point_list = pointer.random_sampling(1);
            excited = point_list(2, :);  % point_list's first point is the centre
        end
        
        function excited = straight_line(obj, idx, direction_vec)
            % Find n points that are orthogonal to direction vector and are
            % within delta, where n is the number of dimensions
            
            vec_dim = length(direction_vec);
            excited = zeros(vec_dim, size(obj.centre, 2));
            for element_no = 1:vec_dim*2  % Include negative direction vectors
                % Get orthogonal vector of e.g. [1, 1, x]
                curr_element = fix((element_no+1)/2);  % changes every 2 iterations
                vec_magnitude = 0;
                for i = 1:vec_dim
                    if i == curr_element
                        continue
                    end
                    vec_magnitude = vec_magnitude + direction_vec(i);
                end
                last_element = -vec_magnitude / direction_vec(curr_element);
                vec_ortho = ones(1, vec_dim);
                vec_ortho(curr_element) = last_element;

                % Find magnitude of direction vector that will produce the
                % furthest excited point
                squared_x = 1 / sum(vec_ortho .^ 2 ./ obj.delta(idx, :) .^ 2);
                x = sqrt(squared_x) * (-1)^element_no;  % Changes between positive and negative

                excited(element_no, :) = obj.centre(idx, :) + x * vec_ortho;
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
            obj.excited = [obj.excited; zeros(iter, 1)];
            
            % Initialise
            % pointer = Pointer(obj.centre(rows, :), obj.delta(rows, :));
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
                
                % Excite or optimise
                if obj.should_excite(i-1) && ~obj.excited(i-1)  % Avoid excitation if previous iter was excited
                    points = obj.get_excited_points(i-1, obj.dir_vec);
                    feasible = zeros(0, 2);
                    for j = 1:size(points, 1)
                        [ineq, eq] = obj.model_con(points(j, :));
                        if obj.system_violated(ineq, eq)
                            % Violate system constraints
                            continue
                        elseif sum(points(j, :) > obj.ub) || sum(points(j, :) < obj.lb)
                            % Violate lower/upper bounds
                            continue
                        else
                            % Feasible point
                            rw = size(feasible, 1) + 1;
                            feasible(rw, 1) = j;
                            feasible(rw, 2) = func_obj(points(j, :));
                        end
                    end
                    if ~isempty(feasible)
                        % Find optimal feasible point
                        [rw, ~] = find(feasible(:, 2)==max(feasible(:, 2)));
                        rw = rw(1);  % If there are duplicate optimas
                        obj.fval_min(i) = feasible(rw, 2);
                        obj.opt_min(i, :) = points(feasible(rw, 1), :);
                        obj.excited(i) = true;
                    end
                end
                if ~logical(obj.fval_min(i))
                    % Optimise using fmincon
                    [opt, fval] = fmincon( ...
                        func_obj, obj.centre(i-1, :), obj.linear_con_A, ...
                        obj.linear_con_b, obj.lineq_con_A, obj.lineq_con_b, ...
                        obj.lb, obj.ub, nonlin_con, obj.options, par ...
                        );
                    obj.fval_min(i) = fval;
                    obj.opt_min(i, :) = opt;
                    obj.excited(i) = false;
                end
                
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
                new_is_worse = (true_curr - true_last) > 0;
            
                if obj.system_violated(ineq, eq)
                    % Current point violates system constraints
                    obj.centre(i, :) = obj.centre(i - 1, :);
                    obj.delta(i, :) = obj.delta(i - 1, :) * obj.delta_reduction;
                    obj.rho(i) = NaN;
                elseif new_is_worse && obj.excited(i)
                    obj.centre(i, :) = obj.centre(i - 1, :);
                    obj.delta(i, :) = obj.delta(i - 1, :);                    
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
                % pointer.update(obj.centre(i, :), obj.delta(i, :));
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
                fimplicit(ellipse, [obj.lb(1) obj.ub(1) obj.lb(2) obj.ub(2)], '--b');
                hold on;
            end
            centres = plot(points(:, 1), points(:, 2), '-b*');
            training = scatter(obj.training_starter(:,1), obj.training_starter(:,2), '+m');
            for i = 1:length(obj.excited)
                if ~obj.excited(i)
                    continue
                end
                xval = [points(i-1, 1), obj.opt_min(i, 1), points(i+1, 1)];
                yval = [points(i-1, 2), obj.opt_min(i, 2), points(i+1, 2)];
                optimas = plot(xval, yval, '-ro');
            end
            legend([centres training optimas], {'Centres', 'Training inputs', 'Excited points'});
        end
    end
end