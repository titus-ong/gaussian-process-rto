classdef PeaksFunction  < matlab.mixin.Copyable
    properties (Constant)
        feasible_point = struct( ...
            "x", -0.3, ...
            "y", 1 ...
            );
        delta_point = struct( ...
            "x", 0.6, ...
            "y", 0.6 ...
            );
        output_fields = {'z'};
        constraints_ineq = [];
        constraints_eq = [];
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [-3, -3];                             % Lower bounds
        ub = [3, 3];                               % Upper bounds
        options = optimset('disp','off');          % Options for GP
        
        min_TR = 0.01                              % Minimum trust region as percentage of original delta
        max_TR = 10                                % Maximum trust region as percentage of original delta
        eta_low = 0.1                              % Rho constant
        eta_high = 0.9                             % Rho constant
        delta_reduction = 0.5                      % Reduction in delta when Rho < eta_low
        delta_expansion = 1.5                      % Expansion in delta when Rho > eta_high
        forgetting_factor = 1.5                    % Allowance for inaccuracies in GP due to outdated data        
    end
    properties
        feasible_point_mat                         % Matrix form
        delta_mat                                  % Matrix form
        decay                                      % Decay boolean
        forget                                     % Forgetting factor
        time                                       % Pseudo time for decay
        op_region_script                           % Functionf or plotting operation region
    end
    methods
        function obj = PeaksFunction()
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
            obj.delta_mat = cell2mat(struct2cell(obj.delta_point))';
            obj.decay = false;
            obj.forget = false;
            obj.time = 0;
            obj.op_region_script = @op_region_plot_peaks;
        end
        
        function objective = objective_value(obj, inputs, outputs)
            % Calculate objective function
            objective = outputs;

        end

        function obj = time_increment(obj)
            if obj.decay
                obj.time = obj.time + 1;
            end
        end
        
        function [c,ceq] = nonlin_con(obj, x, par)  % centre, delta, model, mean_x, std_x, mean_y, std_y)
            % Nonlinear inequality (c<=0) and equality (ceq=0) constraints on x
            c = zeros(1 + length(obj.constraints_ineq), 1);
            
            % Inequality constraints
            for i = 1:length(obj.constraints_ineq)
                c(i) = predict(par.model.(obj.constraints_ineq{i}), ( ...
                    x - par.values_adj.input.mean) ./ par.values_adj.input.std ...
                    ) ...
                    * par.values_adj.(obj.constraints_ineq{i}).std ...
                    + par.values_adj.(obj.constraints_ineq{i}).mean;
            end
            
            % Point is within delta - must be the final equation!
            c(end) = sum((x - par.centre) .^ 2 ./ par.delta .^ 2) - 1;
            
            ceq = zeros(length(obj.constraints_eq), 1);
            % Equality constraints
            for i = 1:length(obj.constraints_eq)
                ceq(i) = predict(par.model.(obj.constraints_eq{i}), ( ...
                    x - par.values_adj.input.mean) ./ par.values_adj.input.std ...
                    ) ...
                    * par.values_adj.(obj.constraints_eq{i}).std ...
                    + par.values_adj.(obj.constraints_eq{i}).mean;
            end
        end
        
        function par = create_par(~, GPobj, idx)
            % Create parameters for fmincon (used in nonlin_fn)
            par.model = GPobj.model(end);
            par.values_adj = GPobj.values_adj;
            par.centre = GPobj.centre(idx, :);
            par.delta = GPobj.delta(idx, :);
%             par.delta_norm = cell2mat(struct2cell(obj.delta_norm))';
        end
        
        function bool = system_violated(obj, con_ineq, con_eq)
            % Test if output violates system constraints
            % Return true if violated, false if not
            tol = 1e-4;
            bool = false;
            
            % Check inequality constraints
            for i = 1:length(con_ineq)
                if con_ineq(i) > 0
                    bool = true;
                end
            end
            
            % Check equality constraints
            for i = 1:length(con_eq)
                if abs(con_eq(i)) > tol
                    bool = true;
                end
            end
            
            obj.time_increment();
            
        end
        
        
        function [objective, con_ineq, con_eq] = get_output(obj, inputs)
            % Get objective function and constraint values using given inputs
            n = size(inputs, 1);
            m = length(obj.constraints_ineq);
            k = length(obj.constraints_eq);
            objective = zeros(n, 1);
            con_ineq = zeros(n, m);
            con_eq = zeros(n, k);
            
            % Solve reactor
            for point = 1:n
                curr_inputs = inputs(point, :);
                z = peaks(curr_inputs(1), curr_inputs(2)+min(max((obj.time-10)/5,0),1));
                objective(point) = obj.objective_value(curr_inputs, z);
                for idx = 1:m
                    val = obj.(obj.constraints_ineq{idx})(curr_inputs, z);
                    con_ineq(point, idx) = val;
                end
                for idx = 1:k
                    val = obj.(obj.constraints_eq{idx})(curr_inputs, z);
                    con_eq(point, idx) = val;
                end
            end
        end
    end
end