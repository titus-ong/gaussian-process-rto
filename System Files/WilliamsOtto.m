classdef WilliamsOtto  < matlab.mixin.Copyable
    properties (Constant)
        feasible_point = struct( ...
            "flowrate_b", 6.9, ...
            "t_reactor", 83 ...
            );
%         delta = 0.15;
%         delta_lb = 0.005;
        delta_point = struct( ...
            "flowrate_b", 0.4, ...
            "t_reactor", 4 ...
            );
%         delta_norm = struct( ...
%             "flowrate_b", 9, ...
%             "t_reactor", 900 ...
%             );
        constraints_ineq = { ...
            'x_a_con', ...
            'x_g_con' ...
            };
        constraints_eq = [];
        constants = struct( ...
            "flowrate_a", 1.8275, ...              % Inlet flowrate of A
            "total_molar_flowrate", 2105.2 ...     % Total molar flowrate
            );
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [4, 70];                              % Lower bounds
        ub = [7, 100];                             % Upper bounds
        options = optimset('disp','off');          % Options for GP

        min_TR = 0.01                              % Minimum trust region as percentage of original delta
        max_TR = 10                                % Maximum trust region as percentage of original delta
        eta_low = 0.1                              % Rho constant
        eta_high = 0.9                             % Rho constant
        delta_reduction = 0.8                      % Reduction in delta when Rho < eta_low
        delta_expansion = 1.2                      % Expansion in delta when Rho > eta_high
        forgetting_factor = 1.5                    % Allowance for inaccuracies in GP due to outdated data
        constraint_tol = 1e-3                      % Tolerance when system constraint is violated
        align_tol = 1e-3                           % Tolerance of points being aligned for excitation
        region_tol = 5e-4                          % Tolerance (fraction of max TR) of points in same region for excitation
    end
    properties
        init_var = struct( ...
            "x_a", 0.2, ...                        % Outlet concentration of A
            "x_b", 0.2, ...                        % Outlet concentration of B
            "x_c", 0.2, ...                        % Outlet concentration of C
            "x_e", 0.2, ...                        % Outlet concentration of E
            "x_g", 0.2, ...                        % Outlet concentration of G
            "x_p", 0.2 ...                         % Outlet concentration of P
            );
        
        feasible_point_mat                         % Matrix form
        delta_mat                                  % Matrix form
        decay                                      % Decay boolean
        time                                       % Pseudo time for decay
        op_region_script                           % Function for plotting operating region
        forget                                     % Switch for forgetting factor
        x_a_val                                    % x_a constraint value (for moving constraint)
    end
    methods
        function obj = WilliamsOtto()
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
            obj.delta_mat = cell2mat(struct2cell(obj.delta_point))';
            obj.decay = false;
            obj.time = 0;
            obj.forget = false;
            obj.op_region_script = @op_region_plot_WO;
            obj.x_a_val = 0.12;
        end
        
        function objective = objective_value(obj, inputs, outputs)
            % Calculate objective function
            flowrate_b = inputs(fieldnames(obj.feasible_point)=="flowrate_b");
            flowrate_out = flowrate_b + obj.constants.flowrate_a;
            x_p = outputs(fieldnames(obj.init_var)=="x_p");
            x_e = outputs(fieldnames(obj.init_var)=="x_e");
            objective = -( ...
                1043.38 * x_p * flowrate_out + ...
                20.92 * x_e * flowrate_out - ...
                79.23 * obj.constants.flowrate_a - 118.34 * flowrate_b ...
                );
        end
        
        function constraint = x_a_con(obj, ~, outputs)
            % x_a < 0.12
            x_a = outputs(fieldnames(obj.init_var)=="x_a");
            constraint = x_a - obj.x_a_val;
        end
        
        function constraint = x_g_con(obj, ~, outputs)
            % x_g < 0.08
            x_g = outputs(fieldnames(obj.init_var)=="x_g");
            constraint = x_g - 0.08;
        end
        
        function obj = time_increment(obj)
            if obj.decay
                obj.time = obj.time + 1;
            end
        end
        
        function [c,ceq] = nonlin_con(obj, x, par)
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
            par.values_adj = GPobj.values_adj(end);
            par.centre = GPobj.centre(idx, :);
            par.delta = GPobj.delta(idx, :);
%             par.delta_norm = cell2mat(struct2cell(obj.delta_norm))';
        end
        
        function [objective, con_ineq, con_eq] = get_output(obj, inputs)
            % Get objective function and constraint values using given inputs
            x0 = squeeze(cell2mat(struct2cell(obj.init_var)));
            option = optimset( ...
                'TolFun', 1e-8, 'TolX', 1e-8, 'Algorithm', ...
                'trust-region-dogleg', 'Display', 'off' ...
                );  % or 'Display','iter'
            
            n = size(inputs, 1);
            m = length(obj.constraints_ineq);
            k = length(obj.constraints_eq);
            objective = zeros(n, 1);
            con_ineq = zeros(n, m);
            con_eq = zeros(n, k);
            
            % Solve reactor
            for point = 1:n
                curr_inputs = inputs(point, :);
                [x, ~] = fsolve(@obj.system_of_eqns, x0, option, curr_inputs);
                objective(point) = obj.objective_value(curr_inputs, x);
                for idx = 1:m
                    val = obj.(obj.constraints_ineq{idx})(curr_inputs, x);
                    con_ineq(point, idx) = val;
                end
                for idx = 1:k
                    val = obj.(obj.constraints_eq{idx})(curr_inputs, x);
                    con_eq(point, idx) = val;
                end
            end
        end
        
        function [eqns] = system_of_eqns(obj, x, inputs)
            % Mass fractions
            x_a = x(1);
            x_b = x(2);
            x_c = x(3);
            x_e = x(4);
            x_g = x(5);
            x_p = x(6);
            M_total = obj.constants.total_molar_flowrate;
            flow_a = obj.constants.flowrate_a;
            
            % Logical arrays for getting indices
            temp = (fieldnames(obj.feasible_point)=="t_reactor");
            flow_b = (fieldnames(obj.feasible_point)=="flowrate_b");
            
            % Total outlet flowrate
            flow_total = flow_a + inputs(flow_b);
            
            % Rate constants
            k1 = 1.6599e6 * exp(-6666.7 / (inputs(temp) + 273.15));
            k2 = 7.2117e8 * exp(-8333.3 / (inputs(temp) + 273.15));
            k3 = 2.6745e12 * exp(-11111 / (inputs(temp) + 273.15));
            
            % Rates of reaction
            r1 = k1 * x_a * x_b * M_total;
            r2 = k2 * x_b * x_c * M_total;
            r3 = k3 * x_c * x_p * M_total;
            
            % System of equations
            eqns(1,1) = (flow_a - r1 - flow_total * x_a)/M_total;
            eqns(2,1) = (inputs(flow_b) - r1 - r2 - flow_total * x_b)/M_total;
            eqns(3,1) = (2 * r1 - 2 * r2 - r3 - flow_total * x_c)/M_total;
            eqns(4,1) = (2 * r2 - flow_total * x_e)/M_total;
            eqns(5,1) = (1.5 * r3 - flow_total * x_g)/M_total;
            eqns(6,1) = (r2 - 0.5 * r3 - flow_total * x_p)/M_total;
        end
    end
end