classdef WilliamsOtto
    properties (Constant)
        feasible_point = struct( ...
            "flowrate_b", 6, ...
            "t_reactor", 79 ...
            );
        delta = struct( ...
            "flowrate_b", 0.12, ...
            "t_reactor", 0.9 ...
            );
        output_fields = { ...
            'flowrate_out', ...
            'x_a', ...
            'x_b', ...
            'x_c', ...
            'x_e', ...
            'x_g', ...
            'x_p' ...
            };
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [4, 70];                              % Lower bounds
        ub = [7, 100];                             % Upper bounds
        options = optimset('disp','off');          % Options for GP
    end
    properties
        init_var = struct( ...
            "flowrate_a", 1.8275, ...                % Inlet flowrate of A
            "x_a_in", 1, ...                       % Inlet concentration of A
            "x_b_in", 1, ...                       % Inlet concentration of B
            "volume", 10, ...                      % Reactor volume
            "total_molar_flowrate", 2105.2, ...    % Total molar flowrate
            "flowrate_out", NaN, ...               % Outlet flowrate
            "x_a", 0.2, ...                        % Outlet concentration of A
            "x_b", 0.2, ...                        % Outlet concentration of B
            "x_c", 0.2, ...                        % Outlet concentration of C
            "x_e", 0.2, ...                        % Outlet concentration of E
            "x_g", 0.2, ...                        % Outlet concentration of G
            "x_p", 0.2 ...                         % Outlet concentration of P
            );
        
        feasible_point_mat                         % Matrix form
        delta_mat                                  % Matrix form
        nonlin_con                                 % Nonlinear constraints
        objective_model                            % Objective function for model
        objective_true                             % Objective function for system
        system_violation                           % Test for system constraint violation
    end
    methods
        function obj = WilliamsOtto()
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
            obj.delta_mat = cell2mat(struct2cell(obj.delta))';
            obj.nonlin_con = @(x, centre, delta, model, mean_x, std_x, mean_y, std_y) ...
                obj.nonlin_fn(x, centre, delta, model, mean_x, std_x, mean_y, std_y);
            obj.objective_model = @(x, model, mean_x, std_x, mean_y, std_y) ...
                obj.model_obj_fn(x, model, mean_x, std_x, mean_y, std_y);
            obj.objective_true = @(x) obj.true_obj_fn(x);
            obj.system_violation = @(y) obj.system_violated(y);
        end
        
        function [c,ceq] = nonlin_fn(~, x, centre, delta, model, mean_x, std_x, mean_y, std_y)
            % Nonlinear inequality (c<=0) and equality (ceq=0) constraints on x
            
            % x_a to be less than 0.12
            function percentage = constraint_x_a(x, model)
                predicted = predict(model.x_a, x);
                percentage = predicted - 0.12;
            end
            c(1) = constraint_x_a(x, model);
            
            % x_g to be less than 0.08
            function percentage = constraint_x_g(x, model)
                predicted = predict(model.x_g, x);
                percentage = predicted - 0.08;
            end
            c(2) = constraint_x_g(x, model);

            ceq = [];
        end
        
        function bool = system_violated(obj, output)
            % Test if output violates system constraints
            % Return true if violated, false if not
            bool = false;
        end
        
        function objective = model_obj_fn(obj, x, model, mean_x, std_x, mean_y, std_y)
            % Calculate objective function for model
            outputs = zeros(1, length(obj.output_fields));
            for i = 1:length(obj.output_fields)
%                 outputs(i) = predict(model(i), (x - mean_x) ./ std_x) .* std_y(i) + mean_y(i);
                outputs(i) = predict(model.(obj.output_fields{i}), x);
            end
            objective = obj.calc_objective(outputs);
        end
        
        function objective = true_obj_fn(obj, x)
            % Calculate objective function for system
            [outputs, ~] = obj.get_output(x);
            objective = obj.calc_objective(outputs);
        end
        
        function objective = calc_objective(obj, outputs)
            % Calculate objective function from outputs
            flowrate_out = (obj.output_fields=="flowrate_out");
            x_p = (obj.output_fields=="x_p");
            x_e = (obj.output_fields=="x_e");
            flowrate_b = outputs(flowrate_out) - obj.init_var.flowrate_a;
            objective = -( ...
                1043.38 * outputs(x_p) * outputs(flowrate_out) + ...
                20.92 * outputs(x_e) * outputs(flowrate_out) - ...
                79.23 * obj.init_var.flowrate_a - 118.34 * flowrate_b ...
                );
        end
        
        function [outputs, outputs_struct] = get_output(obj, inputs)
            % Get output values from HYSYS spreadsheet using given inputs
            outputs_struct = struct();
            x0_fields = { ...
                'x_a', ...
                'x_b', ...
                'x_c', ...
                'x_e', ...
                'x_g', ...
                'x_p', ...
                'flowrate_out' ...  % Has to be calculated for each input
                };
            x0 = [ ...
                obj.init_var.x_a, ...
                obj.init_var.x_b, ...
                obj.init_var.x_c, ...
                obj.init_var.x_e, ...
                obj.init_var.x_g, ...
                obj.init_var.x_p ...
                ];
            option = optimset('TolFun',1e-8,'TolX',1e-8,'Algorithm','trust-region-dogleg','Display','off');%'Display','iter',

            % Logical arrays for getting indices
            flow_b = (fieldnames(obj.feasible_point)=="flowrate_b");
            
            % Solve reactor
            for point = 1:size(inputs, 1)
                [x, ~] = fsolve(@obj.system_of_eqns, x0, option, inputs(point, :));
                x(end + 1) = inputs(flow_b) + obj.init_var.flowrate_a;
                for idx = 1:length(obj.output_fields)
                    field = convertCharsToStrings(obj.output_fields{idx});
                    outputs_struct(point).(obj.output_fields{idx}) = x(x0_fields==field);
                end
            end
            outputs = squeeze(cell2mat(struct2cell(outputs_struct)));
            if size(outputs, 2) ~= length(obj.output_fields)
                outputs = outputs';
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
            M_total = obj.init_var.total_molar_flowrate;
            flow_a = obj.init_var.flowrate_a;
            
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
            eqns(5,1) = (1.5*r3 - flow_total * x_g)/M_total;
            eqns(6,1) = (r2 - 0.5*r3 - flow_total * x_p)/M_total;
        end
    end
end