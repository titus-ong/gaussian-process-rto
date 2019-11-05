classdef PeaksFunction  < matlab.mixin.Copyable
    properties (Constant)
        feasible_point = struct( ...
            "x", -0.33, ...
            "y", 1 ...
            );
        delta = struct( ...
            "x", 0.6, ...
            "y", 0.6 ...
            );
        output_fields = {'z'};
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [-3, -3];                             % Lower bounds
        ub = [3, 3];                               % Upper bounds
        options = optimset('disp','off');          % Options for GP
    end
    properties
        feasible_point_mat                         % Matrix form
        delta_mat                                  % Matrix form
        nonlin_con                                 % Nonlinear constraints
        objective_model                            % Objective function for model
        objective_true                             % Objective function for system
        system_violation                           % Test for system constraint violation
        decay                                      % Decay boolean
        time                                       % Pseudo time for decay
    end
    methods
        function obj = PeaksFunction()
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
            obj.delta_mat = cell2mat(struct2cell(obj.delta))';
            obj.nonlin_con = @(x, centre, delta, model, mean_x, std_x, mean_y, std_y) ...
                obj.nonlin_fn(x, centre, delta, model, mean_x, std_x, mean_y, std_y);
            obj.objective_model = @(x, model, mean_x, std_x, mean_y, std_y) ...
                obj.model_obj_fn(x, model, mean_x, std_x, mean_y, std_y);
            obj.objective_true = @(x) obj.true_obj_fn(x);
            obj.system_violation = @(y) obj.system_violated(y);
            obj.decay = false;
            obj.time = 0;
        end
        
        function [c,ceq] = nonlin_fn(~, x, centre, delta, model, mean_x, std_x, mean_y, std_y)
            % Nonlinear inequality (c<=0) and equality (ceq=0) constraints on x
            
            % Point is within delta
            c(1) = sum((x - centre).^2 ./ delta.^2) - 1;

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
                outputs(i) = predict(model(end).(obj.output_fields{i}), x);
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
            objective = outputs(1);
        end
        
        function [outputs, outputs_struct] = get_output(obj, inputs)
            % Get output values from HYSYS spreadsheet using given inputs
            outputs_struct = struct();

            % Logical arrays for getting indices
            x = (fieldnames(obj.feasible_point)=="x");
            y = (fieldnames(obj.feasible_point)=="y");
            
            % Solve reactor
            for point = 1:size(inputs, 1)
                curr_inputs = inputs(point, :);
                z = peaks(curr_inputs(x), curr_inputs(y));
                outputs_struct(point).z = z;
            end
            outputs = squeeze(cell2mat(struct2cell(outputs_struct)));
            if size(outputs, 2) ~= length(obj.output_fields)
                outputs = outputs';
            end
        end
    end
end