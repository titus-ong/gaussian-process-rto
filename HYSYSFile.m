classdef HYSYSFile  < handle
    properties (SetAccess=private)
        simulation                                 % Simulation object
        flowsheet                                  % Simulation flowsheet
        operations                                 % Simulation operations
        solver                                     % Simulation solver
        spreadsheet                                % Simulation spreadsheet
    end
    properties (Constant)
        cells = struct( ...                        % Cell struct of spreadsheet
            'solvent_flowrate', {'A1', 'B1'}, ...  % Var_name, {input, output}
            'reboiler_duty', {'A2', 'B2'}, ...
            'inlet_co2_comp', {'A3', 'B3'}, ...
            'inlet_gas_flowrate', {'A4', 'B4'}, ...
            'clean_gas_co2', {'A5', 'B5'}, ...
            'p_absorber', {'A6', 'B6'}, ...
            'pout_absorber', {'A7', 'B7'}, ...
            'p_stripper', {'A8', 'B8'}, ...
            'pout_stripper', {'A9', 'B9'}, ...
            't_condenser', {'A10', 'B10'}, ...
            't_mea_recycle', {'A11', 'B11'}, ...
            't_inlet_stripper', {'A12', 'B12'}, ...
            't_inlet_gas', {'A13', 'B13'}, ...
            'recycle_co2', {'A14', 'B14'}, ...
            'recycle_mea', {'A15', 'B15'}, ...
            'recycle_h2o', {'A16', 'B16'}, ...
            'e_cooling_water', {'A17', 'B17'}, ...
            'e_c100', {'A18', 'B18'}, ...
            'e_c200', {'A19', 'B19'}, ...
            'e_j100', {'A20', 'B20'}, ...
            'e_j101', {'A21', 'B21'}, ...
            'objective_true', 'A22' ...
            );
        feasible_point = struct( ...
            "reboiler_duty", 154400, ...
            "solvent_flowrate", 250 ...
            );
        delta = struct( ...                        % 10% of operating range
            "reboiler_duty", 300, ...
            "solvent_flowrate", 15 ...
            );
        output_fields = ["objective_true"];
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [150000, 150];                        % Lower bounds
        ub = [200000, 300];                        % Upper bounds
        options = optimset('disp','iter');         % Options for GP
    end
    properties
        feasible_point_mat                     % Matrix form
        delta_mat                              % Matrix form
        nonlin_con                             % Nonlinear constraints
        objective_model                        % Objective function for model
        objective_true                         % Objective function for system
    end
    methods
        function obj = HYSYSFile(filepath, spreadsheet_name) % Constructor method
            hysys = actxserver('HYSYS.APPLICATION');
            obj.simulation = hysys.SimulationCases.Open(filepath);
            obj.flowsheet = obj.simulation.Flowsheet;
            obj.operations = obj.flowsheet.Operations;
            obj.solver = obj.simulation.Solver;
            obj.spreadsheet = obj.operations.Item(spreadsheet_name);
            
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
            obj.delta_mat = cell2mat(struct2cell(obj.delta))';
            obj.nonlin_con = @(x, centre, delta, model, mean_x, std_x, mean_y, std_y) ...
                obj.nonlin_fn(x, centre, delta, model, mean_x, std_x, mean_y, std_y);
            obj.objective_model = @(x, model, mean_x, std_x, mean_y, std_y) ...
                obj.model_obj_fn(x, model, mean_x, std_x, mean_y, std_y);
            obj.objective_true = @(x) obj.true_obj_fn(x);
        end
        
        function [c,ceq] = nonlin_fn(~, x, centre, delta, model, mean_x, std_x, mean_y, std_y)
            % Nonlinear inequality (c<=0) and equality (ceq=0) constraints on x
            
            % Point is within delta
            c(1) = sum((x - centre).^2) - sum((delta).^2);
            
            % CO2 recovery > 99
            function recovery = co2_fn(x, model, mean_x, std_x, mean_y, std_y)
                predicted = predict(model, (x - mean_x) ./ std_x) .* std_y + mean_y;
                recovery = predicted(1) - 1;  % Replace (1) to find co2 recovery according to output field
            end
            % TO BE ADDED WHEN CO2 WORKS IN HYSYS
            c(2) = co2_fn(x, model, mean_x, std_x, mean_y, std_y);
            
            ceq = [];
        end
        
        function bool = system_violation(obj, output)
            % Test if output violates system constraints
            % Return true if violated, false if not
            bool = false;
        end
        
        function objective = model_obj_fn(obj, x, model, mean_x, std_x, mean_y, std_y)
            % Calculate objective function for model
            outputs = predict(model, (x - mean_x) ./ std_x) .* std_y + mean_y;
            objective = obj.calc_objective(outputs);
        end
        
        function objective = true_obj_fn(obj, x)
            % Calculate objective function for system
            [outputs, ~] = obj.get_output(x);
            objective = obj.calc_objective(outputs);
        end
        
        function objective = calc_objective(~, outputs)
            % Calculate objective function from outputs
            objective = outputs(1);
        end
        
        function value = get_param(obj, parameter)
            % Get parameter from HYSYS spreadsheet
            ref = obj.cells(2).(parameter);
            value = obj.spreadsheet.Cell(ref).CellValue;
        end
        
        function set_param(obj, parameter, value)
            % Set parameter on HYSYS spreadsheet
            ref = obj.cells(1).(parameter);
            obj.spreadsheet.Cell(ref).CellValue = value;
        end
        
        function stop_solver(obj)
            obj.solver.CanSolve = 0;
        end
        
        function start_solver(obj)
            obj.solver.CanSolve = 1;
        end
        
        function [outputs, outputs_struct] = get_output(obj, inputs)
            % Get output values from HYSYS spreadsheet using given inputs
            outputs_struct = struct();
            for point = 1:size(inputs, 1)
                obj.stop_solver();
                
                % Set parameters
                field = fields(obj.feasible_point);
                for idx = 1:length(field)
                    obj.set_param(field{idx}, inputs(point, idx));
                end
                
                obj.start_solver();
                
                % Get parameters
                for idx = 1:length(obj.output_fields)
                    outputs_struct(point).(obj.output_fields{idx}) = obj.get_param(obj.output_fields{idx});
                end
            end
            outputs = squeeze(cell2mat(struct2cell(outputs_struct)));
            if size(outputs, 2) ~= length(obj.output_fields)
                outputs = outputs';
            end
        end
    end
end