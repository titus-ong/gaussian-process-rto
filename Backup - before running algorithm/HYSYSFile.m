classdef HYSYSFile  < handle
    properties (SetAccess=private)
        simulation                                 % Simulation object
        flowsheet                                  % Simulation flowsheet
        operations                                 % Simulation operations
        solver                                     % Simulation solver
        ss_input                                   % Simulation input spreadsheet
        ss_output                                  % Simulation output spreadsheet
    end
    properties (Constant)
        cells = struct( ...                        % Cell struct of spreadsheet
            'solvent_flowrate', 'A1', ...  % Var_name, {input, output}
            'reboiler_duty', 'A2', ...
            'inlet_co2_comp', 'A3', ...
            'inlet_gas_flowrate', 'A4', ...
            'clean_gas_co2', 'A5', ...
            'p_absorber', 'A6', ...
            'pout_absorber', 'A7', ...
            'p_stripper', 'A8', ...
            'pout_stripper', 'A9', ...
            't_condenser', 'A10', ...
            't_mea_recycle', 'A11', ...
            't_inlet_stripper', 'A12', ...
            't_inlet_gas', 'A13', ...
            'recycle_co2', 'A14', ...
            'recycle_mea', 'A15', ...
            'recycle_h2o', 'A16', ...
            'e_cooling_water', 'A17', ...
            'e_c100', 'A18', ...
            'e_c200', 'A19', ...
            'e_j100', 'A20', ...
            'e_j101', 'A21', ...
            'objective_true', 'A22' ...
            );
        feasible_point = struct( ...
            "reboiler_duty", 115000, ...
            "solvent_flowrate", 250 ...
            );
        delta = struct( ...                        % 10% of operating range
            "reboiler_duty", 300, ...
            "solvent_flowrate", 15 ...
            );
        output_fields = { ...
            'solvent_flowrate',  ... 
            'reboiler_duty',  ...
            'inlet_co2_comp', ...
            'inlet_gas_flowrate',  ...
            'clean_gas_co2',  ...
            'p_absorber',  ...
            'pout_absorber',  ...
            'p_stripper',  ...
            'pout_stripper',  ...
            't_condenser',  ...
            't_mea_recycle',  ...
            't_inlet_stripper', ...
            't_inlet_gas',  ...
            'recycle_co2',  ...
            'recycle_mea',  ...
            'recycle_h2o',  ...
            'e_cooling_water', ...
            'e_c100',  ...
            'e_c200',  ...
            'e_j100',  ...
            'e_j101' ...
            };
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [108000, 150];                        % Lower bounds
        ub = [120000, 300];                        % Upper bounds
        options = optimset('disp','iter');         % Options for GP
    end
    properties
        feasible_point_mat                     % Matrix form
        delta_mat                              % Matrix form
        nonlin_con                             % Nonlinear constraints
        objective_model                        % Objective function for model
        objective_true                         % Objective function for system
        system_violation                       % Test for system constraint violation
    end
    methods
        function obj = HYSYSFile(filepath, spreadsheet_input, spreadsheet_output) % Constructor method
            hysys = actxserver('HYSYS.APPLICATION');
            obj.simulation = hysys.SimulationCases.Open(filepath);
            obj.flowsheet = obj.simulation.Flowsheet;
            obj.operations = obj.flowsheet.Operations;
            obj.solver = obj.simulation.Solver;
            obj.ss_input = obj.operations.Item(spreadsheet_input);
            obj.ss_output = obj.operations.Item(spreadsheet_output);
            
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
        
        function bool = system_violated(obj, output)
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
            ref = obj.cells.(parameter);
            value = obj.ss_output.Cell(ref).CellValue;
        end
        
        function set_param(obj, parameter, value)
            % Set parameter on HYSYS spreadsheet
            ref = obj.cells.(parameter);
            obj.ss_input.Cell(ref).CellValue = value;
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
                tic
                obj.stop_solver();
                
                % Set parameters
                field = fields(obj.feasible_point);
                for idx = 1:length(field)
                    obj.set_param(field{idx}, inputs(point, idx));
                end
                toc
                obj.start_solver();
                toc
                % Get parameters
                for idx = 1:length(obj.output_fields)
%                     a = obj.get_param(obj.output_fields{idx});
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