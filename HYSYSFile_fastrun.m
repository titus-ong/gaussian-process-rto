classdef HYSYSFile_fastrun  < matlab.mixin.Copyable
    properties (SetAccess=private)
        simulation                                 % Simulation object
        flowsheet                                  % Simulation flowsheet
        operations                                 % Simulation operations
        solver                                     % Simulation solver
        ss_input                                   % Simulation input spreadsheet
        ss_output                                  % Simulation output spreadsheet
        energy_stream                              % Simulation energy streams
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
            'sweet_co2_flowrate', 'A22', ...
            'recycle_co2_flowrate', 'A23', ...
            'co2_recovery_stripper', 'A24', ...
            'objective_true', 'A25' ...
            );
        feasible_point = struct( ...
            "reboiler_duty", 100000, ...
            "inlet_gas_flowrate", 66 ...
            );
        delta = struct( ...                        % 10% of operating range
            "reboiler_duty", 12000, ...
            "inlet_gas_flowrate", 26 ...
            );
        output_fields = { ...
            'solvent_flowrate', ... 
            'reboiler_duty', ...
            'inlet_co2_comp', ...
            'inlet_gas_flowrate', ...
            'clean_gas_co2', ...
            'p_absorber', ...
            'pout_absorber', ...
            'p_stripper', ...
            'pout_stripper', ...
            't_condenser', ...
            't_mea_recycle', ...
            't_inlet_stripper', ...
            't_inlet_gas', ...
            'recycle_co2', ...
            'recycle_mea', ...
            'recycle_h2o', ...
            'e_cooling_water', ...
            'e_c100', ...
            'e_c200', ...
            'e_j100', ...
            'e_j101', ...
            'sweet_co2_flowrate', ...
            'recycle_co2_flowrate', ...
            'co2_recovery_stripper', ...
            'objective_true' ...
            };
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [80000, 40];                     % Lower bounds
        ub = [200000, 300];                  % Upper bounds
        options = optimset('disp','off');         % Options for GP - iter or off
    end
    properties
        feasible_point_mat                     % Matrix form
        delta_mat                              % Matrix form
        nonlin_con                             % Nonlinear constraints
        objective_model                        % Objective function for model
        objective_true                         % Objective function for system
        system_violation                       % Test for system constraint violation
        op_region_script                       % Function for plotting operating region
    end
    methods
        function obj = HYSYSFile_fastrun(filepath, spreadsheet_input, spreadsheet_output) % Constructor method
            hysys = actxserver('HYSYS.APPLICATION');
            obj.simulation = hysys.SimulationCases.Open(filepath);
            obj.flowsheet = obj.simulation.Flowsheet;
            obj.operations = obj.flowsheet.Operations;
            obj.solver = obj.simulation.Solver;
            obj.ss_input = obj.operations.Item(spreadsheet_input);
            obj.ss_output = obj.operations.Item(spreadsheet_output);
            obj.energy_stream = obj.flowsheet.EnergyStreams;
            
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
            obj.delta_mat = cell2mat(struct2cell(obj.delta))';
            obj.nonlin_con = @(x, centre, delta, model, mean_x, std_x, mean_y, std_y) ...
                obj.nonlin_fn(x, centre, delta, model, mean_x, std_x, mean_y, std_y);
            obj.objective_model = @(x, model, mean_x, std_x, mean_y, std_y) ...
                obj.model_obj_fn(x, model, mean_x, std_x, mean_y, std_y);
            obj.objective_true = @(x) obj.true_obj_fn(x);
            obj.system_violation = @(y) obj.system_violated(y);
            obj.op_region_script = @op_region_plot_hysys;
        end
        
        function [c,ceq] = nonlin_fn(~, x, centre, delta, model, mean_x, std_x, mean_y, std_y)
            % Nonlinear inequality (c<=0) and equality (ceq=0) constraints on x
            
            % Point is within delta
            c(1) = sum((x - centre).^2 ./ delta.^2) - 1;
            
            % Clean [CO2] from absorber
%             function percentage = co2_fn(x, model, mean_x, std_x, mean_y, std_y)
% %                 logic_arr = (obj.output_fields=="clean_gas_co2");
%                 predicted = predict(model(end).clean_gas_co2, x);
% %                 predicted = predict( ...
% %                     model.clean_gas_co2, (x - mean_x) ./ std_x ...
% %                         ) .* std_y(logic_arr) + mean_y(logic_arr);
%                 percentage = predicted - 0.003;
%             end
% 
%             c(2) = co2_fn(x, model, mean_x, std_x, mean_y, std_y);
%             
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
            objective = outputs(obj.output_fields=="objective_true");
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
                obj.stop_solver();
                
                % Set parameters
                field = fields(obj.feasible_point);
                for idx = 1:length(field)
                    obj.set_param(field{idx}, inputs(point, idx));
                end
                
                obj.energy_stream.Item('Cooling water').HeatFlowvalue=-500/3600;
                
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