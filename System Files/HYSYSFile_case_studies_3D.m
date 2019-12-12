classdef HYSYSFile_case_studies_3D  < matlab.mixin.Copyable
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
            'solvent_flowrate', 'A1', ...          % Var_name, {input, output}
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
            'j100_flowrate', 'A25', ...
            'j101_flowrate', 'A26', ...
            'sourgas_CO2_massflow', 'A27', ...
            'objective_true', 'A28' ...
            );
        feasible_point = struct( ...
            "reboiler_duty", 180000, ...            
            'solvent_flowrate', 300, ...
            'inlet_gas_flowrate', 100 ...            
            );

        output_fields = { ...
            'solvent_flowrate', ...          % Var_name, {input, output}
            'reboiler_duty', ...
            'inlet_co2_comp',  ...
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
            'j100_flowrate', ...
            'j101_flowrate', ...
            'sourgas_CO2_massflow', ...
            'objective_true'...
            };
        lb = [80000, 40];                          % Lower bounds
        ub = [200000, 300];                        % Upper bounds
    end
    
    properties
        feasible_point_mat                         % Matrix form
    end
    
    methods
        function obj = HYSYSFile_case_studies_3D(filepath, spreadsheet_input, spreadsheet_output) % Constructor method
            hysys = actxserver('HYSYS.APPLICATION');
            obj.simulation = hysys.SimulationCases.Open(filepath);
            obj.flowsheet = obj.simulation.Flowsheet;
            obj.operations = obj.flowsheet.Operations;
            obj.solver = obj.simulation.Solver;
            obj.ss_input = obj.operations.Item(spreadsheet_input);
            obj.ss_output = obj.operations.Item(spreadsheet_output);
            obj.energy_stream = obj.flowsheet.EnergyStreams;
            
            obj.feasible_point_mat = cell2mat(struct2cell(obj.feasible_point))';
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
        
        function [outputs] = get_output(obj, inputs)
            % Get output values from HYSYS spreadsheet using given inputs
            
            % Pre-allocation
            outputs = zeros(length(inputs(:,1)), length(obj.output_fields));
            
            for point = 1:size(inputs, 1)
                obj.stop_solver();
                
                % Set parameters
                field = fields(obj.feasible_point);
                for idx = 1:length(field)
                    obj.set_param(field{idx}, inputs(point, idx));
                end
                

                obj.energy_stream.Item('Cooling water').HeatFlowvalue = -500/3600;

                
                obj.start_solver();
                
                % Get parameters
                for idx = 1:length(obj.output_fields)
                    try
                        outputs(point, idx) = obj.get_param(obj.output_fields{idx});
                    catch
                        outputs(point, idx) = NaN;
                    end
                end
            end
        end
    end
end