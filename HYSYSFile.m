classdef HYSYSFile
    properties (SetAccess=private)
        simulation       % Simulation object
        flowsheet        % Simulation flowsheet
        operations       % Simulation operations
        solver           % Simulation solver
        spreadsheet      % Simulation spreadsheet
    end
    properties (Constant)
        cells = struct( ...
            'solvent_flowrate', 'C2', ...
            'reboiler_duty', 'C3', ...
            'co2_recovery', 'C4' ...
            );           % Cell struct of spreadsheet
    end
    methods
        function obj = HYSYSFile(filepath, spreadsheet_name) % Constructor method
            hysys = actxserver('HYSYS.APPLICATION');
            obj.simulation = hysys.SimulationCases.Open(filepath);
            obj.flowsheet = obj.simulation.Flowsheet;
            obj.operations = obj.flowsheet.Operations;
            obj.solver = obj.simulation.Solver;
            obj.spreadsheet = obj.operations.Item(spreadsheet_name);
        end
        
        function value = get_param(obj, parameter)
            ref = obj.cells.(parameter);
            value = obj.spreadsheet.Cell(ref).CellValue;
        end
        
        function set_param(obj, parameter, value)
            ref = obj.cells.(parameter);
            obj.spreadsheet.Cell(ref).CellValue = value;
        end
        
        function stop_solver(obj)
            obj.solver.CanSolve = 0;
        end
        
        function start_solver(obj)
            obj.solver.CanSolve = 1;
        end
        
        function outputs = get_output(obj, sample_points, output_fields)
            outputs = struct();
            for point = 1:length(sample_points)
                obj.stop_solver();
                
                % Set parameters
                field = fields(sample_points);
                for idx = 1:length(field)
                    obj.set_param(field{idx}, sample_points(point).(field{idx}));
                end
                
                obj.start_solver();
                
                % Get parameters
                for idx = 1:length(output_fields)
                    outputs(point).(output_fields{idx}) = obj.get_param(output_fields{idx});
                end
            end
        end
    end
end