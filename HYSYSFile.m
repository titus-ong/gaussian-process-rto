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
        feasible_point = struct( ...
            "reboiler_duty", 170000, ...
            "solvent_flowrate", 250 ...
            );
        delta = struct( ...  10% of operating range
            "reboiler_duty", 300, ...
            "solvent_flowrate", 15 ...
            );
        output_fields = ["co2_recovery"];
        objective = @(x, mean_x, std_x, model) x(1);  % Minimise reboiler duty
        linear_con_A = [];                     % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                     % Linear inequality constraints LHS
        lineq_con_A = [];                      % Linear equality constraints LHS
        lineq_con_b = [];                      % Linear equality constraints LHS
        lb = [150000, 150];                    % Lower bounds
        ub = [200000, 300];                    % Upper bounds
        options = optimset('disp','iter');     % Options for GP
    end
    properties
        feasible_point_mat % Matrix form
        delta_mat          % Matrix form
        nonlin_con         % Nonlinear constraints
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
            obj.nonlinfn(x, centre, delta, model, mean_x, std_x, mean_y, std_y);
        end
        
        function [c,ceq] = nonlinfn(obj, x, centre, delta, model, mean_x, std_x, mean_y, std_y)
            % Nonlinear inequality (c<=0) and equality (ceq=0) constraints on x
            
            % Point is within delta
            c(1) = sum((x - centre).^2) - sum((delta).^2);
            
            % CO2 recovery > 0.99
            c(2) = obj.co2_fn(x, mean_x, std_x, mean_y, std_y, model);
            
            ceq = [];
        end
        
        function recovery = co2_fn(~, x, mean_x, std_x, mean_y, std_y, model)
            predicted = predict(model, (x - mean_x ./ std_x)) .* std_y + mean_y;
            recovery = 0.99 - predicted(1);
        end
        
        function value = get_param(obj, parameter)
            % Get parameter from HYSYS spreadsheet
            ref = obj.cells.(parameter);
            value = obj.spreadsheet.Cell(ref).CellValue;
        end
        
        function set_param(obj, parameter, value)
            % Set parameter on HYSYS spreadsheet
            ref = obj.cells.(parameter);
            obj.spreadsheet.Cell(ref).CellValue = value;
        end
        
        function stop_solver(obj)
            obj.solver.CanSolve = 0;
        end
        
        function start_solver(obj)
            obj.solver.CanSolve = 1;
        end
        
        function [outputs, outputs_mat] = get_output(obj, inputs)
            % Get output values from HYSYS spreadsheet using given inputs
            outputs = struct();
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
                    outputs(point).(obj.output_fields{idx}) = obj.get_param(obj.output_fields{idx});
                end
            end
            outputs_mat = squeeze(cell2mat(struct2cell(outputs)));
            if size(outputs_mat, 2) ~= length(obj.output_fields)
                outputs_mat = outputs_mat';
            end
        end
    end
end