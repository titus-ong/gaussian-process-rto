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
            'objective_true', 'A25' ...
            );
        feasible_point = struct( ...
            "reboiler_duty", 180000, ...
            "inlet_gas_flowrate", 66 ...
            );
        delta = struct( ...                        % 10% of operating range
            "reboiler_duty", 12000, ...
            "inlet_gas_flowrate", 26 ...
            );
        constraints_ineq = { ...
            'sweetgas_CO2_comp', ...
            };
        
        constraints_eq = [];
        
        output_fields = { ...
            'objective_true', ...
            'reboiler_duty', ...
            'e_cooling_water' ...
            'clean_gas_co2' ...
            };
        linear_con_A = [];                         % Linear inequality constraints LHS (Ax = b)
        linear_con_b = [];                         % Linear inequality constraints LHS
        lineq_con_A = [];                          % Linear equality constraints LHS
        lineq_con_b = [];                          % Linear equality constraints LHS
        lb = [80000, 40];                          % Lower bounds
        ub = [200000, 300];                        % Upper bounds
        options = optimset('disp','off');          % Options for GP - iter or off
    
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
        op_region_script                           % Function for plotting operating region
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
            obj.op_region_script = @op_region_plot_hysys;
        end
        
        function objective = objective_value(obj, outputs)
            % Calculate objective function from outputs
            objective = outputs(obj.output_fields=="objective_true");
        end
        
        function constraint = sweetgas_CO2_comp(obj, outputs)
            % sweetgas_CO2_comp < 0.03
            sweetgas_CO2_comp = outputs(obj.output_fields=="clean_gas_co2");
            constraint = sweetgas_CO2_comp - 0.03;
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
            par.values_adj = GPobj.values_adj;
            par.centre = GPobj.centre(idx, :);
            par.delta = GPobj.delta(idx, :);
        end
        
        function bool = system_violated(~, con_ineq, con_eq)
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
        
        function [objective, con_ineq, con_eq] = get_output(obj, inputs)
            % Get output values from HYSYS spreadsheet using given inputs
            % and calculate objective function and constraint values
            
            % Pre-allocation
            outputs = zeros(length(obj.output_fields), 1);
            n = size(inputs, 1);
            m = length(obj.constraints_ineq);
            k = length(obj.constraints_eq);
            objective = zeros(n, 1);
            con_ineq = zeros(n, m);
            con_eq = zeros(n, k);
            
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
                    outputs(idx) = obj.get_param(obj.output_fields{idx});
                end

                % Get objective and constraint values
                objective(point) = obj.objective_value(outputs);
                for idx = 1:m
                    val = obj.(obj.constraints_ineq{idx})(outputs);
                    con_ineq(point, idx) = val;
                end
                for idx = 1:k
                    val = obj.(obj.constraints_eq{idx})(curr_inputs, x);
                    con_eq(point, idx) = val;
                end
            end
        end
    end
end