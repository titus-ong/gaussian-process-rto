classdef GPCreator
    properties
        system      % Either HYSYSFile or system object
        model       % GP model
        
        init_input  % Initial training input data
        mean_input  % Mean of input
        std_input   % Standard deviation of input
        norm_input  % Normalised input
        
        init_output % Initial training output data
        mean_output % Mean of output
        std_output  % Standard deviation of output
        norm_output %Normalised output
    end
    methods
        function obj = GPCreator(system, training_input, training_output)
            obj.system = system;
            obj.init_input = training_input;
            obj.init_output = training_output;
            
            % Convert data to matrix
            cell_input = struct2cell(obj.init_input);
            mat_input = squeeze(cell2mat(cell_input))';
            cell_output = struct2cell(obj.init_output);
            mat_output = squeeze(cell2mat(cell_output))';
            
            % Normalise data
            obj.mean_input = mean(mat_input);
            obj.std_input = std(mat_input);
            obj.norm_input = normalize(mat_input);
            obj.mean_output = mean(mat_output);
            obj.std_output = std(mat_output);
            obj.norm_output = normalize(mat_output);
            
            % Initialise GP model
            obj.model = fitrgp(obj.norm_input, obj.norm_output);
        end
    end
end