classdef GPCreator
    properties
        model       % Either HYSYSFile or model object
        init_input  % Initial training input data
        init_output % Initial training output data
    end
    methods
        function obj = GPCreator(model, training_input, training_output)
            obj.model = model;
            obj.init_input = training_input;
            obj.init_output = training_output;
        end
    end
end