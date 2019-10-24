classdef Pointer
    properties
        feasible_point  % Feasible point (struct)
        delta           % Delta of point (struct)
        delta_values    % Delta of point (matrix)
        current_point   % Current point  (struct)
        current_values  % Current point  (matrix)
        dimension       % Dimensions of feasible point
    end
    properties (Constant)
        n_default = 5   % Number of sample points
    end
    methods
        function obj = Pointer(feasible_point, delta)
            obj.feasible_point = feasible_point;
            obj.delta = delta;
            obj.current_point = feasible_point;
            obj.dimension = length(fields(obj.current_point));
            obj.delta_values = zeros(1, obj.dimension);
            obj.current_values = zeros(1, obj.dimension);
        end
        
        function normaliser(obj)
            % do something to obj.feasible_point?
        end
        
        function samples = random_sampling(obj, varargin)
            % Return matrix of current point + n number of random samples
            % around the current point within delta
            
            % Get n from either default value or input
            p = inputParser();
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            addOptional(p, 'n', obj.n_default, validScalarPosNum);
            parse(p, varargin{:});
            
            samples = struct();
            field = fields(obj.current_point);
            
            % Include current point
            for col = 1:obj.dimension
                samples(1).(field{col}) = obj.current_point.(field{col});
            end
            
            % Random sampling within delta
            for point = 2:p.Results.n+1  % Leave a space for current point
                u = -1 + 2 * rand(1, obj.dimension);
                norm = (sum(u.^2))^0.5;
                r = rand(1, 1)^(1 / obj.dimension);
                for col = 1:obj.dimension
                    samples(point).(field{col}) = ( ...
                    obj.delta.(field{col}) * u(col) * r / norm ...
                    + obj.current_point.(field{col}) ...  % shift by current point
                    );
                end
            end
        end
    end
end