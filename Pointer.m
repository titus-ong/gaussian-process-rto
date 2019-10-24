classdef Pointer
    properties
        feasible_point  % Feasible point matrix
        delta           % Delta of point matrix
        current_point   % Current point matrix
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
            obj.dimension = length(feasible_point);
        end
        
        function normaliser(obj)
            % Normalise delta/dataset
        end
        
        function obj = update(obj, current, delta_)
            obj.current_point = current;
            obj.delta = delta_;
        end
        
        function samples = random_sampling(obj, varargin)
            % Return matrix of current point + n number of random samples
            % around the current point within delta
            
            % Get n from either default value or input
            p = inputParser();
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            addOptional(p, 'n', obj.n_default, validScalarPosNum);
            parse(p, varargin{:});
            
            samples = zeros(p.Results.n + 1, obj.dimension);
            
            % Include current point
            samples(1, :) = obj.current_point;
            
            % Random sampling within delta
            for point = 2:p.Results.n+1  % Leave a space for current point
                u = -1 + 2 * rand(1, obj.dimension);
                norm = (sum(u.^2))^0.5;
                r = rand(1, 1)^(1 / obj.dimension);
                samples(point, :) = obj.delta .* u * r / norm + obj.current_point;
            end
        end
    end
end