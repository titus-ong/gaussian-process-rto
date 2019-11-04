%% Load Williams-Otto
tic
WO = WilliamsOtto();
% WO.decay = true;  % Set to see effect of decay on optimisation

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
% sample_points = pointer.random_sampling(5);
Fb = [6.9000, 6.6936, 7.0925, 7.0459, 7.8023, 7.1124, 6.8496, 6.3801, 7.1698, 6.6472, 6.2339]';
Tr = [83.000, 75.119, 88.384, 80.046, 78.908, 77.458, 74.059, 76.123, 86.019, 85.985, 89.158]';
sample_points = [Fb, Tr];
[outputs, outputs_struct] = WO.get_output(sample_points);

GP = GPCreator(WO, sample_points, outputs);

%% Optimisation
iter = 10;
GP.optimise(iter);
toc
GP.plot();
toc