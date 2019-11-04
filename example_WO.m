%% Load Williams-Otto
tic
WO = WilliamsOtto();

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
sample_points = pointer.random_sampling();
% Fb   = [6.9, 5.7, 6.35, 6.6, 6.75]';%
% Tr   = [83, 74, 74.9, 75, 79]';%
% sample_points = [Fb, Tr];
[outputs, outputs_struct] = WO.get_output(sample_points);

GP = GPCreator(WO, sample_points, outputs);

%% Optimisation
iter = 50;
GP.optimise(iter);
toc
GP.plot();
toc