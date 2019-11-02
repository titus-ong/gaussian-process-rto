%% Load Williams-Otto
tic
WO = WilliamsOtto();
toc

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
sample_points = pointer.random_sampling(10);
toc
[outputs, outputs_struct] = WO.get_output(sample_points);
toc

GP = GPCreator(WO, sample_points, outputs);
toc

%% Optimisation
iter = 20;
GP.optimise(iter);
GP.plot();