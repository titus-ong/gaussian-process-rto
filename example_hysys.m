%% Connect to HYSYS
tic
hysys_path = pwd + "\testing.hsc";
hysys = HYSYSFile(hysys_path, "SPRDSHT-1");
toc

%% Initialising Points for Training
pointer = Pointer(hysys.feasible_point_mat, hysys.delta_mat);
sample_points = pointer.random_sampling(10);
[outputs, outputs_mat] = hysys.get_output(sample_points);
toc

GP = GPCreator(hysys, sample_points, outputs_mat);

%% Optimisation
iter = 15;
GP.optimise(iter);
GP.plot();