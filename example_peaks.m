%% Load Peaks
tic
Pks = PeaksFunction();

%% Initialising Points for Training
pointer = Pointer(Pks.feasible_point_mat, Pks.delta_mat);
sample_points = pointer.random_sampling();
[outputs, outputs_struct] = Pks.get_output(sample_points);

GP = GPCreator(Pks, sample_points, outputs);

%% Optimisation
iter = 1;
GP.optimise(iter);
toc
GP.plot();
toc