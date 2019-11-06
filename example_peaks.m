%% Load Peaks
tic
Pks = PeaksFunction();

%% Initialising Points for Training
pointer = Pointer(Pks.feasible_point_mat, Pks.delta_mat);
sample_points = pointer.random_sampling();
[objective, con_ineq, con_eq] = Pks.get_output(sample_points);
outputs = struct();
outputs.objective = objective;
outputs.con_ineq = con_ineq;
outputs.con_eq = con_eq;

GP = GPCreator(Pks, sample_points, outputs);

%% Optimisation
iter = 1;
GP.optimise(iter);
toc
GP.plot();
toc