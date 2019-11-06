%% Load Peaks
tic
Pks = PeaksFunction();
%Pks.decay = true;     % Set to see effect of decay on optimisation
%Pks.forget = true;     % Set to see effect of forgetting factor on optimisation
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
iter = 50;
GP.optimise(iter);
toc
GP.plot();
toc