%% Connect to HYSYS
tic
hysys_path = pwd + "\testing.hsc";
hysys = HYSYSFile(hysys_path, "SPRDSHT-1");
toc

%% Initialising Points for Training
pointer = Pointer(hysys.feasible_point, hysys.delta);
sample_points = pointer.random_sampling(10);
outputs = hysys.get_output(sample_points);
toc

GPer = GPCreator(hysys, sample_points, outputs);

%% Optimisation
iter = 15;
GPer.optimise(iter);
GPer.plot();