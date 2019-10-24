%% Connect to HYSYS
tic
hysys_path = pwd + "\testing.hsc";
hysys = HYSYSFile(hysys_path, "SPRDSHT-1");
toc

%% Initialising Points for Training
feasible_point = struct( ...
    "reboiler_duty", 250, ...
    "solvent_flowrate", 170000 ...
    );
delta = struct( ...  % 10% of operating range
    "reboiler_duty", 15, ...
    "solvent_flowrate", 300 ...
    );
output_fields = ["co2_recovery"];

pointer = Pointer(feasible_point, delta);
sample_points = pointer.random_sampling();
outputs = hysys.get_output(sample_points, output_fields);
toc

GPer = GPCreator(hysys, sample_points, outputs);

%% Optimisation
iter = 15;
GPer.optimise(iter);
GPer.plot();