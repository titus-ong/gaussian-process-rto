%% User name variables
hysys_name = "\fast run.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";

%% Connect to HYSYS
tic
hysys_path = pwd + hysys_name;
hysys = HYSYSFile_fastrun(hysys_path, spreadsheet_input, spreadsheet_output);
toc

%% Initialising Points for Training
pointer = Pointer(hysys.feasible_point_mat, hysys.delta_mat);
sample_points = pointer.random_sampling();
toc
[objective, con_ineq, con_eq] = hysys.get_output(sample_points);
outputs = struct();
outputs.objective = objective;
outputs.con_ineq = con_ineq;
outputs.con_eq = con_eq;
toc

GP = GPCreator(hysys, sample_points, outputs);
toc

%% Optimisation
iter = 10;
GP.optimise(iter);
GP.plot();