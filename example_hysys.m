%% User name variables
hysys_name = "\testing.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";

%% Connect to HYSYS
tic
hysys_path = pwd + hysys_name;
hysys = HYSYSFile(hysys_path, spreadsheet_input, spreadsheet_output);
toc

%% Initialising Points for Training
pointer = Pointer(hysys.feasible_point_mat, hysys.delta_mat);
sample_points = pointer.random_sampling(5);
toc
[outputs, outputs_struct] = hysys.get_output(sample_points);
toc

GP = GPCreator(hysys, sample_points, outputs);
toc

%% Optimisation
iter = 10;
GP.optimise(iter);
% GP.plot();