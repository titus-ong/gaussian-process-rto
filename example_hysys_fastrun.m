%% User name variables
hysys_name = "\System Files\fast run new.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";
addpath 'System Files'\
addpath 'System Files\Operating region plots'\

%% Connect to HYSYS
tic
hysys_path = pwd + hysys_name;
hysys = HYSYSFile_fastrun(hysys_path, spreadsheet_input, spreadsheet_output);
toc

%% Initialising Points for Training
pointer = Pointer(hysys.feasible_point_mat, hysys.delta_mat, hysys.lb, hysys.ub);
sample_points = pointer.random_sampling();
[objective, con_ineq, con_eq] = hysys.get_output(sample_points);
outputs = struct();
outputs.objective = objective;
outputs.con_ineq = con_ineq;
outputs.con_eq = con_eq;

GP = GPCreator(hysys, sample_points, outputs);
toc

%% Optimisation
iter = 10;
GP.optimise(iter);
toc
GP.plot();