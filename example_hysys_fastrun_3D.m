%% User name variables
hysys_name = "\System Files\fast run - 3D.hsc";
spreadsheet_input = "Inputs";
spreadsheet_output = "Outputs";
addpath 'System Files'\
addpath 'System Files\Operating region plots'\

%% Connect to HYSYS
tic
hysys_path = pwd + hysys_name;
hysys = HYSYSFile_fastrun_3D(hysys_path, spreadsheet_input, spreadsheet_output);
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
GP.excite = false;
toc

%% Optimisation
iter = 20;
GP.optimise(iter);
toc
% GP.plot();
%[rb, solvent, gas]
% rows = 20;
% cols = 20;
% heights = 20;
% 
% x1 = linspace(60000, 130000, cols);
% x2 = linspace(200, 600, rows);
% x3 = linspace(40, 120, heights);

% centres=GP.centre;
% centres(:,1)=(centres(:,1)-50000)/10000;
% centres(:,2)=(centres(:,2)-150)/50;
% centres(:,3)=(centres(:,3)-30)/10;


centres=GP.centre;
centres(:,1)=(centres(:,1)-95000)/35000;
centres(:,2)=(centres(:,2)-400)/200;
centres(:,3)=(centres(:,3)-80)/40;
centres = (centres+1)*19/2+1;
% plot3(centres(:,1),centres(:,2),centres(:,3),'-','Linewidth',1.5)
% scatter3(centres(:,1),centres(:,2),centres(:,3),[],GP.fval_min,'filled')