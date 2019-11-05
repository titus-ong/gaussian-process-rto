%% Load Williams-Otto
tic
WO = WilliamsOtto();
% WO.decay = true;  % Set to see effect of decay on optimisation

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
sample_points = pointer.random_sampling();
% Fb = [6.9000, 6.9621, 6.7991, 6.8912, 6.7755, 7.0513]';
% Tr = [83.0000, 81.9557, 84.5682, 84.7586, 85.8207, 80.8806]';
% sample_points = [Fb, Tr];
% Fb   = [6.9,5.7,6.35, 6.6, 6.75]';%
% Tr   = [83,74,74.9,75,79]';%
% sample_points = [Fb, Tr];
[objective, con_ineq, con_eq] = WO.get_output(sample_points);
outputs = struct();
outputs.objective = objective;
outputs.con_ineq = con_ineq;
outputs.con_eq = con_eq;

GP = GPCreator(WO, sample_points, outputs);

%% Optimisation
iter = 20;
GP.optimise(iter);
% for i = 1:20
%     GP.optimise(1);
%     GP.system.time = GP.system.time + 1;
% end
toc
GP.plot();
toc