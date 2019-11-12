%% Load Williams-Otto
tic
WO = WilliamsOtto();
% WO.decay = true;  % Set to see effect of decay on optimisation

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
sample_points = pointer.random_sampling();
% Fb = [6.9000, 7.1722, 6.8980, 6.6067, 6.8077, 6.7585]';
% Tr = [83.0000, 84.7575, 80.9491, 83.9705, 80.2874, 81.9378]';
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