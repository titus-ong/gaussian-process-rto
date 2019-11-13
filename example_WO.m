%% Load Williams-Otto
tic
WO = WilliamsOtto();
% WO.decay = true;  % Set to see effect of decay on optimisation

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
sample_points = pointer.random_sampling();
% Sample points that get stuck:
% Fb = [6.9000, 6.6992, 6.9780, 6.9399, 6.7452, 6.7390]';
% Tr = [83.0000, 80.8881, 84.2125, 82.8227, 79.7888, 83.2433]';
% sample_points = [Fb, Tr];
[objective, con_ineq, con_eq] = WO.get_output(sample_points);
outputs = struct();
outputs.objective = objective;
outputs.con_ineq = con_ineq;
outputs.con_eq = con_eq;

GP = GPCreator(WO, sample_points, outputs);

%% Optimisation
iter = 10;
GP.optimise(iter);
% for i = 1:20
%     GP.optimise(1);
%     GP.system.time = GP.system.time + 1;
% end
toc
GP.plot();
toc