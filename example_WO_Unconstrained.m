%% Load Williams-Otto
tic
WO = WilliamsOttoUnconstrained();
% WO.decay = true;  % Set to see effect of decay on optimisation

%% Initialising Points for Training
pointer = Pointer(WO.feasible_point_mat, WO.delta_mat);
sample_points = pointer.random_sampling();
% Fb = [6.9000, 6.7675, 7.0641, 7.1481, 6.9861, 6.6515]';
% Tr = [83.0000, 85.0262, 83.6761, 81.5130, 79.5847, 80.4423]';
% sample_points = [Fb, Tr];
[outputs, outputs_struct] = WO.get_output(sample_points);

GP = GPCreator(WO, sample_points, outputs);

%% Optimisation
for i = 1:50
    GP.optimise(1);
%     GP.system.time = GP.system.time + 1;
end
toc
GP.plot();
toc