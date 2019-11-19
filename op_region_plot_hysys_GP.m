% Run example script first to get GP variable in workspace
rows = 16;
cols = 31;
model_idx = 11;  % Which iteration of GP to display

x1 = linspace(hysys.lb(1), hysys.ub(1), cols);
x2 = linspace(hysys.lb(2), hysys.ub(2), rows);
temp_GP = copy(GP);
temp_GP.model = temp_GP.model(1:model_idx);
temp_GP.values_adj = temp_GP.values_adj(1:model_idx);

% Pre-allocation
data = zeros(rows, cols);
for i = 1:length(hysys.constraints_ineq)
    val_ineq.(hysys.constraints_ineq{i}) = zeros(rows, cols);
end
for i = 1:length(hysys.constraints_eq)
    val_eq.(hysys.constraints_eq{i}) = zeros(rows, cols);
end

% GP model prediction
for i = 1:cols
    for j = 1:rows
        data(j, i) = temp_GP.obj_fn([x1(i), x2(j)]);
        for k = 1:length(hysys.constraints_ineq)
            val_ineq.(hysys.constraints_ineq{k})(j, i) = constraint(temp_GP, [x1(i), x2(j)], hysys.constraints_ineq{k});
        end
        for k = 1:length(hysys.constraints_eq)
            val_eq.(hysys.constraints_eq{k})(j, i) = constraint(temp_GP, [x1(i), x2(j)], hysys.constraints_eq{k});
        end
    end
end

% Plots
f = figure;
hold on;
% surf(x1, x2, data)
contour(x1, x2, data, 30)

for i = 1:length(hysys.constraints_ineq)
    contour(x1, x2, val_ineq.(hysys.constraints_ineq{i}), [0, 0], 'k-');
end
for i = 1:length(hysys.constraints_eq)
    contour(x1, x2, val_eq.(hysys.constraints_eq{i}), [0, 0], 'm-');
end

training = scatter(temp_GP.training_input(1:6,1), temp_GP.training_input(1:6,2), '+r');
centres = plot(temp_GP.centre(1:model_idx, 1), temp_GP.centre(1:model_idx, 2), '-b*');
optimas = plot(temp_GP.opt_min(1:model_idx, 1), temp_GP.opt_min(1:model_idx, 2), '-go');
legend([centres training optimas], {'Centres', 'Training inputs', 'Optimised points'});

function value = constraint(obj, x, constraint_name)
            value = predict( ...
                obj.model(end).(constraint_name), ...
                    (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                    ) * obj.values_adj(end).(constraint_name).std + obj.values_adj(end).(constraint_name).mean;
end