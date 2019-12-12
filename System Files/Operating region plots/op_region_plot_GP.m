% Run example script first to get GP variable in workspace
rows = 16;
cols = 31;
model_idx = 45;  % Which iteration of GP to display

x1 = linspace(GP.lb(1), GP.ub(1), cols);
x2 = linspace(GP.lb(2), GP.ub(2), rows);
temp_GP = copy(GP);
temp_GP.model = temp_GP.model(1:model_idx);
temp_GP.values_adj = temp_GP.values_adj(1:model_idx);
temp_GP.excited = temp_GP.excited(1:model_idx);

% Pre-allocation
data = zeros(rows, cols);
for i = 1:length(temp_GP.system.constraints_ineq)
    val_ineq.(temp_GP.system.constraints_ineq{i}) = zeros(rows, cols);
end
for i = 1:length(temp_GP.system.constraints_eq)
    val_eq.(temp_GP.system.constraints_eq{i}) = zeros(rows, cols);
end

% GP model prediction
for i = 1:cols
    for j = 1:rows
        data(j, i) = temp_GP.obj_fn([x1(i), x2(j)]);
        for k = 1:length(temp_GP.system.constraints_ineq)
            val_ineq.(temp_GP.system.constraints_ineq{k})(j, i) = constraint(temp_GP, [x1(i), x2(j)], temp_GP.system.constraints_ineq{k});
        end
        for k = 1:length(temp_GP.system.constraints_eq)
            val_eq.(temp_GP.system.constraints_eq{k})(j, i) = constraint(temp_GP, [x1(i), x2(j)], temp_GP.system.constraints_eq{k});
        end
    end
end

% Plots
f = figure;
hold on;
% surf(x1, x2, data)
contour(x1, x2, data, 30)

for i = 1:length(temp_GP.system.constraints_ineq)
    contour(x1, x2, val_ineq.(temp_GP.system.constraints_ineq{i}), [0, 0], 'k-');
end
for i = 1:length(temp_GP.system.constraints_eq)
    contour(x1, x2, val_eq.(temp_GP.system.constraints_eq{i}), [0, 0], 'm-');
end

% % Plot last trust region
% syms x y a b h k
% a = temp_GP.delta(model_idx, 1);
% b = temp_GP.delta(model_idx, 2);
% h = temp_GP.centre(model_idx, 1);
% k = temp_GP.centre(model_idx, 2);
% ellipse = (((x-h)^2)/(a^2))+(((y-k)^2)/(b^2))==1;
% fimplicit(ellipse, [temp_GP.lb(1) temp_GP.ub(1) temp_GP.lb(2) temp_GP.ub(2)], '--b');

training = scatter(temp_GP.training_starter(:,1), temp_GP.training_starter(:,2), '+m');
centres = plot(temp_GP.centre(1:model_idx, 1), temp_GP.centre(1:model_idx, 2), '-b*');
legend([centres training], {'Centres', 'Training inputs'});

% Plot excited points
if sum(temp_GP.excited)
    for i = 1:length(temp_GP.excited)
        if ~temp_GP.excited(i)
            continue
        end
        optimas = plot(temp_GP.opt_min(i-1:i+1, 1), temp_GP.opt_min(i-1:i+1, 2), '-ro');
    end
%     plot(temp_GP.opt_min(1:model_idx, 1), temp_GP.opt_min(1:model_idx, 2), ':k+');
    legend([centres training optimas], {'Centres', 'Training inputs', 'Excited points'});
end

function value = constraint(obj, x, constraint_name)
            value = predict( ...
                obj.model(end).(constraint_name), ...
                    (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                    ) * obj.values_adj(end).(constraint_name).std + obj.values_adj(end).(constraint_name).mean;
end