% Run example script first to get GP variable in workspace
rows = 100;
cols = 100;
heights = 100;
model_idx = 20;  % Which iteration of GP to display

% x1 = linspace(hysys.lb(1), hysys.ub(1), cols);
% x2 = linspace(hysys.lb(2), hysys.ub(2), rows);
% x3 = linspace(hysys.lb(3), hysys.ub(3), heights);
x1 = linspace(60000, 130000, cols);
x2 = linspace(200, 550, rows);
x3 = linspace(60, 100, heights);

temp_GP = copy(GP);
temp_GP.model = temp_GP.model(1:model_idx);
temp_GP.values_adj = temp_GP.values_adj(1:model_idx);
temp_GP.excited = temp_GP.excited(1:model_idx);

% Pre-allocation
data = zeros(rows, cols, heights);
[X1, X2, X3] = meshgrid(x1, x2, x3);
for i = 1:length(hysys.constraints_ineq)
    val_ineq.(hysys.constraints_ineq{i}) = zeros(rows, cols);
end
for i = 1:length(hysys.constraints_eq)
    val_eq.(hysys.constraints_eq{i}) = zeros(rows, cols);
end

% GP model prediction
for i = 1:cols
    for j = 1:rows
        for k = 1:heights
            data(j, i, k) = temp_GP.obj_fn([x1(i), x2(j), x3(k)]);
            for l = 1:length(hysys.constraints_ineq)
                val_ineq.(hysys.constraints_ineq{l})(j, i, k) = constraint(temp_GP, [x1(i), x2(j), x3(k)], hysys.constraints_ineq{l});
            end
            for l = 1:length(hysys.constraints_eq)
                val_eq.(hysys.constraints_eq{l})(j, i, k) = constraint(temp_GP, [x1(i), x2(j), x3(k)], hysys.constraints_eq{l});
            end
        end
    end
end

% Plots
f = figure;
hold on;
xslice = [100000 120000];                               % define the cross sections to view
yslice = [500];
zslice = [70 90];

slice(x1, x2, x3, data, xslice, yslice, zslice)    % display the slices

% xlim([80000 120000])
% ylim([200 500])
% zlim([70 100])

cb = colorbar;                                  % create and label the colorbar
cb.Label.String = 'Obj fn';

% % surf(x1, x2, data)
% contour(x1, x2, data, 30)

% for i = 1:length(hysys.constraints_ineq)
%     contour(x1, x2, val_ineq.(hysys.constraints_ineq{i}), [0, 0], 'k-');
% end
% for i = 1:length(hysys.constraints_eq)
%     contour(x1, x2, val_eq.(hysys.constraints_eq{i}), [0, 0], 'm-');
% end

training = scatter3(temp_GP.training_input(1:6,1), temp_GP.training_input(1:6,2), temp_GP.training_input(1:6,3), 'm');
centres = scatter3(temp_GP.centre(1:model_idx, 1), temp_GP.centre(1:model_idx, 2), temp_GP.centre(1:model_idx, 3));
legend([centres training], {'Centres', 'Training inputs'});
if sum(temp_GP.excited)
    for i = 1:length(temp_GP.excited)
        if ~temp_GP.excited(i)
            continue
        end
        optimas = plot(temp_GP.opt_min(i-1:i+1, 1), temp_GP.opt_min(i-1:i+1, 2), '-ro');
    end
    plot(temp_GP.opt_min(1:model_idx, 1), temp_GP.opt_min(1:model_idx, 2), ':k+');
    legend([centres training optimas], {'Centres', 'Training inputs', 'Excited points'});
end

function value = constraint(obj, x, constraint_name)
            value = predict( ...
                obj.model(end).(constraint_name), ...
                    (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                    ) * obj.values_adj(end).(constraint_name).std + obj.values_adj(end).(constraint_name).mean;
end