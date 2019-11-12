% Run example script first to get GP variable in workspace
rows = 16;
cols = 31;
x1 = linspace(hysys.lb(1), hysys.ub(1), cols);
x2 = linspace(hysys.lb(2), hysys.ub(2), rows);
temp_GP = copy(GP);
model_idx = 26;
temp_GP.model = temp_GP.model(1:model_idx);
temp_GP.values_adj = temp_GP.values_adj(1:model_idx);

data = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = temp_GP.obj_fn([x1(i), x2(j)]);
    end
end

data1 = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data1(j, i) = constraint(temp_GP,[x1(i), x2(j)]);
    end
end

f = figure;
hold on;
% surf(x1, x2, data)
%contour(x1, x2, data1,30)
[one_perc, h] = contour(x1, x2, data1, [0 0]);
clabel(one_perc, h)

scatter(temp_GP.training_input(1:6,1), temp_GP.training_input(1:6,2));

function thing = constraint(obj, x)
            thing = predict( ...
                obj.model(end).sweetgas_CO2_comp, ...
                    (x - obj.values_adj(end).input.mean) ./ obj.values_adj(end).input.std ...
                    ) * obj.values_adj(end).objective.std + obj.values_adj(end).objective.mean;
end