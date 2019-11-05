% Run example script first to get GP variable in workspace
rows = 16;
cols = 31;
x1 = linspace(hysys.lb(1), hysys.ub(1), cols);
x2 = linspace(hysys.lb(2), hysys.ub(2), rows);
temp_GP = copy(GP);
model_idx = 11;
temp_GP.model = temp_GP.model(1:model_idx);

data = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = temp_GP.obj_fn([x1(i), x2(j)]);
    end
end

f = figure;
hold on;
surf(x1, x2, data)
scatter(temp_GP.training_input(1:6,1), temp_GP.training_input(1:6,2));