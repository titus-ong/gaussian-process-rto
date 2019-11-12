% Run example_WO first to get GP variable in workspace
% WO = WilliamsOtto();
rows = 16;
cols = 31;
x1 = linspace(WO.lb(1), WO.ub(1), cols);  % Flowrate of B
x2 = linspace(WO.lb(2), WO.ub(2), rows);  % Temp of reactor
temp_GP = copy(GP);
model_idx = 20;
temp_GP.model = temp_GP.model(1:model_idx);

data = zeros(rows, cols);
x_a = zeros(rows, cols);
x_g = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = temp_GP.obj_fn([x1(i), x2(j)]);
        a = predict(temp_GP.model(model_idx).x_a_con, ([x1(i), x2(j)]));
        g = predict(temp_GP.model(model_idx).x_g_con, ([x1(i), x2(j)]));
        x_a(j, i) = a;
        x_g(j, i) = g;
    end
end

f = figure;
hold on;
contour(x1, x2, data, 30)
% surf(x1, x2, data)
[c, h] = contour(x1, x2, x_a, [0 0]);
% clabel(c, h);
[c, h] = contour(x1, x2, x_g, [0 0]);
% clabel(c, h);