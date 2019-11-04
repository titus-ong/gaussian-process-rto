WO = WilliamsOtto();
rows = 16;
cols = 31;
x1 = linspace(WO.lb(1), WO.ub(1), cols);  % Flowrate of B
x2 = linspace(WO.lb(2), WO.ub(2), rows);  % Temp of reactor

data = zeros(rows, cols);
x_a = zeros(rows, cols);
x_g = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = WO.true_obj_fn([x1(i), x2(j)]);
        [~, temp] = WO.get_output([x1(i), x2(j)]);
        x_a(j, i) = temp.x_a;
        x_g(j, i) = temp.x_g;
    end
end

figure
hold on;
contour(x1, x2, data, 30)
[c, h] = contour(x1, x2, x_a, [0.12 0.12]);
% clabel(c, h);
[c, h] = contour(x1, x2, x_g, [0.08 0.08]);
% clabel(c, h);