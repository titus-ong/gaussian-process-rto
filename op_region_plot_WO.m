WO = WilliamsOtto();
% WO.time = 1000;
rows = 16;
cols = 31;
x1 = linspace(WO.lb(1), WO.ub(1), cols);  % Flowrate of B
x2 = linspace(WO.lb(2), WO.ub(2), rows);  % Temp of reactor

data = zeros(rows, cols);
x_a = zeros(rows, cols);
x_g = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        [data(j, i), con_ineq, con_eq] = WO.get_output([x1(i), x2(j)]);
        x_a(j, i) = con_ineq(1);
        x_g(j, i) = con_ineq(2);
    end
end

figure
hold on;
contour(x1, x2, data, 30)
[c, h] = contour(x1, x2, x_a, [0 0]);
% clabel(c, h);
[c, h] = contour(x1, x2, x_g, [0 0]);
% clabel(c, h);