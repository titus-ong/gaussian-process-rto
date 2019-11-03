WO = WilliamsOtto();
rows = 16;
cols = 31;
x1 = linspace(WO.lb(1), WO.ub(1), cols);  % Flowrate of B
x2 = linspace(WO.lb(2), WO.ub(2), rows);  % Temp of reactor

data = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = WO.true_obj_fn([x1(i), x2(j)]);
    end
end

figure
contour(x1, x2, data, 30)