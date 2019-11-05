WO = WilliamsOttoUnconstrained();
WO.time = 0;
rows = 30;
cols = 31;
x1 = linspace(WO.lb(1), WO.ub(1), cols);  % Flowrate of B
x2 = linspace(WO.lb(2), WO.ub(2), rows);  % Temp of reactor

data = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = WO.true_obj_fn([x1(i), x2(j)]);
    end
end

[~,i] = min(data(:));
[row,col] = ind2sub(size(data), i);

figure
hold on;
contour(x1, x2, data, 30)
scatter(x1(col), x2(row),'rx');