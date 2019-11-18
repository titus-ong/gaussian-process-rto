WO = WilliamsOtto();
rows = 16;
cols = 31;
x1 = linspace(WO.lb(1), WO.ub(1), cols);  % Flowrate of B
x2 = linspace(WO.lb(2), WO.ub(2), rows);  % Temp of reactor

% Pre-allocation
data = zeros(rows, cols);
for i = 1:length(WO.constraints_ineq)
    con_ineq.(WO.constraints_ineq{i}) = zeros(rows, cols);
end
for i = 1:length(WO.constraints_eq)
    con_eq.(WO.constraints_eq{i}) = zeros(rows, cols);
end

for i = 1:cols
    for j = 1:rows
        [data(j, i), val_ineq, val_eq] = WO.get_output([x1(i), x2(j)]);
        for k = 1:length(WO.constraints_ineq)
            con_ineq.(WO.constraints_ineq{k})(j, i) = val_ineq(k);
        end
        for m = 1:length(WO.constraints_eq)
            con_eq.(WO.constraints_eq{m})(j, i) = val_eq(m);
        end
    end
end

% Plots
figure
hold on;
contour(x1, x2, data, 30)

for i = 1:length(WO.constraints_ineq)
    contour(x1, x2, con_ineq.(WO.constraints_ineq{i}), [0 0], 'k-');
end
for i = 1:length(WO.constraints_eq)
    contour(x1, x2, con_eq.(WO.constraints_eq{i}), [0 0], 'm-');
end

% Optima
[~,i] = min(data(:));
[row,col] = ind2sub(size(data), i);
scatter(x1(col), x2(row), 'rx');