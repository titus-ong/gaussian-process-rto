Pks = PeaksFunction();
rows = 16;
cols = 31;
x1 = linspace(Pks.lb(1), Pks.ub(1), cols);
x2 = linspace(Pks.lb(2), Pks.ub(2), rows); 

data = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = Pks.true_obj_fn([x1(i), x2(j)]);
    end
end

figure
hold on;
contour(x1, x2, data, 30)