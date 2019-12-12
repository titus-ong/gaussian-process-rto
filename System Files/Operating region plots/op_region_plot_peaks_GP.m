Pks = PeaksFunction();
rows = 16;
cols = 31;
x1 = linspace(Pks.lb(1), Pks.ub(1), cols);  % Flowrate of B
x2 = linspace(Pks.lb(2), Pks.ub(2), rows);  % Temp of reactor
temp_GP = copy(GP);
model_idx = 1;
temp_GP.model = temp_GP.model(model_idx);

data = zeros(rows, cols);
for i = 1:cols
    for j = 1:rows
        data(j, i) = temp_GP.obj_fn([x1(i), x2(j)]);
    end
end

figure
hold on;
surf(x1, x2, data)
scatter(temp_GP.training_starter(:,1), temp_GP.training_starter(:,2));