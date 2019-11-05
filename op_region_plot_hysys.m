excel = "\fast run results - Copy.xlsx";
rows = 21;  % No. of intervals of reboiler_duty
cols = 21;  % No. of intervals of sour_gas

data_grid = readmatrix(pwd + excel);
x = linspace(40, 300, 21);
y = linspace(80000, 200000, 21);

% Objective
z1 = plot_columns(data_grid, rows, cols, 27);
% Min
[~,i] = min(z1(:));
[row,col] = ind2sub(size(z1), i);

% CO2% in clean gas
z2 = plot_columns(data_grid, rows, cols, 7);

figure
hold on;
contour(x, y, z1, 20)
scatter(x(col),y(row),'rx');
% Draw CO2% contour
[one_perc, h] = contour(x, y, z2, [0.01 0.01]);
clabel(one_perc, h)
[two_perc, h] = contour(x, y, z2, [0.02 0.02]);
clabel(two_perc, h)
[three_perc, h] = contour(x, y, z2, [0.03 0.03]);
clabel(three_perc, h)

figure
hold on;
surf(x, y, z1)
h = scatter3(x(col),y(row),z1(i),'r', 'filled');

function [z] = plot_columns(data_grid, rows, cols, column)
z = zeros(rows, cols);
for i = 1:cols
    z(1:rows, i) = data_grid((i-1)*rows + 1:i*rows, column);
end
end