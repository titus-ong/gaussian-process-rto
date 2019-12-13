excel = "\System Files\Operating region plots\case study 80(10%) run2.xlsx";
rows = 15;  % No. of intervals of reboiler_duty
cols = 17;  % No. of intervals of solvent flowrate

data_grid = readmatrix(pwd + excel);
x = linspace(60000, 200000, 15);
y = linspace(200, 1000, 17);

% Objective
z1 = plot_columns(data_grid, rows, cols, 36)';
% Min
[~,i] = min(z1(:));
[row,col] = ind2sub(size(z1), i);

% CO2% in clean gas
z2 = plot_columns(data_grid, rows, cols, 7)';

figure
hold on;
surf(x, y, z1)
h = scatter3(x(col),y(row),z1(i),'r', 'filled');

figure
hold on;
contour(x, y, z1, 20)
scatter(x(col),y(row),'rx');
% Draw CO2% contour
[one_perc, h] = contour(x, y, z2, [0.007 0.007]);
clabel(one_perc, h)
[two_perc, h] = contour(x, y, z2, [0.01 0.01]);
clabel(two_perc, h)
[three_perc, h] = contour(x, y, z2, [0.02 0.02]);
clabel(three_perc, h)

function [z] = plot_columns(data_grid, rows, cols, column)
z = zeros(rows, cols);
for i = 1:cols
    z(1:rows, i) = data_grid((i-1)*rows + 1:i*rows, column);
end
end