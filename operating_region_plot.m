excel = "\op_region.xlsx";
rows = 7;  % No. of intervals of reboiler_duty
cols = 7;  % No. of intervals of solvent_flowrate

data_grid = readmatrix(pwd + excel);

[x1, y1, z1] = plot_columns(data_grid, rows, cols, 25);
[x2, y2, z2] = plot_columns(data_grid, rows, cols, 8);

figure
hold on;
contour(x1, y1, z1)  % Use surf to see 3D plot
contour(x2, y2, z2, [0.003 0.003])  % To draw 0.003 CO2% contour

function [x, y, z] = plot_columns(data_grid, rows, cols, column)
z = zeros(rows, cols);
for i = 1:cols
    z(1:rows, i) = data_grid((i-1)*rows + 1:i*rows, column);
end
x = 150:25:300;
y = 108000:2000:120000;
end