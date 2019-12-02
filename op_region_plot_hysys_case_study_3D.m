excel = "\case study 3D run2.xlsx";
rows = 8;  % No. of intervals of reboiler_duty
cols = 9;  % No. of intervals of solvent flowrate
heights = 9;

data_grid = readmatrix(pwd + excel);
x1 = linspace(60000, 130000, cols);
x2 = linspace(200, 600, rows);
x3 = linspace(40, 120, heights);

% Objective
z1 = plot_columns(data_grid, rows, cols, heights, 31);
% Min
[~,i] = min(z1);
[row,col,height] = ind2sub(size(z1), i);

% CO2% in clean gas
z2 = plot_columns(data_grid, rows, cols, heights, 8);

%use gp to plot
m = mean(data_grid);
sd = std(data_grid);
X = (data_grid-m)./sd;
X = X(:,1:3);
obj_GP = fitrgp(X,data_grid(:,31),'ConstantSigma',true,'Sigma',1e-10);
constraint_GP = fitrgp(X,data_grid(:,8),'ConstantSigma',true,'Sigma',1e-10);
rows = 20;
cols = 20;
heights = 20;

x1 = linspace(60000, 130000, cols);
x2 = linspace(200, 600, rows);
x3 = linspace(40, 120, heights);

x1 = (x1-mean(x1))./std(x1);
x2 = (x2-mean(x2))./std(x2);
x3 = (x3-mean(x3))./std(x3);

obj = zeros(rows, cols, heights);
constraint = zeros(rows, cols, heights);

for i = 1:cols
    for j = 1:rows
        for k = 1:heights
            obj(j, i, k) = predict(obj_GP,[x1(i), x2(j), x3(k)]);
            
            constraint(j, i, k) = predict(constraint_GP, [x1(i), x2(j), x3(k)]);

        end
    end
end

obj_withinconstraint = obj;

for i = 1:cols
    for j = 1:rows
        for k = 1:heights
            if constraint(j, i, k) < 0.01
            else
            obj_withinconstraint(j, i, k) = NaN;
            end
        end
    end
end

sliceomatic(obj);
sliceomatic(constraint);

% thing=gca;
% copyobj(thing.Children,paste);

[~, idx] = min(obj_withinconstraint(:));
[n, m, t] = ind2sub(size(obj_withinconstraint),idx);
scatter3(n,m,t,40,'r','h');


% xslice = [100000 120000];                               % define the cross sections to view
% yslice = [500];
% zslice = [70 90];
% 
% slice(x1, x2, x3, data, xslice, yslice, zslice)    % display the slices

% xlim([80000 120000])
% ylim([200 500])
% zlim([70 100])

% cb = colorbar;                                  % create and label the colorbar

% surf(x, y, z1)
% h = scatter3(x(col),y(row),z1(i),'r', 'filled');
% 
% figure
% hold on;
% contour(x, y, z1, 20)
% scatter(x(col),y(row),'rx');
% % Draw CO2% contour
% [one_perc, h] = contour(x, y, z2, [0.007 0.007]);
% clabel(one_perc, h)
% [two_perc, h] = contour(x, y, z2, [0.01 0.01]);
% clabel(two_perc, h)
% [three_perc, h] = contour(x, y, z2, [0.02 0.02]);
% clabel(three_perc, h)

function [z] = plot_columns(data_grid, rows, cols, height, column)
z = zeros(rows, cols, height);
    for i = 1:rows
        for j = 1:cols
            z(i, j, 1:height) = data_grid((i-1)*height*cols + (j-1)*cols + 1:j*cols +(i-1)*height*cols, column);
        end
    end
end