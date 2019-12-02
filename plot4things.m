set(groot,'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

y_delta = GP.delta(:,1)/GP.delta(1,1);
y_rho = GP.rho;
y_obj = GP.fval_true;
y_con = GP.ineq_true;
x = size(GP.rho,1);
x = linspace(0,x-1,x);
%% 
fig1 = figure();
hold on

yyaxis left
rho = plot(x,y_rho);
ylabel('$\rho$','interpreter','latex')
test = isnan(y_rho);
for i=1:size(y_rho,1)
    if test(i) == 1
        k = i;
        while test(k) == 1
            k = k-1;
            val = y_rho(k);
        end 
        scatter(i-1, val)
    end
end
      
yyaxis right
delta = plot(x,y_delta);
ylabel('$\Delta$','interpreter','latex')

legend([rho, delta], {'$\rho$','$\Delta$'})
xlabel('Iteration','interpreter','latex')
fig1_plot = gca;
fig1_plot.YAxis(1).Color = 'k';
fig1_plot.YAxis(2).Color = 'k';
xlim([0 size(GP.rho,1)]-1)
%%
fig2 = figure();
hold on

yyaxis left
plot(x,y_obj);
ylabel('$G_{0}$','interpreter','latex')

yyaxis right
plot(x,y_con);
ylabel('$G_{1}$','interpreter','latex')

legend('$G_{0}$', '$G{1}$')
xlabel('Iteration')
fig2_plot = gca;
fig2_plot.YAxis(1).Color = 'k';
fig2_plot.YAxis(2).Color = 'k';
xlim([0 size(GP.rho,1)]-1)