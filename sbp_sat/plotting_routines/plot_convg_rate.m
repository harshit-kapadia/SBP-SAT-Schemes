
function [] = plot_convg_rate(convg_MOdd,convg_MEven,fig_id)

hFig = figure(fig_id);
set(hFig, 'Position', [5 5 800 800]);

max_y_value = 1.1 * max([max(convg_MOdd),max(convg_MEven)]);
min_y_value = 0.85 * min([min(convg_MOdd),min(convg_MEven)]);

xticks = {'\rho','v_1','v_2','\theta',...
          '\sigma_{11}','\sigma_{12}','\sigma_{22}','q_1','q_2'};
      
num_quantities = length(convg_MOdd);

plot(1:num_quantities,convg_MOdd,'-bs',...
     1:num_quantities,convg_MEven,'-rs',...
    'markersize',10,'LineWidth',5,'MarkerFaceColor','b','MarkerFaceColor','r');

ylabel('Convergence Rate of e_{\alpha,M}');
ylim([0.8,1.8]);
grid on;
set(gca,'xtick',[1:num_quantities],'xticklabel',...
    {'$\rho$','$v_1$',...
     '$v_2$','$\theta$',...
     '$\sigma_{11}$','$\sigma_{12}$',...
     '$\sigma_{22}$','$q_1$','$q_2$'},...
     'TickLabelInterpreter','latex');
legend('Odd M','Even M');
set(gca, 'FontSize', 20);

end
