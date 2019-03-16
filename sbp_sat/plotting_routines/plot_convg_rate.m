
function [] = plot_convg_rate(convg_MOdd,convg_MEven,fig_id)

hFig = figure(fig_id);
set(hFig, 'Position', [5 5 800 800]);

xticks = {'\rho','v_1','v_2','\theta','\sigma_{11}','\sigma_{12}','\sigma_{22}','q_1','q_2'};
num_quantities = length(convg_MOdd);
plot(1:num_quantities,convg_MOdd,'-bs',1:num_quantities,convg_MEven,'-rs',...
    'markersize',10,'LineWidth',5,'MarkerFaceColor','b','MarkerFaceColor','r');
ylabel('Convergence Rate of e_M^{\alpha}');
grid on;
set(gca,'xtick',[1:num_quantities],'xticklabel',...
    {'$\tilde{\rho}$','$\tilde{v}_1$','$\tilde{v}_2$','$\tilde{\theta}$','$\tilde{\sigma}_{11}$','$\tilde{\sigma}_{12}$','$\tilde{\sigma}_{22}$','$\tilde{q}_1$','$\tilde{q}_2$'},...
    'TickLabelInterpreter','latex');
legend('Odd M','Even M');
set(gca, 'FontSize', 20);

end