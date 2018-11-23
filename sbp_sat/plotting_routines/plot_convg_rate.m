
function [] = plot_convg_rate(convg_MOdd,convg_MEven,fig_id)

hFig = figure(fig_id);
set(hFig, 'Position', [5 5 800 800]);

xticks = {'\rho','v_1','v_2','\theta','\sigma_{11}','\sigma_{12}','\sigma_{22}','q_1','q_2'};
num_quantities = length(convg_MOdd);
plot(1:num_quantities,convg_MOdd,'-bo',1:num_quantities,convg_MEven,'-r+',...
    'markersize',10);
ylabel('Convergence Rate of e_M^{\alpha}');
grid on;
set(gca,'xtick',[1:num_quantities],'xticklabel',xticks);
legend('Odd M','Even M');
set(gca, 'FontSize', 20);

end
