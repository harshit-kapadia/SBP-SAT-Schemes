function[fig] = plot_error(M_values,error,fig_id)

min_x_value = min(M_values)-0.25;
max_x_value = max(M_values)+0.5;

min_y_value = min(error(end,:));
max_y_value = max(error(1,:));

max_y_value = 1.1 * max_y_value;
min_y_value = 0.85 * min_y_value;

hFig = figure(fig_id);
set(hFig, 'Position', [0 0 800 800]);

fig = loglog(M_values,error(:,1),'-ko', ...
       M_values,error(:,2),'--m+', ...
       M_values,error(:,3),'-m*', ...
       M_values,error(:,4),'-r.', ...
       M_values,error(:,5),'-bx', ...
       M_values,error(:,6),'--bs', ...
       M_values,error(:,7),'-.bd', ...
       M_values,error(:,8),'-c^', ...
       M_values,error(:,9),'--c>',...
       'markersize',10,...
       'LineWidth',3);
    
lg = legend({'$\tilde{\rho}$','$\tilde{v}_1$','$\tilde{v}_2$','$\tilde{\theta}$','$\tilde{\sigma}_{11}$',...
        '$\tilde{\sigma}_{12}$','$\tilde{\sigma}_{22}$','$\tilde{q}_1$','$\tilde{q}_2$'},...
        'Interpreter','latex','Location','eastoutside');
    
lg.FontSize = 20;
    
title('L^2 error in macroscopoic quantities');
xlim([min_x_value,max_x_value]);
ylim([min_y_value,max_y_value]);

h = xlabel('M','FontSize',20);
ylabel('e_M^{\alpha}');
xticks(min(M_values):2:max(M_values));
xt = get(gca, 'YTick');
grid on;
set(gca, 'FontSize', 20);

end