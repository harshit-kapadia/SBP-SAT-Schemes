function [] = plot_residual(residual,fig_id)

min_x_value = min(residual(:,1));
max_x_value = max(residual(:,1));

min_y_value = min(residual(:,2));
max_y_value = max(residual(:,2));

max_y_value = 1.5 * max_y_value;
min_y_value = 0.75 * min_y_value;

max_x_value = 1.5 * max_x_value;
min_x_value = 0.75 * min_x_value;


figure(fig_id);

loglog(residual(:,1),residual(:,2),'-r*');
title('residual decay with time');
h = xlabel('time (t)','FontSize',18);
ylabel('Residual');
xlim([min_x_value,max_x_value]);
ylim([min_y_value,max_y_value]);
xt = get(gca, 'YTick');
grid on;
set(gca, 'FontSize', 16);

end