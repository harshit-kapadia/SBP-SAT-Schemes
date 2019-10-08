% Plot the L2 error in moment approximation, ||f-f_M||, corresponding to 
% odd and characteristic boundary implementation. 
function[] = plot_error_entropy(M_values,error_entropy,error_entropy_odd,slope)

min_x_value = min(M_values)-0.25;
max_x_value = max(M_values)+0.5;

min_y_value = min(error_entropy(end));
max_y_value = max(error_entropy(1));

max_y_value = 1.1 * max_y_value;
min_y_value = 0.85 * min_y_value;

hFig = figure(1);
set(hFig, 'Position',[0 0 800 800]);

[reference_line] = exact_order(M_values(3),...
                            error_entropy(3)+0.01,M_values(9),-slope);

loglog(M_values,error_entropy,'-ro',...
    M_values,error_entropy_odd,'-b*',...
    reference_line(1:2),reference_line(3:4),'--pk',...
       'markersize',10,'LineWidth',5);

legend('Characteristic Penalty','Odd Penalty','reference (slope=-1)');
xlim([min_x_value,max_x_value]);
ylim([min_y_value,max_y_value]);
h = xlabel('M','FontSize',20);
ylabel('E_{I,M}');
xticks(min(M_values):2:max(M_values));
xt = get(gca, 'YTick');
grid on;
set(gca, 'FontSize', 20);

end
