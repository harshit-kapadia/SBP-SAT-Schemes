function f = plot_discretization_error()

filename_char = '../discretization_error/char_penalty_M3.txt';
filename_odd = '../discretization_error/odd_penalty_M3.txt';

result_char = dlmread(filename_char,'\t');
result_odd = dlmread(filename_odd,'\t');

[reference_line] = exact_order(10^(-2),0.2*10^(-3),0.5*10^(-1),2);

hFig = figure(1);
set(hFig, 'Position', [0 0 800 800]);

loglog(result_odd(1,:),result_odd(2,:),'-bs',...
       result_char(1,:),result_char(2,:),'-rs',...
       reference_line(1:2),reference_line(3:4),'--kx',...
       'markersize',10,...
       'LineWidth',5,...
       'MarkerFaceColor','b',...
       'MarkerFaceColor','r');
grid on;
legend('Odd Penalty','Characteristic Penalty',...
       'Reference (slope=2)','Location','southeast');
xlabel('h');
ylabel('$\tilde{E}_h$','Interpreter','latex');
title('Discretization Error');
set(gca,'FontSize',20);

output_filename = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/discretisation_error/error';
print(output_filename,'-depsc');

end