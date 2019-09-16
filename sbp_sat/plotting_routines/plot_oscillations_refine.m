clear all;

nx = [40,60,80,100];
M = [3,4,5,6];
filename = cell(length(nx),1);

for i = 1 : length(nx)
    filename{i} = strcat('../results/gaussian_collision_kinetic/result_char_M',...
                         num2str(M(i)),'_n',num2str(nx(i)),'.txt');
    result{i} = dlmread(filename{i});
end

hFig = figure(1);
set(hFig, 'Position', [0 0 800 800]);

plot(result{1}(1,:),result{1}(2,:),'-bo',...
     result{2}(1,:),result{2}(2,:),'-ro',...
     result{3}(1,:),result{3}(2,:),'-go',...
     result{4}(1,:),result{4}(2,:),'-ko',...
     'markersize',10,...
    'LineWidth',3);

xlabel('x_1');
title('$\tilde{\rho}(x_1,x_2=0.5,t=0.3)$','Interpreter','latex');
legend('M=3','M=4','M=5','M=6','M=7','Location','southeast');
grid on;
xt = get(gca, 'YTick');
%ylim([-0.1,0.15]);
set(gca, 'FontSize', 20);
filename = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_odd_T03';
print(filename,'-depsc');