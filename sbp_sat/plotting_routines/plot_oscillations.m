clear all;

kinetic_T03 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_T03.txt');
kinetic_M5_T03 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_M5_T03.txt');
kinetic_T02 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_T02.txt');

odd_T03 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/odd_T03.txt');
odd_T02 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/odd_T02.txt');

%% plot at t = 0.2
hFig = figure(1);
set(hFig, 'Position', [0 0 800 800]);

plot(kinetic_T03(1,:),kinetic_T03(2,:),'-bo',...
    odd_T03(1,:),odd_T03(2,:),'-ro','markersize',10,...
    'LineWidth',3);

xlabel('x_1');
title('$\tilde{\rho}(x_1,x_2=0.5,t=0.3)$','Interpreter','latex');
legend('Weak kinetic boundary (M=3)','Odd penalty matrix (M=3)','Location','southeast');
grid on;
xt = get(gca, 'YTick');
ylim([-0.1,0.15]);
set(gca, 'FontSize', 20);
filename = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_odd_T03';
print(filename,'-depsc');

%% plot at t=0.3
hFig = figure(2);
set(hFig, 'Position', [0 0 800 800]);
plot(kinetic_T02(1,:),kinetic_T02(2,:),'-bo',...
    odd_T02(1,:),odd_T02(2,:),'-ro','markersize',10,...
    'LineWidth',3);

xlabel('x_1');
title('$\tilde{\rho}(x_1,x_2=0.5,t=0.2)$','Interpreter','latex');
legend('Weak kinetic boundary (M=3)','Odd penalty matrix (M=3)','Location','best');
grid on;
xt = get(gca, 'YTick');
ylim([-0.01,0.2]);
set(gca, 'FontSize', 20);

filename = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_odd_T02';
print(filename,'-depsc');

%% plot for M = 5
hFig = figure(3);
set(hFig, 'Position', [0 0 800 800]);
plot(kinetic_M5_T03(1,:),kinetic_M5_T03(2,:),'-bo',...
     'markersize',10,...
    'LineWidth',3);

xlabel('x_1');
title('$\tilde{\rho}(x_1,x_2=0.5,t=0.3)$','Interpreter','latex');
legend('Weak kinetic boundary (M=5)','Location','best');
grid on;
xt = get(gca, 'YTick');
ylim([-0.03,0.11]);
set(gca, 'FontSize', 20);

filename = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_odd_M5_T03';
print(filename,'-depsc');