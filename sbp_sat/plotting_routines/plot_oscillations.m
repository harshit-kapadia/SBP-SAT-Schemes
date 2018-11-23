clear all;

kinetic_T03 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_T03.txt');
kinetic_T02 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_T02.txt');

odd_T03 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/odd_T03.txt');
odd_T02 = dlmread('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/odd_T02.txt');

hFig = figure(1);
set(hFig, 'Position', [0 0 800 800]);

plot(kinetic_T03(1,:),kinetic_T03(2,:),'-bo',...
    odd_T03(1,:),odd_T03(2,:),'-ro');

xlabel('x_1');
title('$\tilde{\rho}(x_1,x_2=0.5,t=0.3)$','Interpreter','latex');
legend('Kinetic weak boundary','Odd penalty matrix','Location','best');
grid on;
xt = get(gca, 'YTick');
ylim([-0.1,0.4]);
set(gca, 'FontSize', 20);

hFig = figure(2);
set(hFig, 'Position', [0 0 800 800]);
plot(kinetic_T02(1,:),kinetic_T02(2,:),'-bo',...
    odd_T02(1,:),odd_T02(2,:),'-ro');
xlabel('x_1');
title('$\tilde{\rho}(x_1,x_2=0.5,t=0.2)$','Interpreter','latex');
legend('Kinetic weak boundary','Odd penalty matrix','Location','best');
grid on;
xt = get(gca, 'YTick');
ylim([-0.1,0.4]);
set(gca, 'FontSize', 20);
