clear all;
shift = 2; % shift because of coordinate data
ID_ux = 2 + shift;
ID_theta = 4 + shift;
%% lid-driven cavity solution

filename_moments = strcat('../lid_driven_cavity/','result_M',num2str(12),'.txt');
result_char = dlmread(filename_moments,'\t');

filename_moments = strcat('../lid_driven_cavity_odd/','result_M',num2str(12),'.txt');
result_odd = dlmread(filename_moments,'\t');

filename_dvm = '../DVM/lid_driven_cavity/result_DVM_20.txt';
result_dvm = dlmread(filename_dvm,'\t');

grid_points = 101;

grid_x = linspace(0,1,grid_points);

mid_point = find(grid_x == 0.5);

ux_lid_char = reshape(result_char(ID_ux,:),grid_points,grid_points);
ux_lid_odd = reshape(result_odd(ID_ux,:),grid_points,grid_points);
ux_lid_dvm = reshape(result_dvm(ID_ux,:),grid_points,grid_points);
ux_wall = ones(length(grid_x),1);

[X,Y] = meshgrid(grid_x,grid_x);

figure(1);
plot(grid_x,ux_lid_dvm(mid_point,:),'--go',...
     grid_x,ux_lid_char(mid_point,:),'--r*',...
     grid_x,ux_lid_odd(mid_point,:),'--b',...
     grid_x,ux_wall,'--k',...
     'LineWidth',4);
title('v_1 along a cross-section'); 
h = xlabel('x_2','FontSize',18);
ylabel('v_1');
ylim([-0.2,1.1]);
grid on;
legend('DVM','characteristic','odd','wall velocity (v_1^{in})','Location','best');
set(gca, 'FontSize', 16);

%% results for the heated cavity


filename_moments = strcat('../heated_cavity/','result_M',num2str(13),'.txt');
result_char = dlmread(filename_moments,'\t');

filename_moments = strcat('../heated_cavity_odd/','result_M',num2str(13),'.txt');
result_odd = dlmread(filename_moments,'\t');

filename_dvm = '../DVM/heated_cavity/result_DVM_20.txt';
result_dvm = dlmread(filename_dvm,'\t');

grid_points = 101;

grid_x = linspace(0,1,grid_points);

mid_point = find(grid_x == 0.5);

theta_lid_char = reshape(result_char(ID_theta,:),grid_points,grid_points);
theta_lid_odd = reshape(result_odd(ID_theta,:),grid_points,grid_points);
theta_lid_dvm = reshape(result_dvm(ID_theta,:),grid_points,grid_points);
theta_wall = ones(length(grid_x),1);

[X,Y] = meshgrid(grid_x,grid_x);

figure(2);
plot(grid_x,theta_lid_dvm(mid_point,:),'--go',...
     grid_x,theta_lid_char(mid_point,:),'--r*',...
     grid_x,theta_lid_odd(mid_point,:),'--b',...
     grid_x,theta_wall,'--k',...
     'LineWidth',4);
title('\theta along a cross-section'); 
h = xlabel('x_2','FontSize',18);
ylabel('v_1');
ylim([0,1.1]);
grid on;
legend('DVM','characteristic','odd','wall temperature (\theta^{in})','Location','best');
set(gca, 'FontSize', 16);
