% does the plotting for L2 error in macroscopic quantities

% test case could either be :
% 1. gaussian_collision
% 2. unsteady_lid_driven_cavity
% 3. lid_driven_cavity
% 4. heated_cavity

function [] = plot_L2_error_macro(test_case)

base_output = strcat('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/',...
              test_case,'/');

M_values = 3:1:13;

filename_moments = cell(length(M_values));
filename_dvm = strcat('../DVM/',test_case,'/result_n100_DVM_20.txt');
result_mom = cell(length(M_values));% result from moments
result_dvm = dlmread(filename_dvm,'\t');

for i = 1 : length(M_values)
    filename_moments = strcat('../results/',test_case,'_','char','/',...
                        'result_M',num2str(M_values(i)),'.txt');
    result_mom{i} = dlmread(filename_moments,'\t');
    
    filename_moments = strcat('../results/',test_case,'_','odd','/',...
                        'result_M',num2str(M_values(i)),'.txt');
    result_mom_odd{i} = dlmread(filename_moments,'\t');
end

X = result_dvm(1,:);
Y = result_dvm(2,:);

grid_points = 101;

X = reshape(X,grid_points,grid_points);
Y = reshape(Y,grid_points,grid_points);

%% error computation
delta_x = X(2,1)-X(1,1);
[~,PX] = sbp_collocated_2(grid_points-1,delta_x);
[error] = compute_error(result_dvm,result_mom,PX);
plot_error(M_values,error,1);

filename = strcat(base_output,'macro_','char');
print(filename,'-depsc');

[error_Odd] = compute_error(result_dvm,result_mom_odd,PX);
plot_error(M_values,error_Odd,2);

filename = strcat(base_output,'macro_','odd');
print(filename,'-depsc');

%% convergence rate in macroscopic quantities
[convg_MOdd,convg_MEven] = convg_rate_macro(error,M_values); % convergence rate in MOdd and in MEven

plot_convg_rate(convg_MOdd,convg_MEven,3);

filename = strcat(base_output,'convg_rate_','char');
print(filename,'-depsc');

[convg_MOdd_penalty_Odd,convg_MEven_penalty_Odd] = convg_rate_macro(error_Odd,M_values); % convergence rate in MOdd and in MEven

plot_convg_rate(convg_MOdd_penalty_Odd,convg_MEven_penalty_Odd,4);

filename = strcat(base_output,'convg_rate_','odd');
print(filename,'-depsc');

%% making pie charts
per_diff_Odd = 100 * (convg_MOdd - convg_MOdd_penalty_Odd)./convg_MOdd_penalty_Odd;
per_diff_Even = 100 * (convg_MEven - convg_MEven_penalty_Odd)./convg_MEven_penalty_Odd;

per_diff = zeros(length(per_diff_Odd),2);
per_diff(:,1) = per_diff_Odd;
per_diff(:,2) = per_diff_Even;

hFig = figure(5);
set(hFig, 'Position', [5 5 800 800]);
b = bar(per_diff,'BarWidth',1);
grid on;
legend('Odd M','Even M');
set(gca,'xtick',1:9,'xticklabel',...
    {'$\tilde{\rho}$','$\tilde{v}_1$','$\tilde{v}_2$','$\tilde{\theta}$','$\tilde{\sigma}_{11}$','$\tilde{\sigma}_{12}$','$\tilde{\sigma}_{22}$','$\tilde{q}_1$','$\tilde{q}_2$'},...
    'TickLabelInterpreter','latex');
set(gca, 'FontSize', 20);
title('% Difference in convergence rates');

filename = strcat(base_output,'per_diff');
print(filename,'-depsc');

end

