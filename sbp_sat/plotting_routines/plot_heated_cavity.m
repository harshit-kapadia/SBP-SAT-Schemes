clear all;
M_values = 3:2:13;

%% ids of macroscopic quantities
shift = 2; % shift because of coordinate data
ID_rho = 1 + shift;
ID_ux = 2 + shift;
ID_uy = 3 + shift;
ID_theta = 4 + shift;
ID_sigma_xx = 5 + shift;
ID_sigma_xy = 6 + shift;
ID_sigma_yy = 7 + shift;
ID_qx = 8 + shift;
ID_qy = 9 + shift;

filename_moments = cell(length(M_values));
filename_dvm = '../DVM/heated_cavity/result_DVM_10.txt';
result_mom = cell(length(M_values));% result from moments
result_dvm = dlmread(filename_dvm,'\t');

for i = 1 : length(M_values)
    filename_moments{i} = strcat('../heated_cavity/','result_M',num2str(M_values(i)),'.txt');
    result_mom{i} = dlmread(filename_moments{i},'\t');
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

%% residual variation
filename = '../heated_cavity/residual_M13.txt';
residual = dlmread(filename,'\t');

plot_residual(residual,2);


%% computes error in all macroscopic quantities
function error = compute_error(result_dvm,result_mom,PX)

num_mom = length(result_mom);
num_macros = size(result_dvm,1)-2;

% num of mom times the num of macroscopic quantitites
error = zeros(num_mom,num_macros);

for i = 1 : num_mom
    for j = 1 : num_macros
        temp = result_dvm(j + 2,:) - result_mom{i}(j  + 2, :);
        temp = reshape(temp,size(PX));
        int_x = dot(transpose(temp),transpose(PX * temp),2);  % integral along x.
        int_xy = sum(PX*int_x);  % integral along xy.
        error(i,j) = sqrt(int_xy);
    end    
end

end

function[] = plot_error(M_values,error,fig_id)

min_x_value = min(M_values)-0.25;
max_x_value = max(M_values)+1;

min_y_value = min(error(end,:));
max_y_value = max(error(1,:));

max_y_value = 1.5 * max_y_value;
min_y_value = 0.75 * min_y_value;

figure(fig_id);

[reference_line] = exact_order(M_values(3),error(3,5)+0.001,M_values(5),-1);

loglog(M_values,error(:,1),'-go', ...
       M_values,error(:,2),'--m+', ...
       M_values,error(:,3),'-m*', ...
       M_values,error(:,4),'-r.', ...
       M_values,error(:,5),'-bx', ...
       M_values,error(:,6),'--bs', ...
       M_values,error(:,7),'-.bd', ...
       M_values,error(:,8),'-c^', ...
       M_values,error(:,9),'--c>',...
       reference_line(1:2),reference_line(3:4),'--pk',...
       'markersize',8);

legend('rho','u_x','u_y','\theta','\sigma_{xx}','\sigma_{xy}','\sigma_{yy}','q_x','q_y','reference');
title('L^2 error in different macroscopoic quantities');
xlim([min_x_value,max_x_value]);
ylim([min_y_value,max_y_value]);
h = xlabel('M','FontSize',18);
ylabel('Error');
xticks(M_values);
xt = get(gca, 'YTick');
grid on;
set(gca, 'FontSize', 16);

end

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

