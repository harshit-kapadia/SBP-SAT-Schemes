% we compare the results from the characteristic penalty matrix and the odd
% penalty matrix
% test case could either be :
% 1. gaussian_collision
% 2. unsteady_lid_driven_cavity
% 3. lid_driven_cavity
% 4. heated_cavity

% bc_type could be 
% 1. odd
% 2. char

function [] = compare_EM(test_case)

base_output = strcat('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/',test_case,'/');

M_values = 3:1:13;

%% ids of macroscopic quantities
filename_dvm = strcat('../DVM/',test_case,'/result_n100_DVM_20.txt');
result_mom = cell(length(M_values));% result from moments
result_mom_odd = cell(length(M_values));% result from moments
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
grid_points = 101;
X = reshape(X,grid_points,grid_points);
delta_x = X(2,1)-X(1,1);
[~,PX] = sbp_collocated_2(grid_points-1,delta_x);

[error] = compute_error(result_dvm,result_mom,PX);
[error_odd] = compute_error(result_dvm,result_mom_odd,PX);

error_entropy = sqrt(dot(error,error,2));
error_entropy_odd = sqrt(dot(error_odd,error_odd,2));

slope = 1;
plot_error_entropy(M_values,error_entropy,error_entropy_odd,slope);

filename = strcat(base_output,'EM');
print(filename,'-depsc');

end

