% compare the L2 error in moment approximation computed through odd and
% characteristic boundary conditions. 

% test case could either be :
% 1. gaussian_collision
% 2. unsteady_lid_driven_cavity
% 3. lid_driven_cavity
% 4. heated_cavity

% bc_type could be 
% 1. odd
% 2. char

clear all;

test_case = 'gaussian_collision';
nx_values1 = 10:5:35;
error_test1 = grid_refine(test_case,nx_values1);
convg1 = -log(error_test1(1:end-1)'./error_test1(2:end)')./log(nx_values1(1:end-1)./nx_values1(2:end));

test_case = 'unsteady_lid_driven_cavity';
nx_values2 = 10:5:40;
error_test2 = grid_refine(test_case,nx_values2);
convg2 = -log(error_test2(1:end-1)'./error_test2(2:end)')./log(nx_values2(1:end-1)./nx_values2(2:end));

test_case = 'lid_driven_cavity';
nx_values3 = 10:5:40;
error_test3 = grid_refine(test_case,nx_values3);
convg3 = -log(error_test3(1:end-1)'./error_test3(2:end)')./log(nx_values3(1:end-1)./nx_values3(2:end));

test_case = 'heated_cavity';
nx_values4 = [10:5:55];
error_test4 = grid_refine(test_case,nx_values4);
convg4 = -log(error_test4(1:end-1)'./error_test4(2:end)')./log(nx_values4(1:end-1)./nx_values4(2:end));

plot(1:length(error_test1)-1,convg1,'-bo',...
       1:length(error_test2)-1,convg2,'-r*',...
       1:length(error_test3)-1,convg3,'-ms',...
       1:length(error_test4)-1,convg4,'-kd',...
       'markersize',10,'LineWidth',5,...
       'MarkerFaceColor','b','MarkerFaceColor','r',...
       'MarkerFaceColor','m','MarkerFaceColor','k');
   
%set(gca,'YScale','log'); 
ylabel('convergence rate');
xlabel('refinement cycle');
legend('Test-1','Test-2','Test-3','Test-4','Location','best');   
grid on;
set(gca, 'FontSize', 20);

base_output = strcat('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/EM_grid_refine');
print(base_output,'-depsc');




function f = grid_refine(test_case,nx_values)

base_output = strcat('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/',test_case,'/');

%nx_values = 10:5:55;

%% ids of macroscopic quantities
filename_dvm = strcat('../DVM/',test_case,'/result_n100_DVM_20.txt');
result_mom = cell(length(nx_values),1);% result from moments
result_dvm = dlmread(filename_dvm,'\t');

X = result_dvm(1,:);
Y = result_dvm(2,:);
grid_points = 101;

X = reshape(X,grid_points,grid_points);
Y = reshape(Y,grid_points,grid_points);

delta_x = X(2,1)-X(1,1);
[~,PX] = sbp_collocated_2(grid_points-1,delta_x);

for i = 1 : length(nx_values)
    filename_moments = strcat('../results/',test_case,'_','char','/',...
                        'result_n',num2str(nx_values(i)),'_M12','.txt');
   % filename_moments = strcat('../results/heated_cavity_char/result_M12.txt');
    result_mom{i} = dlmread(filename_moments,'\t');
    result_mom{i} = interpolate_2D(result_mom{i},nx_values(i),X,Y);   
end

[error] = compute_error(result_dvm,result_mom,PX);
error_entropy = sqrt(dot(error,error,2));

f = error_entropy;

end

function f = interpolate_2D(result,nx,Xq,Yq)
    
    num_moms = size(result,1);
    num_grid = length(Xq(:));
    f = zeros(num_moms,num_grid);
    
    X_mom = result(1,:);
    Y_mom = result(2,:);
    
    f(1,:) = Xq(:);
    f(2,:) = Yq(:);
    
    X_mom = reshape(X_mom,nx + 1,nx + 1);
    Y_mom = reshape(Y_mom,nx + 1,nx + 1);
    
    for j = 3 : size(result,1)
        temp = reshape(result(j,:),size(X_mom));
        temp = interp2(X_mom',Y_mom',temp',Xq',Yq');
        temp = temp';
        f(j,:) = temp(:);
    end
    
end

