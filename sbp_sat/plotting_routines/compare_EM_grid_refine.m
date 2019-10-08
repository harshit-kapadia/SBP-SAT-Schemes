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

function f = compare_EM_grid_refine(test_case)

base_output = strcat('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/',test_case,'/');

%nx_values = 10:5:55;
nx_values = [10:5:55,95,100];

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
% loglog(error_entropy,'-o');
% slope = 1;
% plot_error_entropy(M_values,error_entropy,error_entropy_odd,slope);
% 
% filename = strcat(base_output,'EM');
% print(filename,'-depsc');

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

