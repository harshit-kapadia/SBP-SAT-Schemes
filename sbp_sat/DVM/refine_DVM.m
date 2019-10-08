% refinement study for the DVM method

function error = refine_DVM(test_case)

nx_values = 10:15:100;
nc_values = [5,7,10,12,15,17,20];

result_dvm = cell(length(nc_values),1);

for i = 1 : length(nx_values)
    filename_dvm = strcat('../DVM/',test_case,'/result_n',num2str(nx_values(i)),...
                          '_DVM_',num2str(nc_values(i)),'.txt');
    result_dvm{i} = dlmread(filename_dvm,'\t');
end

error = error_DVM(result_dvm);

fig = plot(1:length(nx_values)-1,error(:,1),'-ko', ...
           1:length(nx_values)-1,error(:,2),'--m+', ...
           1:length(nx_values)-1,error(:,3),'-m*', ...
           1:length(nx_values)-1,error(:,4),'-r.', ...
           1:length(nx_values)-1,error(:,5),'-bx', ...
           1:length(nx_values)-1,error(:,6),'--bs', ...
       1:length(nx_values)-1,error(:,7),'-.bd', ...
       1:length(nx_values)-1,error(:,8),'-c^', ...
       1:length(nx_values)-1,error(:,9),'--c>',...
       'markersize',10,...
       'LineWidth',3);
    
lg = legend({'$\rho$','$v_1$','$v_2$','$\theta$','$\sigma_{11}$',...
        '$\sigma_{12}$','$\sigma_{22}$','$q_1$','$q_2$'},...
        'Interpreter','latex','Location','eastoutside');
    
lg.FontSize = 20;
    
title('L^2 error in macroscopoic quantities');
xlabel('refinement cycle')
grid on;
set(gca, 'FontSize', 20);
set(gca, 'YScale', 'log');

base_output = strcat('/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/refine_DVM');
filename = strcat(base_output,'EM');
print(filename,'-depsc');

end

function error = error_DVM(result)

steps = length(result);         % total refinements steps
num_moments = size(result{1},1);    % total number of moments
error = zeros(steps-1,num_moments-2);

for i = 1 : (steps-1)
   
    X1 = result{i}(1,:);
    Y1 = result{i}(2,:);
    nx1 = sqrt(length(X1(:))) - 1;
    
    X2 = result{i + 1}(1,:);
    Y2 = result{i + 1}(2,:);
    nx2 = sqrt(length(X2(:))) - 1;
    
    X1 = reshape(X1,[nx1 + 1,nx1 + 1]);
    Y1 = reshape(Y1,[nx1 + 1,nx1 + 1]);
    
    X2 = reshape(X2,[nx2 + 1,nx2 + 1]);
    Y2 = reshape(Y2,[nx2 + 1,nx2 + 1]);
    
    for j = 3 : num_moments
        temp1 = result{i}(j,:);
        temp1 = reshape(temp1,nx1 + 1, nx1 + 1);
        temp1 = interp2(X1',Y1',temp1',X2',Y2');
        temp1 = temp1';
        
        temp2 = result{i + 1}(j,:);
        temp2 = reshape(temp2,nx2 + 1, nx2 + 1);
        
        delta_x2 = X2(2,1)-X2(1,1);
        [~,PX2] = sbp_collocated_2(nx2,delta_x2);
        
        temp = temp1-temp2;
        int_x = dot(transpose(temp),transpose(PX2 * temp),2);  % integral along x.
        int_xy = sum(PX2*int_x);  % integral along xy.
        error(i,j-2) = sqrt(int_xy)/max(max(abs(temp1)));
    end
    
end

end