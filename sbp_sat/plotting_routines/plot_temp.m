M_values = 4:2:13;

    filename = strcat('../DVM/gaussian_collision/result_n100_DVM_20.txt');
    result_dvm = dlmread(filename,'\t');
    
for i = 1 : length(M_values)
    filename_moments = strcat('../results/gaussian_collision_','char','/',...
                        'result_M',num2str(M_values(i)),'.txt');
    result_mom{i} = dlmread(filename_moments,'\t');

    X = reshape(result_mom{1}(1,:),101,101);
    Y = reshape(result_mom{1}(2,:),101,101);
    rho = reshape(result_mom{i}(3,:),101,101);
    rho_DVM = reshape(result_dvm(3,:),101,101);
    figure(1);
    %contourf(X,Y,rho);
    %colorbar;
    %title(strcat('M=',num2str(M_values(i))));
    plot(Y(51,:),rho(51,:),'-',Y(51,:),rho_DVM(51,:),'--','LineWidth',3);
    hold on;
end



