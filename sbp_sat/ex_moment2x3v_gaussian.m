% M is the highest tensor degree
function [] = ex_moment2x3v_gaussian(M)
%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','advection equation',... % name of example
    'initial_condition',@initial_condition,... % it is defined below
    'exact_solution',@exact_solution,...
    'ax',[0 1 0 1],... % extents of computational domain
    'n',[50 50],... % numbers of grid cells in each coordinate direction
    't_end',0.3,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',4,...
    'CFL',2,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'var_plot',1,...
    'to_plot',false...
    );


% we need the boundary matrix and the penalty matrix for all the
% boundaries
par.system.penalty_B = cell(par.num_bc,1);
par.system.penalty = cell(par.num_bc,1);
par.system.B = cell(par.num_bc,1);
par.system.rotator = cell(par.num_bc,1);
par.Kn = inf;

par.system.Ax = dvlp_Ax2D(M);
par.system.P = dvlp_Prod2D(M);
par.system.B{1} = dvlp_BInflow2D(M);

par.n_eqn = size(par.system.Ax,1);

par.system.B{1} = stabilize_boundary(par.system.Ax,par.system.B{1},M);

% rotation matrices for hermite polynomials
par.system.rotator = dvlp_RotatorCartesian(M,false);

par.normals_bc = [1,0;0,1;-1,0;0,-1];

par.system.Ay = par.system.rotator{2}' * par.system.Ax * par.system.rotator{2};


% id =1, x = 1
% id = 2 , y = 1
% id = 3, x = 0
% id = 4, y = 0
for i = 1 : par.num_bc
    par.system.B{i} = par.system.B{1}*par.system.rotator{i};
end

for i = 1 : par.num_bc
    An = par.system.Ax * par.normals_bc(i,1) + par.system.Ay * par.normals_bc(i,2);
    par.system.penalty{i} = dvlp_penalty_char(An,par.system.B{i});
end

for i = 1 : par.num_bc
    par.system.penalty_B{i} = par.system.penalty{i}*par.system.B{i};
end

result = solver(par);

filename = strcat('result_',num2str(M),'.txt');

dlmwrite(filename,result(1).X(:)','delimiter','\t');
dlmwrite(filename,result(1).Y(:)','delimiter','\t');

for i = 1 : length(result)
    dlmwrite(filename,result(i).sol(:)','delimiter','\t','-append');
end




% legend('2d code','1d code');
% filename = 'density.txt';
% dlmwrite(filename,result(1).X(:,1)','delimiter','\t');
% dlmwrite(filename,result(1).sol(:,1)','delimiter','\t','-append');
% 
% ========================================================================
% Run solver and study convergence
% ========================================================================
% resolution = [16 32 64 128 256];
% grid_spacing = [];
% 
% for k = 1:length(resolution)                       % Loop over various grid resolutions.
%     par.n = [1 1]*resolution(k);                   % Numbers of grid cells.
%     solution = solver(par);                         % Run solver.
%     error_temp = 0;
%     disp('Resolution :');
%     disp(par.n);
%     for j = 1:par.n_eqn                               % Loop over solution components.
%         X = solution(j).X;
%         Y = solution(j).Y;
%         PX = solution(j).PX;
%         PY = solution(j).PY;
%         U_theo = exact_solution(X,Y,par.t_end);     % Evaluate true solution.
%         error = abs(solution(j).sol-U_theo);               % Difference between num. and true sol.
%         int_x = dot(transpose(error),transpose(PX * error),2);  % integral along x.
%         int_xy = sum(PY*int_x);  % integral along xy.
%         error_temp = int_xy + error_temp;%Sc. L2 error.
%     end
%     error_L2(k) = sqrt(error_temp);
%     grid_spacing = [grid_spacing max(solution(1).h)];
% end
% 
% reference_line = exact_order(grid_spacing(1),error_L2(1),grid_spacing(end),2);
% 
% figure
% loglog(grid_spacing, error_L2, '-o',reference_line(1:2),reference_line(3:4),'-*');
% xlabel('h'), ylabel('l2-error');
% legend('numerical','second order');
% title('Convergence plot');
% 
% convg_order = log(error_L2(end)/error_L2(1))/log(grid_spacing(end)/grid_spacing(1));
% disp('convergence order');
% disp(convg_order);

end
%========================================================================
% Problem Specific Functions
%========================================================================

function f = initial_condition(x,y,j)
% Maxwellian/Gaussian
x0 = 0.5; % centered in the middle of domain
y0 = 0.5;
sigma_x = 0.1; % such that 6*sigma = 0.6 so in 0.2 length strip near 
sigma_y = 0.1;                           % boundary function value = 0

f = x * 0;

switch j
    case 1 
        f = exp( -((x-x0).^2 * 100) - ((y-y0).^2*100) );
        %f = exp( -((x-x0).^2 / (2*sigma_x.^2)) - ((y-y0).^2 / (2*sigma_y.^2)) );
end

end

function f = bc_inhomo(B,bc_id)    

    f = B(:,1)* 0;
end



