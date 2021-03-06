% M is the highest tensor degree
function [] = ex_unsteady_lid_driven_cavity(M)
%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','advection equation',... % name of example
    'initial_condition',@initial_condition,... % it is defined below
    'exact_solution',@exact_solution,...
    'ax',[0 1 0 1],... % extents of computational domain
    'n',[100 100],... % numbers of grid cells in each coordinate direction
    't_end',1.0,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',4,...
    'CFL',1,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'var_plot',1,... % useless
    'to_plot',false,...
    'compute_density',@compute_density,...
    'compute_ux',@compute_ux,...
    'compute_uy',@compute_uy,...
    'compute_theta',@compute_theta,...
    'compute_sigma_xx',@compute_sigma_xx,...
    'compute_sigma_xy',@compute_sigma_xy,...
    'compute_sigma_yy',@compute_sigma_yy,...
    'compute_qx',@compute_qx,...
    'compute_qy',@compute_qy,...
    'write_solution',@write_solution,...
    'steady_state',false...
    );

% file where the output is written
par.output_filename = 'unsteady_lid_driven_cavity';

par.M = M;

% % incase M if greater then 3 then read the written data. (Only read M + 2)
% par.M = M;
% 
% if M > 3
%     %par.previous_M_data = dlmread(filename,'\t');
%     par.previous_M_data = zeros(1,1);
% else
%     par.previous_M_data = zeros(1,1);
% end

% we need the boundary matrix and the penalty matrix for all the
% boundaries
par.system.penalty_B = cell(par.num_bc,1);
par.system.penalty = cell(par.num_bc,1);
par.system.B = cell(par.num_bc,1);
par.system.rotator = cell(par.num_bc,1);
par.Kn = 0.1;

par.system.Ax = dvlp_Ax2D(M);
par.system.P = dvlp_Prod2D(M);
par.system.BIn = dvlp_BInflow2D(M);
par.system.BWall = dvlp_BWall2D(M);

par.n_eqn = size(par.system.Ax,1);

par.system.BIn = stabilize_boundary(par.system.Ax,par.system.BIn,M);
par.system.BWall = stabilize_boundary(par.system.Ax,par.system.BWall,M);

% rotation matrices for hermite polynomials
par.system.rotator = dvlp_RotatorCartesian(M,false);

par.normals_bc = [1,0;0,1;-1,0;0,-1];

par.system.Ay = par.system.rotator{2}' * par.system.Ax * par.system.rotator{2};

% id =1, x = 1
% id = 2 , y = 1
% id = 3, x = 0
% id = 4, y = 0
par.system.B{1} = par.system.BWall;
par.system.B{2} = par.system.BWall * par.system.rotator{2};
par.system.B{3} = par.system.BWall * par.system.rotator{3};
par.system.B{4} = par.system.BWall * par.system.rotator{4};

% first boundary
par.system.penalty{1} = dvlp_penalty_char(par.system.Ax,par.system.B{1});

% rotator for different boundaries
rotator = dvlp_RotatorCartesian(M,false);

for i = 1 : par.num_bc
    par.system.penalty{i} = rotator{i}' * par.system.penalty{1};
end

for i = 1 : par.num_bc
    par.system.penalty_B{i} = par.system.penalty{i}*par.system.B{i};
end

par.previous_M_data = 0; % actually useless for this example (unsteady)
result = solver(par);

end

function [penalty] = dvlp_penalty_odd(Ax,M)

[odd_ID,~] = get_id_Odd(M);

odd_ID = flatten_cell(odd_ID);

% we pick the columns which get multiplied by the odd variables
penalty = Ax(:,odd_ID);
end

% read data contains the already read files
function f = initial_condition(x,y,j,read_data)

f = x * 0;

% moments_read_data = size(read_data,2) - 2;
% 
% % if we have the data from the previous moment system then we initialize
% % from there
% 
% if j <= moments_read_data
%         f = reshape(read_data(j+2,:),size(x));
% end

end

function f = bc_inhomo(B,bc_id,t)    

    f = B(:,1)* 0;
    
    % tangential velocity
    % tangential direction is -x. Hence the minus (? for unsteady)
    vt = sin(4*pi*t) ;

    % the top wall moves
    switch bc_id
        case 2
            % The B matrix has already been multiplied by the rotator. As a
            % result we consider the second column. But in the DG
            % framework, B is not multiplied by the rotator and so we take
            % the third column.
            f = vt * B(:,2); 
    end
end

function [] = write_solution(temp,par,X,Y,M,residual)
density = par.compute_density(temp);
ux = par.compute_ux(temp);
uy = par.compute_uy(temp);
theta = par.compute_theta(temp);
sigma_xx = par.compute_sigma_xx(temp);
sigma_xy = par.compute_sigma_xy(temp);
sigma_yy = par.compute_sigma_yy(temp);
qx = par.compute_qx(temp);
qy = par.compute_qy(temp);

filename = strcat(par.output_filename,'/result_M',num2str(M),'.txt');
dlmwrite(filename,X(:)','delimiter','\t','precision',10);
dlmwrite(filename,Y(:)','delimiter','\t','precision',10,'-append');
dlmwrite(filename,density(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,ux(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,uy(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,theta(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,sigma_xx(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,sigma_xy(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,sigma_yy(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,qx(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,qy(:)','delimiter','\t','-append','precision',10);

filename_residual = strcat(par.output_filename,'/residual_M',num2str(M),'.txt');
dlmwrite(filename_residual,residual(:)','delimiter','\t','-append','precision',10);
end

function f = compute_density(data)
f = data{1};
end

function f = compute_ux(data)
f = data{2};
end

function f = compute_uy(data)
f = data{3};
end

function f = compute_theta(data)
f = sqrt(2) * (data{4} + data{6} + data{7})/3;
end

function f = compute_sigma_xx(data)
f = sqrt(2) * (2 * data{4} - data{6} - data{7})/3;
end

function f = compute_sigma_xy(data)
% 110
f = data{5};
end

function f = compute_sigma_yy(data)
f = sqrt(2) * (-data{4} + 2 * data{6} - data{7})/3;
end

function f = compute_qx(data)
shift =7;   % shift for tensor degrees which come before 3 
% (3,0,0) + (1,2,0) + (1,0,2)
f = (sqrt(3) * data{shift + 1} + data{shift + 3} + data{shift + 4})/sqrt(2);
end

function f = compute_qy(data)
shift =7;   % shift for tensor degrees which come before 3 
% (0,3,0) + (2,1,0) + (0,1,2)
f = (sqrt(3) * data{shift + 5} + data{shift + 2} + data{shift + 6})/sqrt(2);
end

