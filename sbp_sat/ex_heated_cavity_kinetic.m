% M is the highest tensor degree
function [] = ex_heated_cavity_kinetic(M)
%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','advection equation',... % name of example
    'initial_condition',@initial_condition,... % it is defined below
    'exact_solution',@exact_solution,...
    'ax',[0 1 0 1],... % extents of computational domain
    'n',[50 50],... % numbers of grid cells in each coordinate direction
    't_end',0.5,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',4,...
    'CFL',1,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'var_plot',1,...
    'to_plot',true,...
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
    'steady_state',true...
    );

% file where the output is written
par.output_filename = strcat('heated_cavity_kinetic/result_M',num2str(M),'.txt');

% incase M if greater then 3 then read the written data. (Only read M + 2)
par.M = M;

if M > 3
    %par.previous_M_data = dlmread(filename,'\t');
    par.previous_M_data = zeros(1,1);
else
    par.previous_M_data = zeros(1,1);
end

% we need the boundary matrix and the penalty matrix for all the
% boundaries
par.system.penalty_B = cell(par.num_bc,1);
par.system.penalty = cell(par.num_bc,1);
par.system.B = cell(par.num_bc,1);
par.system.rotator = cell(par.num_bc,1);
par.Kn = 0.1;

par.system.Ax = dvlp_Ax2D(M);
par.system.P = dvlp_Prod2D(M);

par.n_eqn = size(par.system.Ax,1);

% rotation matrices for hermite polynomials
par.system.rotator = dvlp_RotatorCartesian(M,false);

% flux along the y direction
par.system.Ay = par.system.rotator{2}' * par.system.Ax * par.system.rotator{2};

% first boundary
par.system.penalty{1} = dvlp_penalty_kinetic(par.system.Ax,M);


for i = 1 : par.num_bc
    par.system.penalty{i} = par.system.rotator{i}' * par.system.penalty{1} * par.system.rotator{i};
    par.system.B{i} = eye(par.n_eqn); % prescribing B here is just a formality and 
                                      % is being done only to preserve the code structure
end

for i = 1 : par.num_bc
    par.system.penalty_B{i} = par.system.penalty{i}*par.system.B{i};
end

par.previous_M_data = 0;
result = solver(par);

end

% read data contains the already read files
function f = initial_condition(x,y,j,read_data)

f = x * 0;

moments_read_data = size(read_data,2) - 2;

% if we have the data from the previous moment system then we initialize
% from there

if j <= moments_read_data
        f = reshape(read_data(j+2,:),size(x));
end

end

function f = bc_inhomo(B,bc_id,t)    

    f = B(:,1)* 0;
    
    if t <= 1
        thetaIn = exp(-1/(1-(t-1)^2)) * exp(1);
    else
        thetaIn = 1;
    end

    switch bc_id
        case 4
            f(4) = thetaIn/sqrt(2); 
            f(6) = thetaIn/sqrt(2);
            f(7) = thetaIn/sqrt(2);
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

filename = par.output_filename;
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

filename = strcat('heated_cavity_odd/residual_M',num2str(M),'.txt');
dlmwrite(filename,residual(:)','delimiter','\t','-append','precision',10);
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

