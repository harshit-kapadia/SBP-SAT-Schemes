% we solve the Heat Conduction problem using discrete velocity method
function ex_heated_cavity(nc)

par = struct(...
    'name','Inflow DVM',... % name of example
    'ic',@ic,... % initial conditions
    'bc_inhomo_inflow',@bc_inhomo_inflow,... % boundary inhomogeneity for inflow
    'bc_inhomo_wall',@bc_inhomo_wall,... % boundary inhomogeneity for wall 
    'ax',[0 1 0 1],... % coordinates of computational domain
    't_end',2.0,... % the end time of the computation
    'CFL',2.0,...      % the crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'RK_order',4,...
    'output',@output,... % problem-specific output routine (defined below)
    'save_during',false, ... % should we save during the computation
    'compute_during', @compute_during, ...
    'compute_density',@compute_density, ...
    'compute_velocity',@compute_velocity, ...
    'compute_theta',@compute_theta, ...
    'compute_rhoW_prep',@compute_rhoW_prep,...
    'compute_thetaW',@compute_thetaW,...
    'steady_state',true...
    );

par.wall_boundary = [true,true,true,true]; % which of the boundaries are wall boundaries
par.Kn = 0.1;
par.t_plot = false;
par.n_eqn = (2 * nc) * (2 * nc);
% number of points in the spatial discretization
par.n = [100 100];

[par.x_m,par.w_m] = gauss_quadrature(nc,-5,0);
[par.x_p,par.w_p] = gauss_quadrature(nc,0,5);

% velocity grid in one dimension
[v_grid1D,perm] = sort([par.x_m' par.x_p']);
w_grid1D = [par.w_m',par.w_p'];
w_grid1D = w_grid1D(perm);

% velocity grid in two dimensions
[vx_grid2D,vy_grid2D] = meshgrid(v_grid1D,v_grid1D);

% the quadrature weights grid in two dimensions
[wx_grid2D,wy_grid2D] = meshgrid(w_grid1D,w_grid1D);

w = wx_grid2D.*wy_grid2D;

% Convert the matrix of grid points into an array and store as a hyperbolic
% system. Do the same for both the directions.
par.system.Ax = diag(vx_grid2D(:));
par.system.Ay = diag(vy_grid2D(:));
par.all_w = w(:);


% at x = 1 we prescribe boundary to negative velocities and at x = 0 we
% precribe boundary conditions for the positive velocities.
par.system.B{1} = diag(double(diag(par.system.Ax)<0));
par.system.B{3} = diag(double(diag(-par.system.Ax)<0));

% same as the boundaries along x but for the y direction.
par.system.B{2} = diag(double(diag(par.system.Ay)<0));
par.system.B{4} = diag(double(diag(-par.system.Ay)<0));

% prescribe a value to the negative velocities
par.system.penalty{1} = (par.system.Ax-abs(par.system.Ax))/2;
par.system.penalty{3} = (-par.system.Ax-abs(par.system.Ax))/2;

par.system.penalty{2} = (par.system.Ay-abs(par.system.Ay))/2;
par.system.penalty{4} = (-par.system.Ay-abs(par.system.Ay))/2;

% develop the penalty * B matrix
for i = 1 : par.num_bc
    par.system.penalty_B{i} = par.system.penalty{i} * par.system.B{i};
end

% inner product matrix of hermite polynomials computed with gauss
% quadrature integration
par.mass_matrix = mass_matrix_quad(par.system.Ax,par.system.Ay,par.all_w);
par.inv_mass_matrix = inv(par.mass_matrix);
% store the value of f0 at all the quadrature points
par.value_f0 = arrayfun(@(x,y) f0(x,y),diag(par.system.Ax),diag(par.system.Ay));

[par.rhoW_vect, par.rhoW_value, par.pos_U]...
    = compute_rhoW_prep_simple(par.system.Ax, par.system.Ay, par.all_w) ;

result = solver_DVM_2x3v(par);

temp = cell(2,par.n_eqn);

for i = 1 : 2
    for j = 1 : par.n_eqn
        temp{i,j} = result(i,j).sol;
    end
end

density = compute_density(temp,par.system.Ax,par.system.Ay,par.all_w);
[ux,uy] = compute_velocity(temp,par.system.Ax,par.system.Ay,par.all_w);
theta = compute_theta(temp,par.system.Ax,par.system.Ay,par.all_w);
sigma_xx = compute_sigma_xx(temp,par.system.Ax,par.system.Ay,par.all_w);
sigma_xy = compute_sigma_xy(temp,par.system.Ax,par.system.Ay,par.all_w);
sigma_yy = compute_sigma_yy(temp,par.system.Ax,par.system.Ay,par.all_w);
qx = compute_qx(temp,par.system.Ax,par.system.Ay,par.all_w);
qy = compute_qy(temp,par.system.Ax,par.system.Ay,par.all_w);

filename = strcat('heated_cavity/result_DVM_',num2str(nc),'.txt');
dlmwrite(filename,result(1,1).X(:)','delimiter','\t','precision',10);
dlmwrite(filename,result(1,1).Y(:)','delimiter','\t','precision',10,'-append');
dlmwrite(filename,density(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,ux(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,uy(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,theta(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,sigma_xx(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,sigma_xy(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,sigma_yy(:)','delimiter','\t','-append','precision',10);

dlmwrite(filename,qx(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,qy(:)','delimiter','\t','-append','precision',10);

end

function thetaW = compute_thetaW(bc_id,t)
thetaW = 0;

if bc_id == 4 % bottom boundary
    if t <= 1
        thetaW = exp(-1/(1-(t-1)^2)) * exp(1);
    else
        thetaW = 1;
    end
end

end

% wall boundary condition
function f = bc_inhomo_wall(B,bc_id,Ax,Ay,rhoW,id_sys,t)

% tangential velocity and normal velocity of the wall
ux = 0;
uy = 0;

thetaW = compute_thetaW(bc_id,t);

% we initialiase with a cell because of the data structure used in the code
f = cell(1,length(B));
% we find the location of velocities which have a negative normal velocity i.e. the
% velocities to which we need to prescribe a boundary condition
id = find(diag(B) == 1);

for i = 1 : length(id)
    % assumption is that rhoW changes along the grid boundary while thetaW
    % does not
    f{id(i)} = compute_fM(Ax,Ay,rhoW,ux,uy,thetaW,id(i),id_sys);
end

end

function f = bc_inhomo_inflow(B,bc_id,Ax,Ay,id_sys,U,all_weights,t)

    rho = 0;
    ux = 0;
    uy = 0;
    
    if t <= 1
        thetaIn = exp(-1/(1-(t-1)^2)) * exp(1);
    else
        thetaIn = 1;
    end
  
    id = find(diag(B) == 1);
    
    f = diag(B) * 0; 
end

% Compute the Maxwellian, there are two maxellians, one corresponding to
% g and one for h.
function f = compute_fM(Ax,Ay,rho,ux,uy,theta,id,id_sys)

vx = Ax(id,id);
vy = Ay(id,id);

switch id_sys
    
    % maxwellian for g
    case 1
        f = (vx * ux + vy * uy + ((vx^2 + vy^2)/2 - 1) * theta + rho) * f0(vx,vy);
        % maxwellian for h
    case 2
        f = theta * f0(vx,vy)/sqrt(2);
end

end

% we have two systems of the same type
function f = ic(x,y,Ax,Ay,mass_matrix,inv_mass_matrix,value_f0)

%rho = exp(-(y-0.5).*(y-0.5)*100);
%rho = exp(-(x-0.5).*(x-0.5)*100-(y-0.5).*(y-0.5)*100);
rho = x * 0;
ux = x * 0;
uy = x * 0;
theta = x * 0;

f = store_minimized_entropy(Ax,Ay, ...
     mass_matrix,inv_mass_matrix,rho,ux,uy,theta,value_f0);
 
end

% f0 corresponding to the 2d velocity space
function f = f0(vx,vy)
f = exp(-(vx^2+vy^2)/2)/(2 * pi);
end

% compute the density on the complete grid
function f = compute_density(U,Ax,Ay,all_weights)

[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

for k = 1 : length(all_weights)
    temp = U{1,k};
    
    f = f + all_weights(k)*temp;
    
end

end

% compute velocity in the x and y direction
function [vx,vy] = compute_velocity(U,Ax,Ay,all_weights)

[points_x,points_y] = size(U{1,1});

vx = zeros(points_x,points_y);
vy = zeros(points_x,points_y);

% scale by velocity
all_weights_x = diag(Ax).*all_weights;
all_weights_y = diag(Ay).*all_weights;

for k = 1 : length(all_weights_x)
    temp = U{1,k};
    
    vx = vx + (all_weights_x(k)*temp);
    vy = vy + (all_weights_y(k)*temp);
    
end

end

% compute the temperature
function f = compute_theta(U,Ax,Ay,all_weights)

[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

% scale by the He2(xi_1)
all_weights_x = (diag(Ax).*diag(Ax)-1).*all_weights/sqrt(2);

% scale by He2(xi_2)
all_weights_y = (diag(Ay).*diag(Ay)-1).*all_weights/sqrt(2);


for k = 1 : length(all_weights)
    temp_g = U{1,k};
    temp_h = U{2,k};
    
    f = f + sqrt(2) * ((all_weights_x(k) + all_weights_y(k))*temp_g + ...
        all_weights(k)*temp_h)/3;
    
end

end

% compute the sigma xx
% compute the temperature
function f = compute_sigma_xx(U,Ax,Ay,all_weights)

[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

% scale by the He2(xi_1)
all_weights_x = (diag(Ax).*diag(Ax)-1).*all_weights/sqrt(2);

% scale by He2(xi_2)
all_weights_y = (diag(Ay).*diag(Ay)-1).*all_weights/sqrt(2);


for k = 1 : length(all_weights)
    f = f + sqrt(2) * ((2 * all_weights_x(k) - all_weights_y(k))*U{1,k} - ...
        all_weights(k)*U{2,k})/3;
end

end

function f = compute_sigma_xy(U,Ax,Ay,all_weights)

[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

% scale by the He1(xi_1)He1(xi_2)
all_weights_xy = (diag(Ax).*diag(Ay)).*all_weights;

for k = 1 : length(all_weights)
    f = f +  all_weights_xy(k)*U{1,k};
end

end

function f = compute_sigma_yy(U,Ax,Ay,all_weights)

[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

% scale by the He2(xi_1)
all_weights_x = (diag(Ax).*diag(Ax)-1).*all_weights/sqrt(2);

% scale by He2(xi_2)
all_weights_y = (diag(Ay).*diag(Ay)-1).*all_weights/sqrt(2);


for k = 1 : length(all_weights)
    f = f + sqrt(2) * ((-all_weights_x(k) + 2 * all_weights_y(k))*U{1,k} - ...
        all_weights(k)*U{2,k})/3;
end

end

function f = compute_qx(U,Ax,Ay,all_weights)
[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

vx = diag(Ax);
vy = diag(Ay);

% scale the weights by He_3(xi_1)
all_weights_He_3_x = He_3(vx).*all_weights;

all_weights_He_1_x_He_2_y = (He_1(vx).*He_2(vy)).*all_weights;

all_weights_He_1_x = He_1(vx).*all_weights;

for i = 1 : length(all_weights)
    temp_g = U{1,i};
    temp_h = U{2,i};
    
    f = f + (sqrt(3) * all_weights_He_3_x(i) * temp_g + all_weights_He_1_x_He_2_y(i) * temp_g + ...
            all_weights_He_1_x(i) * temp_h)/sqrt(2);
end

end

function f = compute_qy(U,Ax,Ay,all_weights)
[points_x,points_y] = size(U{1,1});

f = zeros(points_x,points_y);

vx = diag(Ax);
vy = diag(Ay);

% scale the weights by He_3(xi_1)
all_weights_He_3_y = He_3(vy).*all_weights;

all_weights_He_2_x_He_1_y = (He_2(vx).*He_1(vy)).*all_weights;

all_weights_He_1_y = He_1(vy).*all_weights;

for i = 1 : length(all_weights)
    temp_g = U{1,i};
    temp_h = U{2,i};
    
    f = f + (sqrt(3) * all_weights_He_3_y(i) * temp_g + all_weights_He_2_x_He_1_y(i) * temp_g + ...
            all_weights_He_1_y(i) * temp_h)/sqrt(2);
end

end

% hermite polynomial of degree 3
function f = He_3(c)
f = (c.*c-3).*c/sqrt(6);
end

function f = He_2(c)
f = (c.*c - 1)/sqrt(2);
end

function f = He_1(c)
f = c;
end




