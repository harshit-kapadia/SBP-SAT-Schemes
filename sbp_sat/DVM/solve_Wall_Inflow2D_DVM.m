% we solve the Heat Conduction problem using discrete velocity method
function solve_Wall_Inflow2D_DVM(nc)

par = struct(...
    'name','Wall+Inflow DVM',... % name of example
    'ic',@ic,... % initial conditions
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'bc_inhomo_2',@bc_inhomo_2,... % source term (defined below)
    'ax',[0 1 0 1],... % coordinates of computational domain
    't_end',5.0,... % the end time of the computation
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
    'compute_thetaW',@compute_thetaW...
    );

par.bc_indicator = [1 0 1 0] ; % 1 for Wall b.c.

par.Kn = 0.1;
par.t_plot = true;
par.n_eqn = (2 * nc) * (2 * nc);
% number of points in the spatial discretization
par.n = [50 50];

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

par.system.Ax = diag(vx_grid2D(:));
par.system.Ay = diag(vy_grid2D(:));
par.all_w = w(:);


% at x = 1 we prescribe boundary to negative velocities
par.system.B{1} = diag(double(diag(par.system.Ax)<0));
par.system.B{3} = diag(double(diag(-par.system.Ax)<0));

par.system.B{2} = diag(double(diag(par.system.Ay)<0));
par.system.B{4} = diag(double(diag(-par.system.Ay)<0));

% prescribe a value to the negative velocities
par.system.penalty{1} = (par.system.Ax-abs(par.system.Ax))/2;
par.system.penalty{3} = (-par.system.Ax-abs(par.system.Ax))/2;

par.system.penalty{2} = (par.system.Ay-abs(par.system.Ay))/2;
par.system.penalty{4} = (-par.system.Ay-abs(par.system.Ay))/2;

% develop the penalty*B matrix
for i = 1 : par.num_bc
    par.system.penalty_B{i} = par.system.penalty{i} * par.system.B{i};
end

% inner product matrix of hermite polynomials computed with gauss
% quadrature integration
par.mass_matrix = mass_matrix_quad(par.system.Ax,par.system.Ay,par.all_w);
par.inv_mass_matrix = inv(par.mass_matrix);
par.value_f0 = arrayfun(@(x,y) f0(x,y),diag(par.system.Ax),diag(par.system.Ay));

[par.rhoW_vect, par.rhoW_value, par.pos_U]...
    = compute_rhoW_prep(par.system.Ax, par.system.Ay, par.all_w) ;

result = solver_DVM_2x3v_g_varying(par);

%two different systems and so two rows
temp = cell(2,par.n_eqn);

for i = 1 : 2
    for j = 1 : par.n_eqn
        temp{i,j} = result(i,j).sol;
    end
end

density = compute_density(temp,par.system.Ax,par.system.Ay,par.all_w);
[ux,uy] = compute_velocity(temp,par.system.Ax,par.system.Ay,par.all_w);
theta = compute_theta(temp,par.system.Ax,par.system.Ay,par.all_w);
sigma_xx = compute_sigmaxx(temp,par.system.Ax,par.system.Ay,par.all_w);

% filename = 'result_Wall2D_Comp/result_DVM.txt';
% dlmwrite(filename,result(1,1).X(:)','delimiter','\t','precision',10);
% dlmwrite(filename,result(1,1).Y(:)','delimiter','\t','precision',10,'-append');
% dlmwrite(filename,density(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,ux(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,uy(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,theta(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,sigma_xx(:)','delimiter','\t','-append','precision',10);

end


% Wall boundary condition
function f = bc_inhomo(B,bc_id,Ax,Ay,thetaW,rhoW,id_sys)

% tangential velocity and normal velocity of the wall
ux = 0;
uy = 0;

f = cell(1,length(B));

id = find(diag(B) == 1);

for i = 1 : length(id)
    f{id(i)} = compute_fM(Ax,Ay,rhoW,ux,uy,thetaW,id(i),id_sys);
end

end

% Inflow boundary condition
function f = bc_inhomo_2(B,bc_id,Ax,Ay,id_sys,U,all_weights,t)

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

% we have two systems of the same type
function f = ic(x,y,Ax,Ay,mass_matrix,inv_mass_matrix,value_f0)

rho = exp(-(x-0.5).*(x-0.5)*100-(y-0.5).*(y-0.5)*100);
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
function f = compute_sigmaxx(U,Ax,Ay,all_weights)

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

% compute the Maxwellian, there are two maxellians, one corresponding to
% g and one for h
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

% Temperature at the wall, needed for Wall boundary implementation
function thetaW = compute_thetaW(t)

thetaW{1} = 0 ;
thetaW{2} = 0 ;
thetaW{3} = 0 ;
thetaW{4} = 0 ;

% if t_temp(RK) <= 1
%     temp_thetaW = exp(-1/(1-(t-1)^2)) * exp(1) ;
% else
%     temp_thetaW = 1 ;
% end
% thetaW{1} = temp_thetaW ;
% thetaW{2} = temp_thetaW ;
% thetaW{3} = -temp_thetaW ;
% thetaW{4} = -temp_thetaW ;

end

function [vect,value,pos_U] = compute_rhoW_prep(Ax,Ay,all_weights)

all_vx = diag(Ax);
all_vy = diag(Ay);

pos_vx_p = find(all_vx > 0); % positive x-vel
pos_vx_m = find(all_vx < 0); % negative x-vel

pos_vy_p = find(all_vy > 0);
pos_vy_m = find(all_vy < 0);

vx_p_x = all_vx(pos_vx_p);
vx_m_x = all_vx(pos_vx_m);
vy_p_x = all_vy(pos_vx_p); % y-vel corresponding to positive x-vel
vy_m_x = all_vy(pos_vx_m); % y-vel corresponding to negative x-vel

vx_p_y = all_vx(pos_vy_p);
vx_m_y = all_vx(pos_vy_m);
vy_p_y = all_vy(pos_vy_p);
vy_m_y = all_vy(pos_vy_m);

num_pos_x = length(vx_p_x);
num_neg_x = length(vx_m_x);

num_pos_y = length(vy_p_y);
num_neg_y = length(vy_m_y);

weight_vx_m = all_weights(pos_vx_m);
weight_vx_p = all_weights(pos_vx_p);

weight_vy_m = all_weights(pos_vy_m);
weight_vy_p = all_weights(pos_vy_p);


vect = cell(1,4);
value = cell(1,4);

int_f0 = 0; % actually will store int_f0_vn (vn is normal velocity at boundary)
int_f0_vx_sq = 0; % integral(f0 * vx_sq * vn)
int_f0_vy_sq = 0; % integral(f0 * vy_sq * vn)
% boundary at x = 1
% integrate the maxwellian with vx_m_x (incoming normal velocity at boundary)
for i = 1 : num_neg_x
    int_f0 = int_f0 + vx_m_x(i) * f0(vx_m_x(i),vy_m_x(i)) * weight_vx_m(i);
    int_f0_vx_sq = int_f0_vx_sq + vx_m_x(i) * vx_m_x(i)^2 * ...
        f0(vx_m_x(i),vy_m_x(i)) * weight_vx_m(i);
    int_f0_vy_sq = int_f0_vy_sq + vx_m_x(i) * vy_m_x(i)^2 * ...
        f0(vx_m_x(i),vy_m_x(i)) * weight_vx_m(i);
end
value{1} = ( (int_f0_vx_sq + int_f0_vy_sq)/2-int_f0 ) / int_f0;
% prepration to
% integrate the kinetic solution with vx_p_x (outgoing normal velocity at boundary)
% int_f_vn = sum((vx_p_x.*U(pos_vx_p)).*weight_vx_p); (actual integration)
vect{1} = ( vx_p_x.*weight_vx_p ) / int_f0;


int_f0 = 0;
int_f0_vx_sq = 0;
int_f0_vy_sq = 0;
% integrate the maxwellian with vy_m_y (incoming normal velocity at boundary)
for i = 1 : num_neg_y
    int_f0 = int_f0 + vy_m_y(i) * f0(vx_m_y(i),vy_m_y(i)) * weight_vy_m(i);
    int_f0_vx_sq = int_f0_vx_sq + vy_m_y(i) * vx_m_y(i)^2 * ...
        f0(vx_m_y(i),vy_m_y(i)) * weight_vy_m(i);
    int_f0_vy_sq = int_f0_vy_sq + vy_m_y(i) * vy_m_y(i)^2 * ...
        f0(vx_m_y(i),vy_m_y(i)) * weight_vy_m(i);
end
value{2} = ( (int_f0_vx_sq + int_f0_vy_sq)/2-int_f0 ) / int_f0;
% prepration to
% integrate the kinetic solution with vy_p_y (outgoing normal velocity at boundary)
vect{2} = ( vy_p_y.*weight_vy_p ) / int_f0;


int_f0 = 0;
int_f0_vx_sq = 0;
int_f0_vy_sq = 0;
% integrate the maxwellian with vx_p_x (incoming normal velocity at boundary)
for i = 1 : num_pos_x
    int_f0 = int_f0 + vx_p_x(i) * f0(vx_p_x(i),vy_p_x(i)) * weight_vx_p(i);
    int_f0_vx_sq = int_f0_vx_sq + vx_p_x(i) * vx_p_x(i)^2 * ...
        f0(vx_p_x(i),vy_p_x(i)) * weight_vx_p(i);
    int_f0_vy_sq = int_f0_vy_sq + vx_p_x(i) * vy_p_x(i)^2 * ...
        f0(vx_p_x(i),vy_p_x(i)) * weight_vx_p(i);
end
value{3} = ( (int_f0_vx_sq + int_f0_vy_sq)/2-int_f0 ) / int_f0;
% prepration to
% integrate the kinetic solution with vx_m_x (outgoing normal velocity at boundary)
vect{3} = ( vx_m_x.*weight_vx_m ) / int_f0;


int_f0 = 0;
int_f0_vx_sq = 0;
int_f0_vy_sq = 0;
% integrate the maxwellian with vy_p_y (incoming normal velocity at boundary)
for i = 1 : num_pos_y
    int_f0 = int_f0 + vy_p_y(i) * f0(vx_p_y(i),vy_p_y(i)) * weight_vy_p(i);
    int_f0_vx_sq = int_f0_vx_sq + vy_p_y(i) * vx_p_y(i)^2 * ...
        f0(vx_p_y(i),vy_p_y(i)) * weight_vy_p(i);
    int_f0_vy_sq = int_f0_vy_sq + vy_p_y(i) * vy_p_y(i)^2 * ...
        f0(vx_p_y(i),vy_p_y(i)) * weight_vy_p(i);
end
value{4} = ( (int_f0_vx_sq + int_f0_vy_sq)/2-int_f0 ) / int_f0;
% prepration to
% integrate the kinetic solution with vy_m_y (outgoing normal velocity at boundary)
vect{4} = ( vy_m_y.*weight_vy_m ) / int_f0;

% % final target
% rhoW = -int_f_vn-thetaW * ((int_f0_vx_sq + int_f0_vy_sq)/2-int_f0);
% rhoW = rhoW/int_f0; % division is performed while computing value{} and vect{}

pos_U{1} = pos_vx_p ;
pos_U{2} = pos_vy_p ;
pos_U{3} = pos_vx_m ;
pos_U{4} = pos_vy_m ;
end
