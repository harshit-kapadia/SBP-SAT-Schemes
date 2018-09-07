% we solve the Heat Conduction problem using discrete velocity method
function testing_DVM(nc)

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

par.t_plot = true;
par.n_eqn = 1;
par.n = [50 2];


par.system.Ax = sparse(1,1,1);
par.system.Ay = sparse(1,1,1);

result = solver_DVM_2x3v_upwind(par);

temp = cell(2,par.n_eqn);

for i = 1 : 2
    for j = 1 : par.n_eqn
        temp{i,j} = result(i,j).sol;
    end
end

% density = compute_density(temp,par.system.Ax,par.system.Ay,par.all_w);
% [ux,uy] = compute_velocity(temp,par.system.Ax,par.system.Ay,par.all_w);
% theta = compute_theta(temp,par.system.Ax,par.system.Ay,par.all_w);
% sigma_xx = compute_sigma_xx(temp,par.system.Ax,par.system.Ay,par.all_w);
% sigma_xy = compute_sigma_xy(temp,par.system.Ax,par.system.Ay,par.all_w);
% sigma_yy = compute_sigma_yy(temp,par.system.Ax,par.system.Ay,par.all_w);
% qx = compute_qx(temp,par.system.Ax,par.system.Ay,par.all_w);
% qy = compute_qy(temp,par.system.Ax,par.system.Ay,par.all_w);
% 
% filename = strcat('heated_cavity/result_DVM_',num2str(nc),'.txt');
% dlmwrite(filename,result(1,1).X(:)','delimiter','\t','precision',10);
% dlmwrite(filename,result(1,1).Y(:)','delimiter','\t','precision',10,'-append');
% dlmwrite(filename,density(:)','delimiter','\t','-append','precision',10);
% 
% dlmwrite(filename,ux(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,uy(:)','delimiter','\t','-append','precision',10);
% 
% dlmwrite(filename,theta(:)','delimiter','\t','-append','precision',10);
% 
% dlmwrite(filename,sigma_xx(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,sigma_xy(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,sigma_yy(:)','delimiter','\t','-append','precision',10);
% 
% dlmwrite(filename,qx(:)','delimiter','\t','-append','precision',10);
% dlmwrite(filename,qy(:)','delimiter','\t','-append','precision',10);

end

function thetaW = compute_thetaW(bc_id,t)
thetaW = 0;

% if bc_id == 4 % bottom boundary
%     if t <= 1
%         thetaW = exp(-1/(1-(t-1)^2)) * exp(1);
%     else
%         thetaW = 1;
%     end
% end

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
function f = ic(x,y)
f = cell(2,1);

f{1} = sin(pi * x);
f{2} = sin(pi * x);
 
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




