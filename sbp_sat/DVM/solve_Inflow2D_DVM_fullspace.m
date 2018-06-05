% we solve the Heat Conduction problem using discrete velocity method 
function solve_Inflow2D_DVM_fullspace(nc)

par = struct(...
'name','Inflow DVM',... % name of example
'ic',@ic,... % initial conditions
'bc_inhomo',@bc_inhomo,... % source term (defined below)
'ax',[0 1 0 1],... % coordinates of computational domain
 't_end',0.3,... % the end time of the computation
 'CFL',2.0,...      % the crude cfl number
 'num_bc',4,... % number of boundaries in the domain
 'RK_order',4,... 
'output',@output,... % problem-specific output routine (defined below)
'save_during',false, ... % should we save during the computation
'compute_during', @compute_during, ...
'compute_density',@compute_density, ...
'compute_velocity',@compute_velocity, ...
'compute_theta',@compute_theta...
);

par.Kn = 0.1;
par.t_plot = false;
par.n_eqn = nc^2;
% number of points in the spatial discretization
par.n = [50 50];

[par.x,par.w] = gauss_quadrature(nc,-5,5);

% velocity grid in one dimension
[v_grid1D,perm] = sort(par.x');
w_grid1D = par.w';
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
result= solver_DVM_2x3v(par);

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

filename = 'result_Inflow_Comp/result_DVM_fullspace.txt';
dlmwrite(filename,result(1,1).X(:)','delimiter','\t','precision',10);
dlmwrite(filename,result(1,1).Y(:)','delimiter','\t','precision',10,'-append');
dlmwrite(filename,density(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,ux(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,uy(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,theta(:)','delimiter','\t','-append','precision',10);
dlmwrite(filename,sigma_xx(:)','delimiter','\t','-append','precision',10);

end


% working on inflow boundaries, we consider vacuum boundary conditions
function f = bc_inhomo(B,bc_id,Ax,Ay,id_sys,U,all_weights,t)

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
