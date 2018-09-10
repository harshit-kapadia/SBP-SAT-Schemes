function output = solver_DVM_2x3v(par)

if ~isfield(par,'save_during'), par.save_during = false; end % default value of save during the computation

if ~isfield(par,'ic'),       par.ic = @zero; end% Default: no init. cond.

if ~isfield(par,'source'),  par.source = @zero; end % Default: no source.


if par.num_bc ~=4
    assert(1 == 0, 'not valid num bc'); 
end   

% if we are not looking for the steady state then take the residual to be
% simply zero. 
if ~isfield(par,'steady_state')
    par.steady_state = false;
end

% find which of the given data is time dependent
time_dep = [nargin(par.source)>2 nargin(par.bc_inhomo_inflow)>7];
        
% if the number of arguments are greater than 3 then definitely we have 
% an anisotropic source term.
if ~isfield(par,'source_ind') 
    par.source_ind = 1:par.n_eqn;  
end

%% Grid Setup

% size of the grid. h1 for the x direction and h2 for the y direction
h = (par.ax([2 4])-par.ax([1 3]))./par.n; % ./ does right-array division

% Set of grid nodes in x and y direction
x{1} = par.ax(1):h(1):par.ax(2);
y{1} = par.ax(3):h(2):par.ax(4);

% Generate both grids. X contains the x coordinates and Y contains
% the y coordinates.


[X,Y] = cellfun(@ndgrid,x,y,'Un',0);


%% Time-step

% a crude approximation for delta_t
par.dt = min(h)/max(abs(diag(par.system.Ax)))/par.CFL; % 1 eigenv.

disp('delta t');
disp(par.dt);

% with largest magnitude

%% Setup difference operators
% differential opeartor for the X grid, i.e. dx matrix
[DX{1},PX{1}] = sbp_collocated_2(par.n(1),h(1));

% differential operators for the Y grid, i.e. dy matrix
[DY{1},PY{1}] = sbp_collocated_2(par.n(2),h(2));

% we transpose DY for the data structure
DY = cellfun(@transpose,DY,'Un',0);

%% initialisation of solution variables
% we have two systems now, one corresponding to g and the other to h
U = cell(2,par.n_eqn); dxU = U; dyU = U; UTemp = U; 
fM = cell(2,par.n_eqn);

% data structure for storing the values at the boundaries. The one coming
% from Sigma * B and the one coming from Sigma * g are both stored in
% bc_values. 
k_RK = cell(4,2);

for i = 1 : 2
    for j = 1:par.n_eqn
        fM{i,j} = X{1}*0;
    end
end

% prescribe the initial conditions
if isfield(par,'ic')
    U = capargs(par.ic,X{1},Y{1},par.system.Ax,par.system.Ay, ...
                    par.mass_matrix,par.inv_mass_matrix,par.value_f0);
end

%% Related to boundary treatment

% data structure for storing the values at the boundaries
bc_values = cell(2,par.n_eqn);

% we need to know which elements are coupled with which one at the
% boundaries. The 
% ID = 1
% a loop over all the boundaries
% consider the coupling for the Sigma * B term
bc_coupling_penalty_B = cell(par.num_bc,1);

% consider the coupling for the Sigma * g term
bc_coupling_penalty = cell(par.num_bc,1);

% we need to do this fixing for the machine error
for i = 1 : par.num_bc
    bc_coupling_penalty_B{i} = cellfun ( @(a) find(abs(a) > 1e-14),...
                                num2cell(par.system.penalty_B{i},2), 'Un', 0 );
    bc_coupling_penalty{i} = cellfun ( @(a) find(abs(a)> 1e-14),...
                                 num2cell(par.system.penalty{i},2), 'Un', 0 );
end

% scaling for the boundary conditions
bc_scaling = [1/PX{1}(1,1) 1/PY{1}(1,1) 1/PX{1}(1,1) 1/PY{1}(1,1)];


bc_g = cell(2,par.num_bc);
rhoW = cell(par.num_bc,1);

% compute the boundary inhomogeneity
t = 0;

% compute the boundary inhomogeneity at t = 0
for i = 1 : 2
    for j = 1:par.num_bc
       % a-priori computations can only be done for the inflow boundary
       if ~par.wall_boundary(j)
        % need to convert to cell for the computations which follow
        bc_g{i,j} = num2cell(capargs(par.bc_inhomo_inflow,par.system.B{j},j,par.system.Ax, ...
                                     par.system.Ay,i,U,par.all_w,t));
       end
        
    end
end

%% Time Loop
cputime = zeros(1,3);
t = 0; step_count = 0;
rho = par.compute_density(U,par.system.Ax,par.system.Ay,par.all_w);
initial_total_rho =integrate_xy(rho,PX{1},PY{1});

residual = 0;

while t < par.t_end || residual > 10^(-10)
   
    if ~par.steady_state
        if t+par.dt > par.t_end
            par.dt = par.t_end-t;
        end
    end
    
    % the ode sytem can be written as U_t = Op.
    % RK = 4 implementation
     tic
     switch par.RK_order
         case 2
             t_temp = [t t + par.dt/2];
             dt_temp = [0 par.dt/2];
             weight = [0 1];
         case 4
             t_temp = [t t + par.dt/2 t + par.dt/2 t + par.dt];
             dt_temp = [0 par.dt/2 par.dt/2 par.dt];
             weight = [1/6 2/6 2/6 1/6];
         otherwise
             assert(1==0, 'Should run with Runge Kutta order 2 or 4');
     end
     
     UTemp = U; 
     for RK = 1 : par.RK_order
         evaluate = time_dep & (t_temp(RK) > 0);
         
         for i = 1 : 2
             for j = 1:par.n_eqn
                 bc_values{i,j} = X{1} * 0;
             end
         end
         
         if evaluate(2)
             for i = 1 : 2
                 for j = 1:par.num_bc
                     if ~par.wall_boundary(j)
                     % need to convert to cell for the computations which follow
                        bc_g{i,j} = num2cell(capargs(par.bc_inhomo_inflow,par.system.B{j}, ...
                                            j,par.system.Ax,par.system.Ay,i,UTemp,par.all_w,t_temp(RK)));
                     end
                     
                 end
             end
         end
         
         rho = par.compute_density(UTemp,par.system.Ax,par.system.Ay,par.all_w);
         [ux,uy] = par.compute_velocity(UTemp,par.system.Ax,par.system.Ay,par.all_w);
         theta = par.compute_theta(UTemp,par.system.Ax,par.system.Ay,par.all_w);
        
         % we need to reconstruct the maxwellian which is independent of
         % the coming computation and therefore not influenced.
         fM = store_minimized_entropy(par.system.Ax,par.system.Ay, ...
                            par.mass_matrix,par.inv_mass_matrix,rho,ux,uy,theta,par.value_f0);
         
                      
         % compute the derivatives for g and h
         for i = 1 : 2
             for j = 1 : par.n_eqn
                 dxU{i,j} = DX{1} * UTemp{i,j};
                 dyU{i,j} = UTemp{i,j} * DY{1}; % DY is actually transpose of dy
             end
             
          % extract all the value at x = x_end.
              bc_ID = 1;
              values = cellfun(@(a) a(end,:),UTemp(i,:),'Un',0);
             
             if par.wall_boundary(bc_ID) % if the we have a gas wall then we recompute bc_g
                rhoW{bc_ID} = -sumcell(cellfun(@(a) a(end,:),...
                              UTemp(1,par.pos_U{bc_ID}),'Un',0),par.rhoW_vect{bc_ID})...
                             - par.compute_thetaW(bc_ID,t_temp(RK)) * par.rhoW_value{bc_ID};
                % for a wall bc_g{i,bc_ID} is a cell with every element being a vector. For the in
                %-flow boundary, every element of the cell is a scalar
                %quantity. This is because in case of wall, the density at
                %the wall changes along the entire boundary. 
                bc_g{i,bc_ID} = capargs(par.bc_inhomo_wall, par.system.B{bc_ID}, bc_ID,...
                        par.system.Ax, par.system.Ay, rhoW{bc_ID}, i,t_temp(RK));
             end
             
             for j = 1 : par.n_eqn
                 bc_values{i,j}(end,:) = bc_values{i,j}(end,:) + bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g{i,bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
             end
             
             
             
          bc_ID = 2;
          values = cellfun(@(a) a(:,end),UTemp(i,:),'Un',0);
          
          if par.wall_boundary(bc_ID) % if the we have a gas wall then we recompute bc_g
              rhoW{bc_ID} = -sumcell(cellfun(@(a) a(:,end),...
                  UTemp(1,par.pos_U{bc_ID}),'Un',0),par.rhoW_vect{bc_ID})...
                  - par.compute_thetaW(bc_ID,t_temp(RK)) * par.rhoW_value{bc_ID};
              
              bc_g{i,bc_ID} = capargs(par.bc_inhomo_wall, par.system.B{bc_ID}, bc_ID,...
                  par.system.Ax, par.system.Ay, rhoW{bc_ID}, i,t_temp(RK));
          end
          
          for j = 1 : par.n_eqn
              bc_values{i,j}(:,end) = bc_values{i,j}(:,end) + bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                  par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                  sumcell(bc_g{i,bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                  par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
          end


             bc_ID = 3;
             values = cellfun(@(a) a(1,:),UTemp(i,:),'Un',0);
             
             if par.wall_boundary(bc_ID) % if the we have a gas wall then we recompute bc_g
                rhoW{bc_ID} = -sumcell(cellfun(@(a) a(1,:),...
                              UTemp(1,par.pos_U{bc_ID}),'Un',0),par.rhoW_vect{bc_ID})...
                             - par.compute_thetaW(bc_ID,t_temp(RK)) * par.rhoW_value{bc_ID};
                bc_g{i,bc_ID} = capargs(par.bc_inhomo_wall, par.system.B{bc_ID}, bc_ID,...
                        par.system.Ax, par.system.Ay, rhoW{bc_ID}, i,t_temp(RK));
             end
             
             
             for j = 1 : par.n_eqn
                 bc_values{i,j}(1,:) = bc_values{i,j}(1,:) + bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g{i,bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})));
             end

             bc_ID = 4;
             values = cellfun(@(a) a(:,1),UTemp(i,:),'Un',0);
             
             if par.wall_boundary(bc_ID) % if the we have a gas wall then we recompute bc_g
                 rhoW{bc_ID} = -sumcell(cellfun(@(a) a(:,1),...
                     UTemp(1,par.pos_U{bc_ID}),'Un',0),par.rhoW_vect{bc_ID})...
                     - par.compute_thetaW(bc_ID,t_temp(RK)) * par.rhoW_value{bc_ID};
                 bc_g{i,bc_ID} = capargs(par.bc_inhomo_wall, par.system.B{bc_ID}, bc_ID,...
                     par.system.Ax, par.system.Ay, rhoW{bc_ID}, i,t_temp(RK));
             end
             
             for j = 1 : par.n_eqn
                 bc_values{i,j}(:,1)= bc_values{i,j}(:,1) + bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g{i,bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})));
             end
      
             for j = 1 : par.n_eqn
                 
                 
                 % assuming that the system is diagonal
                    W = - (dxU{i,j} * par.system.Ax(j,j) +  dyU{i,j} * par.system.Ay(j,j));
                    
                    k_RK{RK,i}{j} = W + bc_values{i,j};
             
                 k_RK{RK,i}{j} = k_RK{RK,i}{j} + ...
                                 (-UTemp{i,j}+fM{i,j})/par.Kn;
             end
         end
         
         % we should not update Utemp, we have not computed the two coupled
         % systems.
         if RK ~= par.RK_order
             for i = 1 : 2
                 for j = 1 : par.n_eqn
                     UTemp{i,j} = U{i,j} + k_RK{RK,i}{j} * dt_temp(RK + 1);
                 end
             end
         end
     end
      
    residual = 0;
    rho_before = par.compute_density(U,par.system.Ax,par.system.Ay,par.all_w);
    [ux_before,uy_before] = par.compute_velocity(U,par.system.Ax,par.system.Ay,par.all_w);
    theta_before = par.compute_theta(U,par.system.Ax,par.system.Ay,par.all_w);
    sigma_xx_before = par.compute_sigma_xx(U,par.system.Ax,par.system.Ay,par.all_w);
    sigma_xy_before = par.compute_sigma_xy(U,par.system.Ax,par.system.Ay,par.all_w);
    sigma_yy_before = par.compute_sigma_yy(U,par.system.Ax,par.system.Ay,par.all_w);
    qx_before = par.compute_qx(U,par.system.Ax,par.system.Ay,par.all_w);
    qy_before = par.compute_qy(U,par.system.Ax,par.system.Ay,par.all_w);
    
    for RK = 1 : par.RK_order
        for i = 1 : 2
            for j = 1 : par.n_eqn
                U{i,j} = U{i,j} + weight(RK) * k_RK{RK,i}{j} * par.dt;
            end
            
        end
    end
    
    rho_after = par.compute_density(U,par.system.Ax,par.system.Ay,par.all_w);
    [ux_after,uy_after] = par.compute_velocity(U,par.system.Ax,par.system.Ay,par.all_w);
    theta_after = par.compute_theta(U,par.system.Ax,par.system.Ay,par.all_w);
    sigma_xx_after = par.compute_sigma_xx(U,par.system.Ax,par.system.Ay,par.all_w);
    sigma_xy_after = par.compute_sigma_xy(U,par.system.Ax,par.system.Ay,par.all_w);
    sigma_yy_after = par.compute_sigma_yy(U,par.system.Ax,par.system.Ay,par.all_w);
    qx_after = par.compute_qx(U,par.system.Ax,par.system.Ay,par.all_w);
    qy_after = par.compute_qy(U,par.system.Ax,par.system.Ay,par.all_w);
    
    if par.steady_state
        residual = norm((rho_after-rho_before)/par.dt)^2 ...
                   +norm((ux_after-ux_before)/par.dt)^2 ...
                   +norm((uy_after-uy_before)/par.dt)^2 ...
                   +norm((theta_after-theta_before)/par.dt)^2 ...
                   +norm((sigma_xx_after-sigma_xx_before)/par.dt)^2 ...
                   +norm((sigma_xy_after-sigma_xy_before)/par.dt)^2 ...
                   +norm((sigma_yy_after-sigma_yy_before)/par.dt)^2 ...
                   +norm((qx_after-qx_before)/par.dt)^2 ...
                   +norm((qy_after-qy_before)/par.dt)^2;
    end
    
    residual = sqrt(residual);
        
    step_count = step_count + 1;
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;


    if mod(step_count,10) == 0
        disp('time: neqn: step_count:');
        disp(t);
        disp(par.n_eqn);
        disp(step_count);
                    
    % change in density
        disp('change in density');
        disp(integrate_xy(rho,PX{1},PY{1})-initial_total_rho);

        disp('residual');
        disp(residual);
    end
    
    if mod(step_count,500) == 0
        temp_residual = residual;
        par.write_solution(U,par,X{1},Y{1});
    end
    
    tic
    
    if par.t_plot
        
        var_plot = par.compute_theta(U,par.system.Ax,par.system.Ay,par.all_w);
        figure(1);
        surf(X{1},Y{1},var_plot), axis xy equal tight;
        
        title(sprintf('t = %0.2f',t));
        colorbar;
        xlabel('x'), ylabel('y');
        
        xlim(par.ax([1 2]));
        ylim(par.ax([3 4]));
        zlim([-0.1 1]);
        
        drawnow

    end
    
    % option to save the data during time step, required for convergence
    % studies
    if par.save_during && mod(step_count,100) == 0
        % we compute the norms of the different features of the solution
        capargs(par.compute_during,U,weight,k_RK,PX,DX,t);
    end
    
    cputime(2) = cputime(2) + toc;
end

fprintf('%0.0f time steps\n',step_count)           % Display test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%15.2fs%5.0f%%\n',... % and CPU times.
    'plotting:%16.2fs%5.0f%%\n'],cputime)
fprintf('residual while writting %0.15e:\n',temp_residual);
output = struct('X',X{1}, ...
                'Y',Y{1}, ...  
                'sol',U, ...
                 'PX',PX{1}, ...
                 'PY',PY{1}, ...
                 'h',h);
end


function z = capargs(fct,varargin)
% Call function fct with as many arguments as it requires (at least 1),
% and ignore further arguments.
narg = max(nargin(fct),1);
z = fct(varargin{1:narg});

end

function f = zero(varargin)
% Zero function.
f = zeros(size(varargin{1}));
end

function S = sumcell(A,w)
% Add vector of cells A, weighted with vector w.
S = 0;

for j = 1:length(w)
    S = S+A{j}*w(j);
end

end


function f = integrate_xy(data,PX,PY)
diag_PX = diag(PX);
diag_PY = diag(PY);

% integral along x at different y points
int_x = zeros(size(data,2),1);

for i = 1 : length(int_x)
    int_x(i) = dot(diag_PX,data(:,i));
end

% integrate along the y direction
f = full(dot(int_x,diag_PY));
end
