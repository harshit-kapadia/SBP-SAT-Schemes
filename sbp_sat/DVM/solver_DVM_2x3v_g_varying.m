function output = solver_DVM_2x3v_g_varying(par)

% default value of save during the computation
if ~isfield(par,'save_during'), par.save_during = false; end 

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
time_dep = [nargin(par.source)>2 nargin(par.bc_inhomo)>8];

% if the number of arguments are greater than 3 then definitely we have
% an anisotropic source term.
if ~isfield(par,'source_ind')
    par.source_ind = 1:par.n_eqn;
end

% corresponding to every row in Ax, stores the non-zero indices
Ix = cellfun(@find,num2cell(par.system.Ax',1),'Un',0);
Iy = cellfun(@find,num2cell(par.system.Ay',1),'Un',0);

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

for i = 1 : 2
    for j = 1:par.n_eqn
        bc_values{i,j} = X{1}*0;
    end
end

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


%% Time Loop
cputime = zeros(1,3);
t = 0; step_count = 0;


residual = 0;


while t < par.t_end || residual > 10^(-6)
    
    residual = 0;
    
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
        %         evaluate = time_dep & (t_temp(RK) > 0);
        %
        %         % no need to compute bc_g
        %         if evaluate(2)
        %             % no need to compute bc_g
        %         end
        
        rho = par.compute_density(UTemp,par.system.Ax,par.system.Ay,par.all_w);
        [ux,uy] = par.compute_velocity(UTemp,par.system.Ax,par.system.Ay,par.all_w);
        theta = par.compute_theta(UTemp,par.system.Ax,par.system.Ay,par.all_w);
        
        % we need to reconstruct the maxwellian which is independent of
        % the coming computation and therefore not influenced.
        fM = store_minimized_entropy(par.system.Ax,par.system.Ay, ...
            par.mass_matrix,par.inv_mass_matrix,rho,ux,uy,theta,par.value_f0);
        
        
        % Temperature at the wall, needed for Wall boundary implementation
        if t_temp(RK) <= 1
            temp_thetaW = exp(-1/(1-(t-1)^2)) * exp(1);
        else
            temp_thetaW = 1;
        end
        thetaW{1} = temp_thetaW ;
        thetaW{2} = temp_thetaW ;
        thetaW{3} = -temp_thetaW ;
        thetaW{4} = -temp_thetaW ;
        
%         par.b_nodes = zeros(4,1) ;
%         for bc_ID = [1 2 3 4]
%             if rem(bc_ID,2)==0
%                 par.b_nodes(bc_ID) = par.n(1)+1 ;
%             else
%                 par.b_nodes(bc_ID) = par.n(2)+1 ;
%             end
%         end
        % compute the derivatives for g and h
        for i = 1 : 2
            for j = 1 : par.n_eqn
                dxU{i,j} = DX{1} * UTemp{i,j};
                dyU{i,j} = UTemp{i,j} * DY{1}; % DY is actually transpose of dy
            end
            
            
%             for bc_ID = [1 2 3 4]
%                 for node = 1:par.b_nodes(bc_ID)
%                     switch bc_ID
%                         case 1
%                             values = cellfun(@(a) a(end,node),UTemp(i,:),'Un',0);
%                         case 2
%                             values = cellfun(@(a) a(node,end),UTemp(i,:),'Un',0);
%                         case 3
%                             values = cellfun(@(a) a(1,node),UTemp(i,:),'Un',0);
%                         case 4
%                             values = cellfun(@(a) a(node,1),UTemp(i,:),'Un',0);
%                     end
%                     
%                     bc_g = num2cell(capargs(par.bc_inhomo,...
%                         par.system.B{bc_ID}, bc_ID,par.system.Ax,...
%                         par.system.Ay,i,UTemp,node,par.all_w,t_temp(RK)));
%                     
%                     
%                     for j = 1 : par.n_eqn
%                         temp_copy = bc_scaling(bc_ID) * ( sumcell( ...
%                             values(bc_coupling_penalty_B{bc_ID}{j}), ...
%                             par.system.penalty_B{bc_ID}(j,...
%                                 bc_coupling_penalty_B{bc_ID}{j}) ) - ...
%                             sumcell(bc_g(bc_coupling_penalty{bc_ID}{j}),...
%                             par.system.penalty{bc_ID}(j,...
%                                 bc_coupling_penalty{bc_ID}{j})) );
%                         switch bc_ID
%                             case 1
%                                 bc_values{i,j}(end,node) = temp_copy;
%                             case 2
%                                 bc_values{i,j}(node,end) = temp_copy;
%                             case 3
%                                 bc_values{i,j}(1,node) = temp_copy;
%                             case 4
%                                 bc_values{i,j}(node,1) = temp_copy;
%                         end
%                     end
%                     
%                 end
%             end
            
            

            % extract all the value at x = x_start.
            bc_ID = 1;
            values = cellfun(@(a) a(end,:),UTemp(i,:),'Un',0);
            rhoW{bc_ID} = -cell2mat(cellfun(@(a) a(end,:)',...
                UTemp(1,par.pos_U{bc_ID}),'Un',0)) * par.rhoW_vect{bc_ID} ...
                - thetaW{bc_ID} * par.rhoW_value{bc_ID} ;
            rhoW{bc_ID} = rhoW{bc_ID}' ;
            bc_g = capargs(par.bc_inhomo, par.system.B{bc_ID}, bc_ID,...
                par.system.Ax, par.system.Ay, thetaW{bc_ID}, rhoW{bc_ID}, i);
            
             for j = 1 : par.n_eqn
                 bc_values{i,j}(end,:) = bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
             end
             
             
             bc_ID = 2;
             values = cellfun(@(a) a(:,end),UTemp(i,:),'Un',0);
             rhoW{bc_ID} = -cell2mat(cellfun(@(a) a(:,end),...
                 UTemp(1,par.pos_U{bc_ID}),'Un',0)) * par.rhoW_vect{bc_ID} ...
                 - thetaW{bc_ID} * par.rhoW_value{bc_ID} ;
             bc_g = capargs(par.bc_inhomo, par.system.B{bc_ID}, bc_ID,...
                 par.system.Ax, par.system.Ay, thetaW{bc_ID}, rhoW{bc_ID}, i);
             
             for j = 1 : par.n_eqn
                 bc_values{i,j}(:,end) = bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
             end
             
             bc_ID = 3;
             values = cellfun(@(a) a(1,:),UTemp(i,:),'Un',0);
             rhoW{bc_ID} = -cell2mat(cellfun(@(a) a(1,:)',...
                 UTemp(1,par.pos_U{bc_ID}),'Un',0)) * par.rhoW_vect{bc_ID} ...
                 - thetaW{bc_ID} * par.rhoW_value{bc_ID} ;
             rhoW{bc_ID} = rhoW{bc_ID}' ;
             bc_g = capargs(par.bc_inhomo, par.system.B{bc_ID}, bc_ID,...
                 par.system.Ax, par.system.Ay, thetaW{bc_ID}, rhoW{bc_ID}, i);
             
             for j = 1 : par.n_eqn
                 bc_values{i,j}(1,:) = bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})));
             end
             
             bc_ID = 4;
             values = cellfun(@(a) a(:,1),UTemp(i,:),'Un',0);
             rhoW{bc_ID} = -cell2mat(cellfun(@(a) a(:,1),...
                 UTemp(1,par.pos_U{bc_ID}),'Un',0)) * par.rhoW_vect{bc_ID} ...
                 - thetaW{bc_ID} * par.rhoW_value{bc_ID} ;
             bc_g = capargs(par.bc_inhomo, par.system.B{bc_ID}, bc_ID,...
                 par.system.Ax, par.system.Ay, thetaW{bc_ID}, rhoW{bc_ID}, i);
             
             for j = 1 : par.n_eqn
                 bc_values{i,j}(:,1)= bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                     par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                     sumcell(bc_g(bc_coupling_penalty{bc_ID}{j}),...
                     par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})));
             end


            
            for j = 1 : par.n_eqn
                
                % multiplication of the derivatives and the system matrices
                W = -sumcell([dxU(i,Ix{j}),dyU(i,Iy{j})],...
                    [par.system.Ax(j,Ix{j}),par.system.Ay(j,Iy{j})]);
                
                k_RK{RK,i}{j} = (W + bc_values{i,j});
                
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
    
    for RK = 1 : par.RK_order
        for i = 1 : 2
            for j = 1 : par.n_eqn
                U{i,j} = U{i,j} + weight(RK) * k_RK{RK,i}{j} * par.dt;
                
                % if we are computing for the steady state
                if par.steady_state
                    residual = residual + norm(weight(RK) * k_RK{RK,i}{j})^2;
                end
            end
        end
    end
    
    % earlier we computed the square of the norm
    residual = sqrt(residual);
    
    step_count = step_count + 1;
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;
    
    
    if mod(step_count,50) == 0
        disp('time: neqn: step_count: residual: ');
        disp(t);
        disp(par.n_eqn);
        disp(step_count);
        disp(residual);
    end
    
    tic
    
    if par.t_plot
        var_plot = par.compute_density(U,par.system.Ax,par.system.Ay,par.all_w);
        
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



