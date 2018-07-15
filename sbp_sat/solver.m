function output = solver(par)

%% Few Checks

% if no initial conditions then zero
if ~isfield(par,'initial_condition')
    par.initial_condition = @zero; % Default: no init. cond.
end
if ~isfield(par,'source')
    par.source = @zero; % Default: no source.
end
if par.num_bc ~=4
    assert(1 == 0, 'not valid num bc');
end

% find which of the given data is time dependent
time_dep = [nargin(par.source)>3 nargin(par.bc_inhomo)>4];


%% Storing index of non-zero entries in system matrices (separately for
% each equation)
% corresponding to every row in Ax, stores the non-zero indices
Ix = cellfun(@find,num2cell(par.system.Ax',1),'Un',0);
Iy = cellfun(@find,num2cell(par.system.Ay',1),'Un',0);
% each row of Ax in separate cell as here num2cell(Ax',1)
% what does 'Un',0 do? and purpose of Ix, Iy? --> now clear
% Error using cellfun
% Non-scalar in Uniform output, at index 1, output 1.
% Set 'UniformOutput' to false.


%% Grid Setup

% size of the grid. h1 for the x direction and h2 for the y direction
h = (par.ax([2 4])-par.ax([1 3]))./par.n; % ./ does right-array division

% Set of grid nodes in x and y direction
x{1} = par.ax(1):h(1):par.ax(2);
y{1} = par.ax(3):h(2):par.ax(4);

% Generate both grids. X contains the x coordinates and Y contains
% the y coordinates.
[X,Y] = cellfun(@ndgrid,[x],[y],'Un',0);


%% Time-step

% a crude approximation for delta_t
par.dt = par.CFL * min(h)/abs(eigs(par.system.Ax,1,'lm')); % 1 eigenv.
% with largest magnitude


%% Setup difference operators

switch par.diff_order
    case 2
        % differential opeartor for the X grid, i.e. dx matrix
        [DX{1},PX{1}] = sbp_collocated_2(par.n(1),h(1));
        
        % differential operators for the Y grid, i.e. dy matrix
        [DY{1},PY{1}] = sbp_collocated_2(par.n(2),h(2));
        
    case 4
        assert(par.diff_order==2, ...
            '4th order difference operator not implemented yet');
    otherwise
        assert(1==0, ...
            'Should run with 2nd or 4th order difference operator');
end

% we transpose DY for the data structure
DY = cellfun(@transpose,DY,'Un',0);


%% Initialize the solution variables

% dxU and dyU are derivatives in x and y direction
% force is forcing term
% set initial conditions else initialize with zero
U = cell(1,par.n_eqn); dxU = U; dyU = U; UTemp = U;

k_RK = cell(4,1); % to store the 4 slopes- k1, k2, k3 and k4, i.e. L(U,t)

for j = 1:par.n_eqn
    U{j} = X{1}*0 + capargs(par.initial_condition,X{1},Y{1},j);
end

t = 0; % setting the current time


%% Related to boundary treatment

% data structure for storing the values at the boundaries
bc_values = cell(1,par.n_eqn);
for j = 1:par.n_eqn
    bc_values{j} = X{1}*0;
end

% scaling for the boundary conditions
bc_scaling = [1/PX{1}(1,1) 1/PY{1}(1,1) 1/PX{1}(1,1) 1/PY{1}(1,1)];

% boundary inhomogeneity
bc_g = cell(par.num_bc,1);


%% Forcing Term

force = U;

for j = 1:par.n_eqn
    force{j} = X{1}*0;
end

for j = 1:par.n_eqn
    force{j} = capargs(par.source,X{1},Y{1},j,t);
end


%% Time Loop
par.t_plot = [par.t_plot par.t_end inf]; plot_count = 1;
cputime = zeros(1,3);
step_count = 0;

while t < par.t_end
    
    if t+par.dt > par.t_end
        par.dt = par.t_end-t;
    end
    
    % Runge Kutta time-stepping (RK2 and RK4)
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
        
        if evaluate(1)
            for j = 1:par.n_eqn
                force{j} = capargs(par.source,X{1},Y{1},j,t_temp(RK));
            end
        end
        
%         if evaluate(2)
%             % COMPUTING bc_g every time-step in the code below
%         end
        
        for i = 1:par.n_eqn
            dxU{i} = DX{1} * UTemp{i};
            dyU{i} = UTemp{i} * DY{1}; % DY is actually transpose of dy
        end
        
        %
        % extract all the values at x = 1, last row of the matrix.
        % id = 1, x = 1, last row of X
        % id = 2, y = 1, last column of Y
        % id = 3, x = 0, first row of X
        % id = 4, y = 0, first column of Y
        
        for bc_ID = [2 4 1 3]
            if rem(bc_ID,2)==0
                b_nodes = par.n(1)+1 ;
            else
                b_nodes = par.n(2)+1 ;
            end
            
            for i=1:b_nodes
                if bc_ID==2
                    y_coord = 1 ;
                elseif bc_ID==4
                    y_coord = 0 ;   
                else
                    y_coord = y{1}(i) ;
                end
                switch bc_ID
                    case 1
                        values = cellfun(@(a) a(end,i),UTemp,'Un',0);
                    case 2
                        values = cellfun(@(a) a(i,end),UTemp,'Un',0);
                    case 3
                        values = cellfun(@(a) a(1,i),UTemp,'Un',0);
                    case 4
                        values = cellfun(@(a) a(i,1),UTemp,'Un',0);
                end
                
                par.sys.B = par.get_boundary_operator(par.system.nx, par.system.ny, bc_ID, y_coord) ;
                
                switch par.penalty_id
                    case 1
                        par.sys.penalty = par.get_penalty_1(par.system.A_n{bc_ID}, par.sys.B, par.system.nx, par.system.ny, bc_ID) ;
                    case 2
                        par.sys.penalty = par.get_penalty_2(par.system.nx, par.system.ny, bc_ID) ;
                end
                
                par.sys.penalty_B = par.sys.penalty * par.sys.B ;
                
                bc_g_2{1} = (capargs(par.bc_inhomo,par.sys.B,bc_ID,y_coord,par.n_eqn,t_temp(RK)));
                
                for j = 1 : par.n_eqn
                    temp_copy = bc_scaling(bc_ID) * ( sumcell(values,...
                        par.sys.penalty_B(j,:)) - ...
                        sumcell(bc_g_2(1), par.sys.penalty(j,:)) );
                    switch bc_ID
                    case 1
                        bc_values{j}(end,i) = temp_copy;
                    case 2
                        bc_values{j}(i,end) = temp_copy;
                    case 3
                        bc_values{j}(1,i) = temp_copy;
                    case 4
                        bc_values{j}(i,1) = temp_copy;
                end
                end
            end
            
        end
        %
        
        for i = 1 : par.n_eqn
            % Multiplication of dxU and dyU by system matrices Ax and Ay
            % [1] dxU(Ix{i}) returns a cell with components of dxU where Ax
            % has non-zero entries for that particular row (component no.)
            % --> each dxU entry is a cell as for each component the values
            % need to be stored at all grid nodes
            % [2] For Ax(i,Ix{i}), lets say i = 1 (component/equ. 1) then
            % Ix{1} is [3 4 5] non-zero entries in 1st row of Ax.
            % Ax(1, [3 4 5]) gives (1,3), (1,4) and (1, 5) entries of Ax
            % [3] W is one vector having a value stored at each grid node
            W = -sumcell([dxU(Ix{i}),dyU(Iy{i})],...
                [par.system.Ax(i,Ix{i}),par.system.Ay(i,Iy{i})]);
            
            % for each RK stage k_RK{RK} contains 1 x n_equ sized cell
            k_RK{RK}{i} = (W  + force{i} + bc_values{i});
            
            if RK ~= par.RK_order
                UTemp{i} = U{i} + k_RK{RK}{i} * dt_temp(RK + 1);
            end
        end
    end
    
    for RK = 1 : par.RK_order
        for i = 1 : par.n_eqn
            U{i} = U{i} + weight(RK) * k_RK{RK}{i} * par.dt;
        end
    end
    
    step_count = step_count + 1;
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;
    
    %% Plotting
    
    % for plotting (new)
    if mod(step_count,50) == 0
        disp('time: neqn: step_count ');
        disp(t);
        disp(par.n_eqn);
        disp(step_count);
    end
    
    if par.to_plot
        
        %         figure(1)
        % %          surf(X{1},Y{1},U{par.var_plot}), axis xy equal tight;
        %         contourf(X{1},Y{1},U{par.var_plot}), axis xy equal tight;
        %
        %         % get rid of lines in surf
        %         colormap summer;
        %         shading interp;
        %         title(sprintf('t = %0.2f',t));
        % %         colorbar;
        %         xlabel('x'), ylabel('y');
        %
        %
        %         xlim(par.ax([1 2]));
        %         ylim(par.ax([3 4]));
        %         zlim([-0.75 0.75]);
        %
        %         drawnow
        
        % plot 2
        
        
        figure(2)
        plot(x{1}, U{1}(:,2), x{1}, sin(4*pi*x{1})*cos(4*pi*sqrt(2)*t), 'k--') ;
        
        title(sprintf('t = %0.2f',t));
        %         colorbar;
        xlabel('x'), ylabel('y');
        
        
        xlim(par.ax([1 2]));
        ylim([-1 1]) ;
        %         ylim(par.ax([3 4]));
        %         zlim([-0.75 0.75]);
        
        drawnow
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
    'PX',PX, ...
    'PY',PY, ...
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
% Add vector of cells A, weighted with (corresponding entries of) vector w.
S = 0; % *S is a vector as A{1} is vector*. Need this separately for
% initialization of S with proper length? NOPE
% A{1} can also be a matrix, right? then mat. vect. multip. *WRONG* A all
% the entries of matrix but A{i} has non-zero ith column entries
for j = 1:length(w)
    S = S + A{j}*w(j);
end
end

function vector = sumcell_2(A,w,var,nodes,n_eqn)
vector = [] ;
for i=1:nodes
    x{i} = [];
end
for i=1:nodes
    for j=1:n_eqn
        x{i} = [x{i} A{j}(i)] ;
        y{i} = w{i}(var,:) ;
    end
end
for i=1:nodes
    vector = [vector sum(x{i} .* y{i})] ;
    %     vector = vector + x{i} .* y{i} ;
end
end

