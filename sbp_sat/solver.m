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

% % find which of the given data is time dependent
% time_dep = [nargin(par.source)>2 0];


%% Storing index of non-zero entries in system matrices (separately for
                                                           % each equation)
% corresponding to every row in Ax, stores the non-zero indices
Ix = cellfun(@find,num2cell(par.system_data.Ax',1),'Un',0);
Iy = cellfun(@find,num2cell(par.system_data.Ay',1),'Un',0);
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
par.dt = par.CFL * min(h)/abs(eigs(par.system_data.Ax,1,'lm')); % 1 eigenv.
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
U = cell(1,par.n_eqn); dxU = U; dyU = U; force = U; UTemp = U;

k_RK = cell(4,1); % to store the 4 slopes- k1, k2, k3 and k4, i.e. L(U,t)

for j = 1:par.n_eqn
    U{j} = X{1}*0 + capargs(par.initial_condition,X{1},Y{1},j);
    
    force{j} = X{1}*0;
    bc_values{j} = X{1}*0;
end

%% Time Loop
par.t_plot = [par.t_plot par.t_end inf]; plot_count = 1;
cputime = zeros(1,3);
t = 0; step_count = 0;

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
        for i = 1:par.n_eqn
            dxU{i} = DX{1} * UTemp{i};
            dyU{i} = UTemp{i} * DY{1}; % DY is actually transpose of dy
        end
        
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
                [par.system_data.Ax(i,Ix{i}),par.system_data.Ay(i,Iy{i})]);
                        
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
    tic
    while t>=par.t_plot(plot_count)-1e-14  % If current time has exceeded
        lambda = (par.t_plot(plot_count)-t+par.dt)/par.dt; %plotting time, define
        Uplot = U{par.mom_output}; % linear interpolation in time.
        if plot_count==length(par.t_plot)-1      % Is final time reached?
            if nargout                   % If yes, save solution at final
                output = struct('x',x{1},'y',y{1},'U',U,'Px',PX{1},'Py',PY{1});     % time.
            end
        else                           % If not, invoke plotting routine.
            xplot = x{1};             % Assign grids at
            yplot = y{1};          % outputted moments.
            %             if length(par.mom_output)==1        % If only a single moment
            %                 xplot = xplot{:}; yplot = yplot{:}; Uplot = Uplot{par.mom_output};% is
            %             end                         % plotted, remove cell structure.
            if nargout(par.output)                     % Call output rou,
                par = par.output(par,xplot,yplot,Uplot,plot_count);%tine-
            else                          % allowing for it to modify the
                par.output(par,xplot,yplot,Uplot,plot_count)% struct par.
            end
        end
        plot_count = plot_count+1;
    end
    cputime(2) = cputime(2) + toc;
end

fprintf('%0.0f time steps\n',step_count)           % Display test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%15.2fs%5.0f%%\n',... % and CPU times.
    'plotting:%16.2fs%5.0f%%\n'],cputime)

U_theo = cell(1,par.n_eqn);
error = cell(1,par.n_eqn);
for j = 1:par.n_eqn
    U_theo{j} = X{1}*0 + capargs(par.theoretical_solution,X{1},Y{1},t,j);
    error{j} = U{j} - U_theo{j};
end
[output.error] = error{:};

error_l2 = zeros(par.n_eqn,1);
for j = 1:par.n_eqn
    error_l2(j) = norm(error{j},1);
end
error_l2 = num2cell(error_l2);
[output.error_l2] = error_l2{:};


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
S = A{1}*w(1); % S is a vector as A{1} is vector. Need this separately for
% initialization of S with proper length? NOPE
for j = 2:length(w)
    S = S+A{j}*w(j);
end
end

