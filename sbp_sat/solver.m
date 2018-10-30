function output = solver(par)

%% Few Checks

% if no initial conditions then zero
if ~isfield(par,'initial_condition')
    par.initial_condition = @zero; % Default: no init. cond.
end

if ~isfield(par,'steady_state')
    par.steady_state = false; % Default value for steady state
end

if ~isfield(par,'source')
    par.source = @zero; % Default: no source.
end 

if par.num_bc ~=4
    assert(1 == 0, 'not valid num bc');
end

% % find which of the given data is time dependent
time_dep = [nargin(par.source)>3 nargin(par.bc_inhomo)>2];


%% Storing index of non-zero entries in system matrices (separately for
                                                           % each equation)
% corresponding to every row in Ax, stores the non-zero indices
Ix = cellfun(@find,num2cell(par.system.Ax',1),'Un',0);
Iy = cellfun(@find,num2cell(par.system.Ay',1),'Un',0);
Ix_prod = cellfun(@find,num2cell(par.system.P,1),'Un',0);  

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


[X,Y] = cellfun(@ndgrid,x,y,'Un',0);


%% Time-step

% a crude approximation for delta_t
par.dt = min(h) * par.CFL/(2 * abs(eigs(par.system.Ax,1,'lm'))); % 1 eigenv.
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
    U{j} = X{1}*0 + capargs(par.initial_condition,X{1},Y{1},j,par.previous_M_data,par.n_eqn);
    
    force{j} = X{1}*0;
end

%% Related to boundary treatment

% data structure for storing the values at the boundaries
bc_values = cell(par.num_bc,par.n_eqn);

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
    % for the boundary id=i, we find the variables to which every variable
    % in the moment system
    % is coupled 
    % the matrix penalty_B{i} contains Sigma * B at the boundary with ID =
    % i. the matrix penalty{i} contains the penalty matrix at the boundary
    % with ID = i.
    bc_coupling_penalty_B{i} = cellfun ( @(a) find(abs(a) > 1e-14), num2cell(par.system.penalty_B{i},2), 'Un', 0 );
        % num2cell(par.system.penalty_B{i},2) splits B{i} into 
        % separate cells where dim specifies which dimensions of A to 
        % include in each cell
    bc_coupling_penalty{i} = cellfun ( @(a) find(abs(a)> 1e-14), num2cell(par.system.penalty{i},2), 'Un', 0 );
end

% scaling for the boundary conditions
bc_scaling = [1/PX{1}(1,1) 1/PY{1}(1,1) 1/PX{1}(1,1) 1/PY{1}(1,1)];

%the boundary inhomogeneity. Every component of the cell will be equal to
%the number of boundar conditions which we need to prescribe. 
bc_g = cell(par.num_bc,1);

% compute the boundary inhomogeneity
t = 0;
for j = 1:par.num_bc
    % need to convert to cell for the computations which follow
    bc_g{j} = num2cell(capargs(par.bc_inhomo,par.system.B{j},j,t));
end

for j = 1:par.n_eqn
    force{j} = capargs(par.source,X{1},Y{1},j,t);
end

%% Time Loop
cputime = zeros(1,3);
step_count = 0;
residual = 0;
plot_counter = 0;

while t < par.t_end || residual > 10^(-8)
    
    if t+par.dt > par.t_end && ~par.steady_state
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
        
       % need to reinitialize
       for i = 1 : par.num_bc
            for j = 1:par.n_eqn
                bc_values{i,j} = X{1}*0;
            end
       end
        
        evaluate = time_dep & (t_temp(RK) > 0);
            
            if evaluate(1)
                for j = 1:par.n_eqn 
                         force{j} = capargs(par.source,X{1},Y{1},j,t_temp(RK)); 
                end
            end
            
            if evaluate(2)
                for j = 1:par.num_bc
                 % need to convert to cell for the computations which follow
                    bc_g{j} = num2cell(capargs(par.bc_inhomo,par.system.B{j},j,t_temp(RK)));
   
                end
            end
            
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
        bc_ID = 1;
        values = cellfun(@(a) a(end,:),UTemp,'Un',0);
        for j = 1 : par.n_eqn
%                 the term, values(bc_coupling_penalty_B{bc_ID}{j}), gives us the
%                 value of all the variables, at the boundary, which are
%                 coupled with the j-th variable. 
%                 par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j})
%                 gives us the j-th row of the penalty matrix and the
%                 entries in all those columns which have no zeros.
            bc_values{j}(end,:) = bc_values{j}(end,:) + bc_scaling(bc_ID) * ( sumcell( values(bc_coupling_penalty_B{bc_ID}{j}),...
                                  par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j}) ) - ...
                                  sumcell(bc_g{bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                                  par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
        end
                    
        bc_ID = 2;
        values = cellfun(@(a) a(:,end),UTemp,'Un',0);
        for j = 1 : par.n_eqn
            bc_values{j}(:,end) = bc_values{j}(:,end) + bc_scaling(bc_ID) * ( sumcell(values(bc_coupling_penalty_B{bc_ID}{j}), ...
                                  par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j})) - ...
                                  sumcell(bc_g{bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                                  par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
        end
%        

        bc_ID = 3;
        values = cellfun(@(a) a(1,:),UTemp,'Un',0);
        for j = 1 : par.n_eqn
            bc_values{j}(1,:) = bc_values{j}(1,:) + bc_scaling(bc_ID) * ( sumcell(values(bc_coupling_penalty_B{bc_ID}{j}), ...
                                par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j})) - ...
                                sumcell(bc_g{bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                                par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
        end
        

        
        bc_ID = 4;
        values = cellfun(@(a) a(:,1),UTemp,'Un',0);
        for j = 1 : par.n_eqn
            bc_values{j}(:,1) = bc_values{j}(:,1) + bc_scaling(bc_ID) * ( sumcell(values(bc_coupling_penalty_B{bc_ID}{j}), ...
                                par.system.penalty_B{bc_ID}(j,bc_coupling_penalty_B{bc_ID}{j})) - ...
                                sumcell(bc_g{bc_ID}(bc_coupling_penalty{bc_ID}{j}),...
                                par.system.penalty{bc_ID}(j,bc_coupling_penalty{bc_ID}{j})) );
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
            
            % add contribution from all the boundaries
%             for j = 1 : par.num_bc
%                 k_RK{RK}{i} = k_RK{RK}{i} + bc_values{j,i};
%             end
            
            k_RK{RK}{i} = k_RK{RK}{i} ...
                           + sumcell(UTemp(Ix_prod{i}),par.system.P(i,Ix_prod{i}))/par.Kn;
        end
        
        for i = 1 : par.n_eqn
            if RK ~= par.RK_order
                UTemp{i} = U{i} + k_RK{RK}{i} * dt_temp(RK + 1);
            end
        end
    end
    
    residual = 0; % need to reinitialise
    for RK = 1 : par.RK_order
        for i = 1 : par.n_eqn
            U{i} = U{i} + weight(RK) * k_RK{RK}{i} * par.dt;
            
                if par.steady_state && i <= 13 % only check the residual in macroscopic quantities 
                    residual = residual + norm(weight(RK) * k_RK{RK}{i})^2;
                end
        end
    end
    residual = sqrt(residual);
    
    t = t + par.dt;
    cputime(1) = cputime(1) + toc;
    
    if mod(step_count,100) == 0
        disp('time: neqn: step_count: ');
        disp(t);
        disp(par.n_eqn);
        disp(step_count);
        disp('residual: ');
        disp(residual);
    end
    
%     if mod(step_count,500) == 0
%         temp_residual = residual;
%         par.write_solution(U,par,X{1},Y{1},par.M,[t,temp_residual]);
%     end
    
    %% Plotting
    if par.to_plot && mod(step_count,10) == 0
        
        surface_plot = surf(X{1},Y{1},par.compute_theta(U)), axis xy equal tight;
       
       
        title(sprintf('t = %0.2f',t));
        colorbar;
        xlabel('x'), ylabel('y');
        
        xlim(par.ax([1 2]));
        ylim(par.ax([3 4]));
        zlim([-0.2 0.5]);
        
        plot_counter = plot_counter + 1;
        filename_figure = strcat('/Users/neerajsarna/Dropbox/my_papers/MPI_ppt/pictures/odd_bc-',...
                                 num2str(plot_counter),'.png');                     
        saveas(surface_plot,filename_figure);                      
        drawnow;
    end
    
    step_count = step_count + 1;
end

fprintf('%0.0f time steps\n',step_count)           % Display test
cputime = reshape([cputime;cputime/sum(cputime)*1e2],1,[]);   % case info
fprintf(['CPU-times\n advection:%15.2fs%5.0f%%\n',... % and CPU times.
    'plotting:%16.2fs%5.0f%%\n'],cputime)
fprintf('final residual while writting: %0.15e\n',temp_residual);

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
% Add vector of cells A, weighted with (corresponding entries of) vector w.
S = 0; % *S is a vector as A{1} is vector*. Need this separately for
% initialization of S with proper length? NOPE
% A{1} can also be a matrix, right? then mat. vect. multip. *WRONG* A all
% the entries of matrix but A{i} has non-zero ith column entries
for j = 1:length(w)
    S = S + A{j}*w(j);
end
end

