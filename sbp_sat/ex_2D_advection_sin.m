clear all

%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','advection equation',... % name of example
    'n_eqn',1,... % number of eq per grid point
    'initial_condition',@initial_condition,... % it is defined below
    'exact_solution',@exact_solution,...
    'ax',[0 1 0 1],... % extents of computational domain
    'n',[100 100],... % numbers of grid cells in each coordinate direction
    't_end',1.0,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',4,...
    'CFL',2,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'var_plot',1,...
    'to_plot',false,...
    'source',@source... % source function to satisfy the manufactured solution
    );

[par.system] = get_system(par.n_eqn);

% we need the boundary matrix and the penalty matrix for all the
% boundaries
par.system.penalty_B = cell(par.num_bc,1);
par.system.penalty = cell(par.num_bc,1);
par.system.B = cell(par.num_bc,1);

% id =1, x = 1
% id = 2 , y = 1
% id = 3, x = 0
% id = 4, y = 0
par.system.B{1} = zeros(par.n_eqn) ;
par.system.B{2} = zeros(par.n_eqn) ;
par.system.B{3} =  eye(par.n_eqn) ;
par.system.B{4} = eye(par.n_eqn) ;

% we now compute the penalty matrix
par.normals_bc = [1,0;0,1;-1,0;0,-1];
for i = 1 : par.num_bc
    An = par.system.Ax * par.normals_bc(i,1) + par.system.Ay * par.normals_bc(i,2);
    par.system.penalty{i} = dvlp_penalty_char(An,par.system.B{i});
end

for i = 1 : par.num_bc
    par.system.penalty_B{i} = par.system.penalty{i}*par.system.B{i};
end

%========================================================================
% Run solver and study convergence
%========================================================================
resolution = [16 32 64 128 256];
grid_spacing = [];

for k = 1:length(resolution)                       % Loop over various grid resolutions.
    par.n = [1 1]*resolution(k);                   % Numbers of grid cells.
    solution = solver(par);                         % Run solver.
    error_temp = 0;
    disp('Resolution :');
    disp(par.n);
    for j = 1:par.n_eqn                               % Loop over solution components.
        X = solution(j).X;
        Y = solution(j).Y;
        PX = solution(j).PX;
        PY = solution(j).PY;
        U_theo = exact_solution(X,Y,par.t_end);     % Evaluate true solution.
        error = abs(solution(j).sol-U_theo);               % Difference between num. and true sol.
        int_x = dot(transpose(error),transpose(PX * error),2);  % integral along x.
        int_xy = sum(PY*int_x);  % integral along xy.
        error_temp = int_xy + error_temp;%Sc. L2 error.
    end
    error_L2(k) = sqrt(error_temp);
    grid_spacing = [grid_spacing max(solution(1).h)];
end

reference_line = exact_order(grid_spacing(1),error_L2(1),grid_spacing(end),2);

figure
loglog(grid_spacing, error_L2, '-o',reference_line(1:2),reference_line(3:4),'-*');
xlabel('h'), ylabel('l2-error');
legend('numerical','second order');
title('Convergence plot');

convg_order = log(error_L2(end)/error_L2(1))/log(grid_spacing(end)/grid_spacing(1));
disp('convergence order');
disp(convg_order);

%========================================================================
% Problem Specific Functions
%========================================================================

function[system] = get_system(n_equ)

value_x = ones(n_equ);
system.Ax = spdiags([value_x], [0], n_equ, n_equ);

value_y = ones(n_equ);
system.Ay = spdiags([value_y], [0], n_equ, n_equ);

end

function f = initial_condition(x,y)
f = sin(pi * x).* sin(pi * y);
end

function f = exact_solution(x,y,t)
f = sin(pi * x).* sin(pi * y) * cos(pi * t);
end

function f = bc_inhomo(B,bc_id)
    switch bc_id
        % east boundary but in matrix - south, i.e. last row -> x-dir #cell
        case 1
            boundary_value = 0;
        % north boundary
        case 2
            boundary_value = 0;
        % west boundary
        case 3
            boundary_value = 0; % for g = 0
        % south boundary
        case 4
            boundary_value = 0; % for g = 0
    end
    
    f = boundary_value * diag(B); % multiplying by diag(B) so the structure
    % we get is cell{#bc_ID}<--{n_eqn}<--vector(size = # nodes at boundary)
end

% function f = regular_unitstep(t) % regularized unitstep function
%     if t < 1
%         f = exp(-1/(1-(t-1)^2)+1);
%     else
%         f = 1;
%     end
% end

function[reference_order] = exact_order(x1,y1,x2,order)

y2 = y1 * exp(order * log(x2/x1));
reference_order = [x1 x2 y1 y2];
end

function f = source(x,y,ind,t)
f = -pi * sin(pi * t) * sin(pi * x).* sin(pi * y) ...
    + pi * cos(pi * t) * cos(pi * x).* sin(pi * y) ...
    + pi * cos(pi * t) * sin(pi * x).* cos(pi * y);
end



