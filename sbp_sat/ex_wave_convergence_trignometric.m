clear all
close all
clc

%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','wave-equation',... % name of example
    'n_eqn',3,... % number of eq per grid point
    'initial_condition',@initial_condition,... % it is defined below
    'theoretical_solution',@theoretical_solution,...
    'source',@source,...
    'ax',[0 1 0 1],... % extents of computational domain
    'n',[100 4],... % numbers of grid cells in each coordinate direction
    't_end',0.2,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',2,...
    'CFL',0.5,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'get_penalty_1',@get_penalty_1,...
    'get_penalty_2',@get_penalty_2,...
    'get_boundary_operator',@get_boundary_operator,...
    'penalty_id',2,...
    'var_plot',1,...
    'to_plot',false,...
    'output',@output... % problem-specific output routine (defined below)
    );

par.t_plot = linspace(0,par.t_end,50);

par.c = 1 ; % wave speed
[par.sys] = get_sys(par.c, par.num_bc);

% [par.sys] = get_penalty(par.sys, par.num_bc, par.c, par.penalty_id);


%========================================================================
% Run solver and study convergence
%========================================================================

resolution = [16 32 64 128];
grid_spacing = [];

for k = 1:length(resolution)                   % Loop over various grid resolutions.
    %     par.n = [1 1]*resolution(k);                   % Numbers of grid cells.
    par.n = [resolution(k) 4];                   % Numbers of grid cells.
    solution = solver(par);                         % Run solver.
    error_temp = 0;
    fprintf('\n Resolution :');
    disp(par.n);
    for j = 1:par.n_eqn                               % Loop over solution components.
        X = solution(j).X;
        Y = solution(j).Y;
        PX = solution(j).PX;
        PY = solution(j).PY;
        U_theo = theoretical_solution(X,Y,j,par.t_end);     % Evaluate true solution.
        error{j} = abs(solution(j).sol-U_theo);               % Difference between num. and true sol.
        int_x = dot(transpose(error{j}),transpose(PX * error{j}),2);  % integral along x.
        int_xy = sum(PY*int_x);  % integral along xy.
        %         int_xy = sum(int_x);  % not integrating along y
        error_temp = int_xy + error_temp;%Sc. L2 error.
    end
    error_L2(k) = sqrt(error_temp);
    grid_spacing = [grid_spacing min(solution(1).h)];
    % grid_spacing = [grid_spacing max(solution(1).h)]; when cells in y or
    % x are not constant at 3 for all degree of refinement
end

% reference_line = exact_order(grid_spacing(1),error_L2(1),grid_spacing(end),2);
%
% % convergence plot
% figure
% loglog(grid_spacing, error_L2, '-o',reference_line(1:2),reference_line(3:4),'-*');
% xlabel('h'), ylabel('l2-error');
% legend('numerical','second order');
% title('Convergence plot');

convg_order = log(error_L2(end)/error_L2(1))/log(grid_spacing(end)/grid_spacing(1));
disp('convergence order');
disp(convg_order);

% % plot error in domain --> all components separately
% for j = 1:par.n_eqn
%     figure
%     % contourf(X,Y,error), axis xy equal tight;
%     surf(X,Y,error{j}), axis xy equal tight;
%     % get rid of lines in surf
%     colormap summer;        ylim(par.ax([3 4]));
%
%     shading interp;
%
%     title(['Error in domain: variable ' num2str(j)]);
%     colorbar;
%     xlabel('x'), ylabel('y')
% end


%========================================================================
% Problem Specific Functions
%========================================================================

function[sys] = get_sys(c, num_bc)
sys.Ax = [0 c^2 0; 1 0 0; 0 0 0];
sys.Ay = [0 0 c^2; 0 0 0; 1 0 0];

sys.nx = [1 0 -1 0] ; % size = num_bc and east-north-west-south order
sys.ny = [0 1 0 -1] ;

sys.A_n = cell(num_bc,1);
for i = 1 : num_bc
    sys.A_n{i} = sys.Ax*sys.nx(i) + sys.Ay*sys.ny(i);
    
    sys.A_n_positive{i} = sqrt(sys.A_n{i}' * sys.A_n{i});
    
    [sys.X_n{i}, sys.Lambda{i}] = eig(sys.A_n{i});
    diag_Lambda = diag(sys.Lambda{i});
    column_index = diag_Lambda<0;
    sys.X_n_neg{i} = sys.X_n{i}(:,column_index);
end
end

function B = get_boundary_operator(nx, ny, bc_id, y)
alpha = y.^2 ;
B = [1 -alpha*nx(bc_id) -alpha*ny(bc_id)] ;
end

function penalty = get_penalty_1(A_n, A_n_positive, X_n_neg, B, bc_id)
if bc_id == 2 || bc_id == 4
    penalty = [0; 0; 0];
else
    penalty = 0.5 * ( (A_n-A_n_positive) * X_n_neg * inv(B*X_n_neg) );
end
end

function penalty = get_penalty_2(nx, ny, bc_id)
if bc_id == 2 || bc_id == 4
    penalty = [0; 0; 0];
else
    penalty = [0; nx(bc_id); ny(bc_id)];
end
end

function f = initial_condition(x,y,var_number)
t = 0;
k = 4*pi;
switch var_number
    case 1
        f = sin(k*x) .* cos(k*sqrt(2)*t);
    case 2
        f = (-1/sqrt(2)) * cos(k*x) .* sin(k*sqrt(2)*t);
    case 3
        f = (-1/sqrt(2)) * sin(k*x) .* sin(k*sqrt(2)*t);
end
end

function f = theoretical_solution(x,y,var_number,t)
k = 4*pi;
switch var_number
    case 1
        f = sin(k*x) .* cos(k*sqrt(2)*t);
    case 2
        f = (-1/sqrt(2)) * cos(k*x) .* sin(k*sqrt(2)*t);
    case 3
        f = (-1/sqrt(2)) * sin(k*x) .* sin(k*sqrt(2)*t);
end
end

function f = source(x,y,var_number,t)

% f = 0;

k = 4*pi;
switch var_number
    case 1
        f = (-1/sqrt(2)) * k * sin(k*x) .* sin(k*sqrt(2)*t);
    case 2
        f = 0;
    case 3
        f = (-1) * k * sin(k*x) .* cos(k*sqrt(2)*t);
end

end

function f = bc_inhomo(B,bc_id,y,n_eqn,t)
switch bc_id
    % east boundary but in matrix - south, i.e. last row -> x-dir #cell
    case 1
        u_theo = [] ;
        for i = 1:n_eqn
            u_theo = [u_theo theoretical_solution(1,y,i,t)] ;
        end
        u_theo = u_theo' ;
        f = B * u_theo ; % = B * u (u is column vector with 3 variable values) (B*u should give a scalar value at the end)
        % north boundary
    case 2
        f = 0 ;
        % west boundary
    case 3
        u_theo = [] ;
        for i = 1:n_eqn
            u_theo = [u_theo theoretical_solution(0,y,i,t)] ;
        end
        u_theo = u_theo' ;
        f = B * u_theo ;
        % south boundary
    case 4
        f = 0 ; % for g = 0
end
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



