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
    'n',[2 100],... % numbers of grid cells in each coordinate direction
    't_end',0.2,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',2,...
    'CFL',0.5,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'bc_inhomo',@bc_inhomo,... % source term (defined below)
    'get_penalty',@get_penalty,...
    'penalty_id',2,...
    'var_plot',1,...
    'to_plot',true,...
    'output',@output... % problem-specific output routine (defined below)
    );

par.t_plot = linspace(0,par.t_end,50);

par.c = 1 ; % wave speed
[par.system] = get_system(par.c);

[par.system] = get_penalty(par.system, par.num_bc, par.c, par.penalty_id);


%========================================================================
% Run solver and study convergence
%========================================================================

resolution = [16 32 64 128];
grid_spacing = [];

for k = 1:length(resolution)                   % Loop over various grid resolutions.
%     par.n = [1 1]*resolution(k);                   % Numbers of grid cells.
    par.n = [2 resolution(k)];                   % Numbers of grid cells.
    solution = solver(par);                         % Run solver.
    error_temp = 0;
    disp('Resolution :');
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

reference_line = exact_order(grid_spacing(1),error_L2(1),grid_spacing(end),2);

figure
loglog(grid_spacing, error_L2, '-o',reference_line(1:2),reference_line(3:4),'-*');
xlabel('h'), ylabel('l2-error');
legend('numerical','second order');
title('Convergence plot');

% convergence plot
convg_order = log(error_L2(end)/error_L2(1))/log(grid_spacing(end)/grid_spacing(1));
disp('convergence order');
disp(convg_order);

% plot error in domain --> all components separately
for j = 1:par.n_eqn
    figure
    % contourf(X,Y,error), axis xy equal tight;
    surf(X,Y,error{j}), axis xy equal tight;
    % get rid of lines in surf
    colormap summer;        ylim(par.ax([3 4]));

    shading interp;
    
    title(['Error in domain: variable ' num2str(j)]);
    colorbar;
    xlabel('x'), ylabel('y')
end


%========================================================================
% Problem Specific Functions
%========================================================================

function[system] = get_system(c)
system.Ax = [0 c^2 0; 1 0 0; 0 0 0];
system.Ay = [0 0 c^2; 0 0 0; 1 0 0];
end

function[system] = get_penalty(system, num_bc, c, penalty_id)
system.nx = [1 0 -1 0] ; % size = num_bc and east-north-west-south order
system.ny = [0 1 0 -1] ;

% we need the boundary matrix and the penalty matrix for all the
% boundaries
system.penalty_B = cell(num_bc,1);
system.penalty = cell(num_bc,1);
system.B = cell(num_bc,1);

system.A_n = cell(num_bc,1);

switch penalty_id
    case 1
        for i = 1 : num_bc
            system.A_n{i} = system.Ax*system.nx(i) + system.Ay*system.ny(i);
            system.B{i} = [1 -c^2*system.nx(i) -c^2*system.ny(i)];
            
            A_n_positive{i} = sqrt(system.A_n{i}' * system.A_n{i});
%             B1 = V * abs(Lambda) * V'
            
            [system.X_n{i}, system.Lambda{i}] = eig(system.A_n{i});
            
            diag_Lambda = diag(system.Lambda{i});
            column_index = find(diag_Lambda<0);
            system.X_n_neg{i} = system.X_n{i}(:,column_index);
            
            system.penalty{i} = 0.5 * ( (system.A_n{i}-A_n_positive{i}) *...
                system.X_n_neg{i} * inv(system.B{i}*system.X_n_neg{i}) );
            
            system.penalty_B{i} = system.penalty{i}*system.B{i};
        end
    case 2
        for i = 1 : num_bc
            system.A_n{i} = system.Ax*system.nx(i) + system.Ay*system.ny(i);
            system.B{i} = [1 -c^2*system.nx(i) -c^2*system.ny(i)];
            
            system.penalty{i} = [0; system.nx(i); system.ny(i)];
            
            system.penalty_B{i} = system.penalty{i}*system.B{i};
        end
end
system.penalty_B{1} = system.penalty_B{1} * 0;
system.penalty_B{3} = system.penalty_B{3} * 0;

system.penalty{1} = system.penalty{1} * 0;
system.penalty{3} = system.penalty{3} * 0;
end

function f = initial_condition(x,y,var_number)
t = 0;
k = 4*pi;
switch var_number
    case 1
        f = sin(k*y) .* cos(k*sqrt(2)*t);
    case 2
        f = (-1/sqrt(2)) * sin(k*y) .* sin(k*sqrt(2)*t);
    case 3
        f = (-1/sqrt(2)) * cos(k*y) .* sin(k*sqrt(2)*t);
end
end

function f = theoretical_solution(x,y,var_number,t)
k = 4*pi;
switch var_number
    case 1
        f = sin(k*y) .* cos(k*sqrt(2)*t);
    case 2
        f = (-1/sqrt(2)) * sin(k*y) .* sin(k*sqrt(2)*t);
    case 3
        f = (-1/sqrt(2)) * cos(k*y) .* sin(k*sqrt(2)*t);
end
end

function f = source(x,y,var_number,t)

% f = 0;

k = 4*pi;
switch var_number
    case 1
        f = (-1/sqrt(2)) * k * sin(k*y) .* sin(k*sqrt(2)*t);
    case 2
        f = (-1) * k * sin(k*y) .* cos(k*sqrt(2)*t);
    case 3
        f = 0;
end

end

function f = bc_inhomo(B,bc_id,U,n_eqn,t)
switch bc_id
    % east boundary but in matrix - south, i.e. last row -> x-dir #cell
    case 1
        f = 0 ;
        % north boundary
    case 2
        u_theo = [] ;
        for i = 1:n_eqn
            u_theo = [u_theo theoretical_solution(0,1,i,t)] ;
        end
        u_theo = u_theo' ;
        f = B * u_theo ; % = B * u (u is column vector with 3 variable values) (B*u should give a scalar value at the end)
        % west boundary
    case 3
        f = 0 ; % for g = 0
        % south boundary
    case 4
        u_theo = [] ;
        for i = 1:n_eqn
            u_theo = [u_theo theoretical_solution(0,0,i,t)] ;
        end
        u_theo = u_theo' ;
        f = B * u_theo ;
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

% function output(par,x,y,U,step)
% 
% % imagesc(x,y,U'), axis xy equal tight;
% 
% % surf(x,y,U'), axis xy equal tight;
% % % get rid of lines in surf
% % colormap summer;
% % shading interp;
% 
% contourf(x,y,U'), axis xy equal tight;
% 
% title(sprintf('t = %0.2f',par.t_plot(step)));
% colorbar;
% xlabel('x'), ylabel('y')
% 
% drawnow
% end



