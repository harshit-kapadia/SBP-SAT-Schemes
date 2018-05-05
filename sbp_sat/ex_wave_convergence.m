clear all
close all
clc

%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','wave-equation',... % name of example
    'n_eqn',10,... % number of eq per grid point
    'initial_condition',@initial_condition,... % it is defined below
    'theoretical_solution',@theoretical_solution,...
    'ax',[0 1 0 1],... % extents of computational domain
    'n',[100 100],... % numbers of grid cells in each coordinate direction
    't_end',0.2,... % end time of computation
    'diff_order',2,... % the difference order in the physical space
    'RK_order',2,...
    'CFL',0.5,...      % crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'output',@output... % problem-specific output routine (defined below)
    );

par.t_plot = linspace(0,par.t_end,50);

[par.system_data] = get_system_data(par.n_eqn);

% % the moment variable to be output
par.mom_output = 2;

%========================================================================
% Run solver
%========================================================================

% % solve the system
% solution = solver(par);

resolution = [16 32 64 128 256];

error_l2 = zeros(par.n_eqn,length(resolution));
for k = 1:length(resolution)                       % Loop over various grid resolutions.
    par.n = [1 1]*resolution(k);                   % Numbers of grid cells.
    solution = solver(par);                         % Run solver.
    for j = 1:1                                % Loop over solution components.
        [X,Y] = ndgrid(solution(j).x,solution(j).y);     % Grid on which solution lives.
        U_theo = theoretical_solution(X,Y,par.t_end);     % Evaluate true solution.
        error = abs(solution(j).U-U_theo);               % Difference between num. and true sol.
        int_x = dot(transpose(error),transpose(solution(j).Px * error),2);  % integral along x.
        int_xy = sum(solution(j).Py*int_x);  % integral along xy.
        error_l2(j,k) = sqrt(int_xy);%Sc. L2 error.
    end
end
figure
loglog( (resolution), error_l2(1,:), '-o' )
xlabel('# cells in each direction'), ylabel('l2-error')
title('Convergence plot')

rate = log(error_l2(1,4)/error_l2(1,1)) / log((resolution(1))/(resolution(4)));
rate

% rate = log(error_l2(1,4)/error_l2(1,1)) / log(sqrt(resolution(4))/sqrt(resolution(1)));
% rate2 = log2( error_l2(1,3)/error_l2(1,4) );

% convg_rate = zeros(par.n_eqn,length(resolution));
% for i = 1:1 % loop over components
%     for j = 2:length(resolution)
%         convg_rate(i,j) = log2( error_l2(i,j) / error_l2(i,j) );
%     end
% end


%========================================================================
% Problem Specific Functions
%========================================================================

function[system_data] = get_system_data(n_equ)

value_x = ones(n_equ);
system_data.Ax = spdiags([value_x], [0], n_equ, n_equ);

value_y = ones(n_equ);
system_data.Ay = spdiags([value_y], [0], n_equ, n_equ);

end

function f = initial_condition(x,y)
% Maxwellian/Gaussian
x0 = 0.5; % centered in the middle of domain
y0 = 0.5;
sigma_x = 0.1; % such that 6*sigma = 0.6 so in 0.2 length strip near 
sigma_y = 0.1;                           % boundary function value = 0
f = exp( -((x-x0).^2 / (2*sigma_x.^2)) - ((y-y0).^2 / (2*sigma_y.^2)) );
end

function f = theoretical_solution(x,y,t)
% Maxwellian/Gaussian
a = 1; % wave speed
x0 = 0.5 + a*t; % centered in the middle of domain
y0 = 0.5 + a*t;
sigma_x = 0.1; % such that 6*sigma = 0.6 so in 0.2 length strip near 
sigma_y = 0.1;                           % boundary function value = 0
f = exp( -((x-x0).^2 / (2*sigma_x.^2)) - ((y-y0).^2 / (2*sigma_y.^2)) );
end

function output(par,x,y,U,step)

imagesc(x,y,U'), axis xy equal tight;
title(sprintf('t = %0.2f',par.t_plot(step)));
colorbar;
xlabel('x'), ylabel('y')

drawnow
end



