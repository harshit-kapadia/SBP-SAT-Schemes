clear variables
close all
clc

%========================================================================
% Problem Parameters
%========================================================================
par = struct(...
    'name','advection',... % name of example
    'n_mom',3,... % order of moment approximation, % number of eq per grid point
    'source',@source,... % source term (defined below)
    'ax',[0 1 0 1],... % coordinates of computational domain
    'n',[10 10],... % numbers of grid cells in each coordinate direction
    't_end',1.0,... % the end time of the computatio
    'diff_order',2,... % the difference order in the physical space
    'CFL',0.5,...      % the crude cfl number
    'num_bc',4,... % number of boundaries in the domain
    'output',@output... % problem-specific output routine (defined below
    );

%par.t_plot = linspace(0,par.t_end,50);

par.n_eqn = (par.n_mom + 1) * (par.n_mom + 2)/2;

[par.system_data.Ax,par.system_data.Ay] = get_system_data(par.n_mom);

% the moment variable to be output
par.mom_output = 1;


% solve the system
solution = solver(par);

%========================================================================
% Problem Specific Functions
%========================================================================

function f = source(x,y)
% Radiation source (only for zeroth moment).
f = 3<x&x<4&3<y&y<4;
end

function output(par,x,y,U,step)
% Plotting routine.
cax = [-7 0];                    % Colormap range used for log10 scaling.
vneg = cax(1)-diff(cax)/254;        % Value assigned where U is negative.
V = log10(max(U,1e-50));% Cap U s.t. U>0 and use logarithmic color scale.
Vcm = max(V,cax(1));                      % Cap colormap plot from below.
Vcm(U<0) = vneg;              % Assign special value where U is negative.
clf, subplot(1,3,1:2)
imagesc(x,y,Vcm'), axis xy equal tight, caxis([vneg cax(2)])
title(sprintf('t = %0.2f',par.t_plot(step)))
cm = jet(256); cm(1,:) = [1 1 1]*.5; % Change lowest color entry to gray.
colormap(cm), colorbar('ylim',cax)
xlabel('x'), ylabel('y')
subplot(1,3,3)
plot(y,interp2(x,y,V',y*0+3.5,y))   % Evaluate solution along line x=3.5.
axis([par.ax(1:2) cax+[-1 1]*.1])
title('Cut at x=3.5'), xlabel('y')
drawnow
end


