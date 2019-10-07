clear all;

nx = [60,80,100,120,140];
M = [3,4,5,6,7];
filename = cell(length(nx),1);

filename = strcat('../DVM/gaussian_collision/result_n100_DVM_20.txt');
result_dvm = dlmread(filename,'\t');
X = reshape(result_dvm(1,:),101,101);
rho_dvm = reshape(result_dvm(3,:),101,101);

for i = 1 : length(nx)
    filename = strcat('../results/gaussian_collision_kinetic/result_M',...
                         num2str(M(i)),'_n',num2str(nx(i)),'.txt');
    result{i} = dlmread(filename);
end

%% long time simulation
filename = strcat('../results/gaussian_collision_kinetic/result_long_time_M',...
                   num2str(5),'_n',num2str(80),'.txt');

result_long_time = dlmread(filename);
    
filename = strcat('../results/gaussian_collision_kinetic/result_long_time_char_M',...
                   num2str(5),'_n',num2str(80),'.txt');

result_long_time_char = dlmread(filename);

hFig = figure(1);
set(hFig, 'Position', [0 0 800 800]);

plot(result_long_time(1,:),result_long_time(2,:),'-bo',...
    result_long_time_char(1,:),result_long_time_char(2,:),'-r*',...
     'markersize',10,...
     'LineWidth',3);

xlabel('x_1');
ylabel('$\rho(x_1,x_2=0.5,t=0.6)$','Interpreter','latex');
lg = legend('kinetic penalty','characteristic penalty','Location','best'); 
lg.FontSize = 20;
grid on;
xt = get(gca, 'YTick');
set(gca, 'FontSize', 20);
output = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/long_time_t0p6';
print(output,'-depsc');

%% plots with kinetic penalty matrix
hFig = figure(2);
set(hFig, 'Position', [0 0 800 800]);
a1 = axes();
plot(result{1}(1,:),result{1}(2,:),'-bo',...
     result{2}(1,:),result{2}(2,:),'-rd',...
     result{3}(1,:),result{3}(2,:),'-g*',...
     result{4}(1,:),result{4}(2,:),'-ks',...
     result{5}(1,:),result{5}(2,:),'-c>',...
     X(:,51),rho_dvm(:,51),'--m',...
     'markersize',10,...
    'LineWidth',3);

xlabel('x_1');
ylabel('$\rho(x_1,x_2=0.5,t=0.3)$','Interpreter','latex');
legend('(M=3,N_x=60)','(M=4,N_x=80)','(M=5,N_x=100)','(M=6,N_x=120)','(M=7,N_x=140)',...
       'DVM','Location','southeast');
grid on;
xt = get(gca, 'YTick');
set(gca, 'FontSize', 20);

a2 = axes();
a2.Position = [0.15 0.6600 0.25 0.25];

plot(result{1}(1,1:5),result{1}(2,1:5),'-bo',...
     result{2}(1,1:7),result{2}(2,1:7),'-rd',...
     result{3}(1,1:8),result{3}(2,1:8),'-g*',...
     result{4}(1,1:10),result{4}(2,1:10),'-ks',...
     result{5}(1,1:12),result{5}(2,1:12),'-c>',...
     'markersize',10,...
    'LineWidth',3);

annotation('rectangle',[.135 .42 .1 .2],'LineWidth',1.5,'LineStyle','--');
annotation('arrow',[0.167 0.2],[0.62 0.66],'LineWidth',1.5);

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])


output = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/kinetic_refine';
print(output,'-depsc');

%% plots with characteristic penalty matrix

for i = 1 : length(nx)
    
    filename = strcat('../results/gaussian_collision_kinetic/result_char_M',...
                         num2str(M(i)),'_n',num2str(nx(i)),'.txt');
    
    result_char{i} = dlmread(filename);
end

hFig = figure(3);
set(hFig, 'Position', [0 0 800 800]);

plot(result_char{1}(1,:),result_char{1}(2,:),'-bo',...
     result_char{2}(1,:),result_char{2}(2,:),'-rd',...
     result_char{3}(1,:),result_char{3}(2,:),'-g*',...
     result_char{4}(1,:),result_char{4}(2,:),'-ks',...
     result_char{5}(1,:),result_char{5}(2,:),'-c>',...
     X(:,51),rho_dvm(:,51),'--m',...
     'markersize',10,...
     'LineWidth',3);

xlabel('x_1');
ylabel('$\rho(x_1,x_2=0.5,t=0.3)$','Interpreter','latex');
legend('(M=3,N_x=60)','(M=4,N_x=80)','(M=5,N_x=100)','(M=6,N_x=120)','(M=7,N_x=140)',...
       'DVM','Location','southeast');
grid on;
xt = get(gca, 'YTick');
set(gca, 'FontSize', 20);

a2 = axes();
a2.Position = [0.15 0.6600 0.25 0.25];

plot(result_char{1}(1,1:5),result_char{1}(2,1:5),'-bo',...
     result_char{2}(1,1:7),result_char{2}(2,1:7),'-rd',...
     result_char{3}(1,1:8),result_char{3}(2,1:8),'-g*',...
     result_char{4}(1,1:10),result_char{4}(2,1:10),'-ks',...
     result_char{5}(1,1:12),result_char{5}(2,1:12),'-c>',...
     'markersize',10,...
    'LineWidth',3);

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

annotation('rectangle',[.135 .42 .1 .2],'LineWidth',1.5,'LineStyle','--');
annotation('arrow',[0.167 0.2],[0.62 0.66],'LineWidth',1.5);

output = '/Users/neerajsarna/Dropbox/my_papers/Publications/Comparitive_BC/results/gaussian_collision/char_refine';
print(output,'-depsc');