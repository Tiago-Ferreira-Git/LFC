
clearvars -except path; close all; clc;
%run('../plot_options');
addpath('../../pstess/')

[g,bus,line] = get_g('IEEE_bus_118_res');
tol = 1e-8;                      % tolerance for convergence
itermax = 50;                    % maximum number of iterations
acc = 1.0;                       % acceleration factor


bus(119:end,10) = 2;

buses_with_generation = bitor(bus(:,10)==2 , bus(:,10)==1 );
bus(buses_with_generation,11) = 9999;
bus(buses_with_generation,12) = -9999;
bus(119:end,11) = 0.001;
bus(119:end,12) = -0.001;


[bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);
bus(119:end,10) = 2;


angles = bus(:,3);
p_gen = bus(buses_with_generation,4);
q_gen = bus(buses_with_generation,5);

p_load = bus(:,6);
q_load = bus(:,7);






bus_initial = bus;

flag_integrator = 1;
flag_ren = 1;
flag_plot_metrics = 0;


n_areas = 250;
n_areas = 30;
[A,B,C,D,W,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,n_areas,flag_ren,flag_integrator);
A_c = A;

%%
h = 2.5;

[A,B,W] = discrete_dynamics(A,B,h,W);



%%
tic


simulation_hours = 24;
t = 0:h:3600*simulation_hours;
%7*24*3600;






hour = 1;
day = 0;
flag = 1;
bus = bus_initial;
k_ = 1:3600/h:length(t);





to_plot = zeros(simulation_hours,n_areas);
plot_angles = zeros(simulation_hours,size(bus,1));
plot_q_gen = zeros(simulation_hours,sum(buses_with_generation));
plot_p_gen = zeros(simulation_hours,sum(buses_with_generation));
plot_q_load = zeros(simulation_hours,size(bus,1));
plot_p_load = zeros(simulation_hours,size(bus,1));



plot_angles(1,:) = angles;
plot_q_gen(1,:) = q_gen;
plot_p_gen(1,:) = p_gen;
plot_q_load(1,:) = q_load;
plot_p_load(1,:) = p_load;

for k = 1:length(t)-1
    
    if(k == k_(hour+1) )
       
       hour = hour + 1;
       flag = 1;

       if(24 <= hour -day*24 )
             day = day + 1
             %K = get_gain(A,B,E,R_,q);
       end
       
        [bus,k_tie,angles,q_gen,p_gen,q_load,p_load] = update_dynamics(areas,bus,line,n_areas,A,bus_ss,g,network,hour,bus_initial);
        
        to_plot(hour,:) = k_tie;
        
        plot_angles(hour,:) = angles;
        plot_q_gen(hour,:) = q_gen;
        plot_p_gen(hour,:) = p_gen;
        plot_q_load(hour,:) = q_load;
        plot_p_load(hour,:) = p_load;


        to_plot(hour,:) = k_tie;
    end


end


Y = diff(plot_p_load);
mask = abs(Y(4,:)) >  0.000001;
sum(mask)

title = sprintf('./fig/angles_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,plot_angles,'LineWidth',1.5);
ylabel('$Angles$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');
title(end-2:end) = 'png';
saveas(gca,title,'png');


title = sprintf('./fig/p_load_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,plot_p_load,'LineWidth',1.5);
ylabel('$P_{load}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');
title(end-2:end) = 'png';
saveas(gca,title,'png');


title = sprintf('./fig/q_load_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,plot_q_load,'LineWidth',1.5);
ylabel('$Q_{load}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');
title(end-2:end) = 'png';
saveas(gca,title,'png');



title = sprintf('./fig/p_gen_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,plot_p_gen,'LineWidth',1.5);
ylabel('$P_{gen}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');
title(end-2:end) = 'png';
saveas(gca,title,'png');



title = sprintf('./fig/q_gen_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,plot_q_gen,'LineWidth',1.5);
ylabel('$Q_{gen}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');
title(end-2:end) = 'png';
saveas(gca,title,'png');


title = sprintf('./fig/Kties_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(2:simulation_hours,to_plot(2:end,:),'LineWidth',1.5);
ylabel('$T_{{tie}_{i}}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');
title(end-2:end) = 'png';
saveas(gca,title,'png');


covariances = zeros(n_areas,1);
means = zeros(n_areas,1);
for i =1:n_areas

    covariances(i) = cov(to_plot(2:end,i));
    means(i) = mean(to_plot(2:end,i));

    % title = sprintf('./fig/K_tie_%d_%d.png',i,simulation_hours);
    %  figure
    %  set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %  hold on
    %  grid on
    %  box on;
    %  plot(1:simulation_hours,to_plot(1:end,i)-to_plot(1,i),'LineWidth',1.5);
    %  ylabel('$T_{{tie}_{i}} - T_{{tie}_{i,0}}$ ','interpreter','latex');
    % xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    %  hold off
    %  set(gcf,'renderer','Painters');
    %  saveas(gca,title,'png');


end

sum(plot_p_gen,2) - sum(plot_p_load,2)