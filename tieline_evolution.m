

clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


flag_ren = 1;
flag_plot_metrics = 0;


[mpc,n_res,idx] = get_g('case118',flag_ren);

mpc = runopf(mpc,opt);
clearvars -except mpc flag_plot_metrics flag_ren n_res idx

mpc_initial = mpc;


n_gen = size(mpc.gen,1);
n_areas = 250;
n_areas = 30;
[A,B,C,D,W,~,~,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(mpc,n_areas,flag_ren);

%%
h = 2.5;

[A,B,W] = discrete_dynamics(A,B,W,h);



%%
tic


simulation_hours = 24;
t = 0:h:3600*simulation_hours;
%7*24*3600;






hour = 1;
day = 0;
flag = 1;
k_ = 1:3600/h:length(t);





to_plot = zeros(simulation_hours,n_areas);
plot_angles = zeros(simulation_hours,size(mpc.bus,1));
plot_q_gen = zeros(simulation_hours,size(mpc.gen,1));
plot_p_gen = zeros(simulation_hours,size(mpc.gen,1));
plot_q_load = zeros(simulation_hours,size(mpc.bus,1));
plot_p_load = zeros(simulation_hours,size(mpc.bus,1));



plot_angles(1,:) = mpc.bus(:,9);
plot_q_gen(1,:) = mpc.gen(:,3);
plot_p_gen(1,:) = mpc.gen(:,2);
plot_q_load(1,:) = mpc.bus(:,4);
plot_p_load(1,:) = mpc.bus(:,3);

for k = 1:length(t)-1
    
    if(k == k_(hour+1) )
       
       hour = hour + 1;
       flag = 1;

       if(24 <= hour -day*24 )
             day = day + 1
             %K = get_gain(A,B,E,R_,q);
       end
       
        [mpc,~,~,~,~,W,~,k_tie] = update_dynamics(mpc,network,flag_ren,hour,mpc_initial,idx,h);
        
        to_plot(hour,:) = k_tie;
        plot_angles(hour,:) = mpc.bus(:,9);
        plot_q_gen(hour,:) = mpc.gen(:,3);
        plot_p_gen(hour,:) = mpc.gen(:,2);
        plot_q_load(hour,:) = mpc.bus(:,4);
        plot_p_load(hour,:) = mpc.bus(:,3);
    end


end


bus_lines = [];
for i =1:n_areas
    bus_lines = [bus_lines; network(i).to_bus(:,2)];
    bus_lines = [bus_lines; network(i).to_bus(:,3)];
end

bus_lines = unique(bus_lines);


title = sprintf('./fig/angles_%d.eps',simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,plot_angles,'LineWidth',1.4);
%p3 = plot(1:simulation_hours,plot_angles(:,bus_lines),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.4);
%p1 = plot(1:simulation_hours,plot_angles(:,ismember(1:size(mpc.bus,1),bus_lines)),'Color',[0 0.4470 0.7410],'LineWidth',1.4);
%legend([p1(1) p3(1)],{'Buses used','Buses Unused'})
ylabel('Angles [º]','interpreter','latex');
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




%%

figure
hold on
p3 = plot(1:simulation_hours,plot_p_gen(:,1:idx(1)-1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.4);
p1 = plot(1:simulation_hours,plot_p_gen(:,idx),'Color',[0 0.4470 0.7410],'LineWidth',1.4);
%p2 = plot(1:simulation_hours,data.data(1:simulation_hours,2:end)'.*100,'Color',[0.9290 0.6940 0.1250],'LineWidth',1.4);
legend([p1(1) p3(1)],{'RES','Conventional Generation','Data'})
%p2(1),'Data'
legend show
hold off
xlabel('Time - [h]');
ylabel('Active power Output (MW)');
set(gcf,'renderer','Painters');
title='./fig/generation_p.png';
saveas(gca,title,'png');


figure
hold on
p3 = plot(1:simulation_hours,plot_q_gen(:,1:idx(1)-1),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.4);
p1 = plot(1:simulation_hours,plot_q_gen(:,idx),'Color',[0 0.4470 0.7410],'LineWidth',1.4);
legend([p1(1) p3(1)],{'RES','Conventional Generation'})
hold off
xlabel('Time - [h]');
ylabel('Reactive power Output (MVar)');
set(gcf,'renderer','Painters');
title='./fig/generation_q.png';
saveas(gca,title,'png');
%%











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
