clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


flag_ren = 0;
flag_plot_metrics = 0;



[mpc,n_res,idx_initial] = get_g('case118',flag_ren);
idx = idx_initial;
mpc = runopf(mpc,opt);

clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  simulation_hours simulation_seconds h

mpc_initial = mpc;

n_gen = size(mpc.gen,1);

n_areas = 30;

[A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss] = get_global_ss(mpc,n_areas,flag_ren);
max(A_c,[],'all')
network_initial = network;
% save('data/sim_118_30')


%%
% clear all
% clearvars -except path; close all; clc;
%load('data/sim_118_30')
h = 0.2;

simulation_hours = 0;
simulation_seconds =  500 +  3600*simulation_hours;
simulation_hours = ceil(simulation_seconds/3600);

[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);



% Controller gain synthesis 
q = zeros(1,size(A,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
q(1,freq_index) = 40;
angle_index = cumsum(bus_ss(:,2));
q(1,angle_index) = 1;

% q(1,freq_index(17)) = 10;
% q(1,freq_index(25)) = 10;
% 
% q(1,angle_index(17)) = 1;
% q(1,angle_index(25)) = 1;


%% Initial conditions

[x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:4:end,1);

 teste = Pgen0 - (PL0 + Pt0 - Ploss); 

%% Simulation Setup



t_L = 0:h:simulation_seconds;

w = zeros(size(W,2),size(t_L,2));

[w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);

P_res = x0(ren_ss,1) + w_ren;
P_load = PL0 + w_load;


P_res = P_res(:,1:3600/h:end);
P_load = P_load(:,1:3600/h:end);



%% Nonlinear Simulation


% 
% tspan = [0 simulation_seconds];
% 
% 
u_index = zeros(length(network)+1,1);
u_index(1) = 1;
for i = 1:length(network) 
    u_index(i+1) = u_index(i) + network(i).machines;
end
% 
% 
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% K = lqr(A_c,B_c,diag(q),0.1*eye(size(B,2)));
% 
% 
% [t_nL,x_nL] = ode45(@(t,x_nL) nonlinear_model(t,x_nL,K,network,bus_ss,x0,u0,P_load,P_res,Pt0,u_index), tspan,x0,opts);
% 
% 
% y_nL = C*x_nL';
% P_mech_nL = (C_mech*x_nL');


%% Nonlinear Discrete simulation

x_nL_d = zeros(size(A,1),size(t_L,2));
delta_u = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));
x_nL_d(:,1) = zeros(size(A,1),1);


K  = LQROneStepLTI(A,B,diag(q),0.1*eye(size(B,2)),E);
%K = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
x_nL_d(:,1) = x0;
k_ = 1+3600/h:3600/h:3600*simulation_hours/h + 1;
t_fault = 600/h; %seconds
mpc_initial = mpc;
%mac_remove =  [18];
mac_remove =  [];
bus_remove = [];
% for i = 1:3:n_areas
%     if(size(network(i).to_bus,1) > 1)
%         bus_remove = [bus_remove ; network(i).to_bus(1,2:3)];
%     end
% end
hour = 1;
for k = 1:length(t_L)
    k
    delta_u(:,k) = -K*(x_nL_d(:,k)-x0);
    delta_u(:,k) = min(max(delta_u(:,k),-0.1),0.1);


    if(k == t_fault)
    
        idx = idx_initial - size(mac_remove,2);
        [A,B,C,W,machine_ss,E,network,mpc,bus_ss,x_nL_d(:,k)] = fault(mpc,network,flag_ren,mac_remove,bus_remove,h,x(end,:)',bus_ss(:,2));
        [x0_,u0_,Pt0_,PL0_,~]  = initial_conditions(size(A_c,1),size(B_c,2),bus_ss(:,2),network,mpc);
        %P_load = PL0 ;
        %+ w_load;
        %P_load = P_load(:,1:3600/h:end);
    
    end 

    [~,x] = ode45(@(t,x) nonlinear_model(t,x,K,network,bus_ss,x0,u0,P_load,P_res,Pt0,u_index,delta_u(:,k)),[0 h],x_nL_d(:,k),opts);

    x_nL_d(:,k+1) = x(end,:)';
    y_nL_d(:,k) = C*(x_nL_d(:,k));




    % 
    % if(k == k_(hour) && k ~= 1)
    % 
    %         [mpc,A,B,C,D,W,~,k_tie] = update_dynamics(mpc,network,flag_ren,hour,mpc_initial,idx,h);
    % 
    %         flag = 1;
    % 
    % 
    %         hour = hour + 1;
    %         [x0,u0_,Pt0_,PL0_,~]  = initial_conditions(size(A_c,1),size(B_c,2),bus_ss(:,2),network,mpc);
    %         P_load = PL0 + w_load;
    %         P_load = P_load(:,1:3600/h:end);
    % 
    % end


end

P_mech_nL_d = (C_mech*x_nL_d);
P_mech_nL_d = P_mech_nL_d(:,1:end-1);


%% Linear simulation
% mpc = mpc_initial;
% 
% x_L = zeros(size(A,1),size(t_L,2));
% delta_u = zeros(size(B,2),size(t_L,2));
% y_L = zeros(size(C,1),size(t_L,2));
% [x0_,u0_,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B_c,2),bus_ss(:,2),network,mpc);
% x_L(:,1) = x0;
% y_L(:,1) = C*x0;
% 
% 
% 
% %K  = LQROneStepLTI(A,B,diag(q),0.1*eye(size(B,2)),E);
% %K = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
% 
% for k = 1:length(t_L)-1
% 
% 
% 
% 
%     delta_u(:,k) = -K*(x_L(:,k)-x0);
%     delta_u(:,k) = min(max(delta_u(:,k),-0.1),0.1);
% 
% 
%     x_L(:,k+1) = x0 + A*(x_L(:,k)-x0) + B*delta_u(:,k)+ W*w(:,k);
%     y_L(:,k+1) = C*(x_L(:,k+1));
% 
%     if(k == t_fault)
% 
%         idx = idx_initial - size(mac_remove,2);
%         [A,B,C,W,machine_ss,E,network,mpc,bus_ss] = fault(mpc,network,flag_ren,mac_remove,bus_remove,h,x(end,:)',bus_ss(:,2));
%         [x0_,u0_,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B_c,2),bus_ss(:,2),network,mpc);
% 
% 
%     end 
% 
%     if(k == k_(hour) && k ~= 1)
% 
%         [mpc,A,B,C,D,W,~,k_tie] = update_dynamics(mpc,network,flag_ren,hour,mpc_initial,idx,h);
% 
%         flag = 1;
% 
% 
%         hour = hour + 1;
%         [x0_,u0_,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B_c,2),bus_ss(:,2),network,mpc);
%     end
% 
% 
% end
% P_mech_L = (C_mech*(x_L+x0));

freq_limit = 0.05/50;

%%

%load('fig/Nonlinear_discrete/sims')


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(1:4:end,:)','LineWidth',1.5);
yline(1+freq_limit,'--');
yline(1-freq_limit,'--');
title('Discrete Nonlinear Simulation')
label = '$\omega$ (pu)';
ylabel(label,'interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
title_= './fig/f.fig';
savefig(title_);
set(gcf,'renderer','Painters');
title_='./fig/f.png';
saveas(gca,title_,'png');

%%
% for i = 1:n_areas
% 
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
% 
%     stairs(t_L,y_L(1+4*(i-1),:)','LineWidth',1.5);
%     %plot(t_nL,y_nL(1+4*(i-1),:),'LineWidth',1.5);
%     stairs(t_L,y_nL_d(1+4*(i-1),:)','LineWidth',1.5);
%     %xlim([0 3600])
%     %ylim([0.99995 1.00005])
%     yline(1+freq_limit,'--');
%     yline(1-freq_limit,'--');
%     %,'$\omega_{{nL}_d}$'
%     legend('$\omega_L$','$\omega_{nL}$','$\omega_{nL_d}$','Interpreter','latex','Location','best')
%     label = sprintf('$\\omega_{%d}$ (pu)',i);
%     ylabel(label,'interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%     title=sprintf('./fig/f_%d.fig',i);
%     savefig(title);
%     set(gcf,'renderer','Painters');
%     title=sprintf('./fig/f_%d.png',i);
%     saveas(gca,title,'png');
% 
% end


%%

% 
% for i = 1:n_gen
% 
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
%     stairs(t_L,P_mech_nL_d(i,:)','LineWidth',1.5);
%     %legend('$P_{{mech}_{nL}}$','$P_{{mech}_{L}}$','Interpreter','latex','Location','best')
%     label = sprintf('$P_{%d}$ (pu)',i);
%     ylabel(label,'interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% 
% end


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(2:4:end,:)','LineWidth',1.5);
legend('$ P_{m_1}$','$ P_{m_2}$','$ P_{m_3}$','Interpreter','latex')
ylabel('$ P_{m}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
hold off




figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(3:4:end,:)','LineWidth',1.5);
legend('$ P_{{tie}_1}$','$ P_{{tie}_2}$','$ P_{{tie}_3}$','Interpreter','latex')
ylabel('$ P_{tie}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(4:4:end,:)','LineWidth',1.5);
legend('$ \delta_{1}$','$ \delta_{2}$','$ \delta_{3}$','Interpreter','latex')
ylabel('$ \delta$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');




% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% plot(t_L,P_mech_L,'LineWidth',1.5);
% plot(t_nL,P_mech_nL,'LineWidth',1.5);
% legend('$ P_{m_1}$','$ P_{m_2}$','$ P_{m_3}$','$ P_{m_4}$','$ P_{m_5}$','Interpreter','latex','Location','southwest')
% ylabel('$P_{m}$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');



