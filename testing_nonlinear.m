clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    

% 
% debug = 2;
% flag_ren = 0;
% flag_plot_metrics = 0;
% 
% 
% [mpc,n_res,idx_initial] = get_g('case118',flag_ren);
% idx = idx_initial;
% mpc = runopf(mpc,opt);
% clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  simulation_hours simulation_seconds h debug opt
% 
% n_areas = 30;
% [A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss,E_fs] = get_global_ss(mpc,n_areas,flag_ren,debug);
% max(A_c,[],'all')
% network_initial = network;



%%%%%%save('data/sim_118_30')


clear all
clearvars -except path; close all; clc;

load('data/sim_118_30')
%%

h = 0.1;

simulation_hours = 0;
simulation_seconds = 400 + 3600*simulation_hours;

[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);

%



% Initial conditions

[x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:3:end,1);

teste = Pgen0 - (PL0 + Pt0 - Ploss); 

% Simulation Setup



t_L = 0:h:simulation_seconds;

w = zeros(size(W,2),size(t_L,2));

[w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);

P_res = x0(ren_ss,1) + w_ren;
P_load = PL0 + w_load;


% Nonlinear Simulation

tspan = [0 simulation_seconds];


u_index = zeros(length(network)+1,1);
u_index(1) = 1;
for i = 1:length(network) 
    u_index(i+1) = u_index(i) + network(i).machines;
end
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
% K = lqr(A_c,B_c,diag(q),0.1*eye(size(B,2)));
% 
% 
% [t_nL,x_nL] = ode45(@(t,x_nL) nonlinear_model(t,x_nL,K,network,bus_ss,x0,u0,P_load,P_res,Pt0,u_index), tspan,x0,opts);
% 
% 
% y_nL = C*x_nL';
% P_mech_nL = (C_mech*x_nL');


% Nonlinear Discrete simulation

x_nL_d = zeros(size(A,1),size(t_L,2));
delta_u_nld = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));
x_nL_d(:,1) = zeros(size(A,1),1);

t_sh = 100*h;
tic



% Controller gain synthesis 
q = zeros(1,size(A,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
angle_index = cumsum(bus_ss(:,2));

q(1,freq_index) = 10; 
q(1,angle_index) = 0.1;


meas = zeros(3,length(t_L));

x_nL_d = zeros(size(A,1),size(t_L,2));

delta_u_nld = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));

x_nL_d(:,1) = x0;
s_h = 0.0;

y_increment = C*(x_nL_d(:,1)-x0);


R_= 100;

tic
[K,~,trace_records]  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E,NaN);
toc
%[K,E_fs] = slow_ss(mpc,debug,network,h);


K_local = zeros(size(K));
K_local(logical(E_fs)) = K(logical(E_fs));
K_neighbour = zeros(size(K));
K_neighbour(~logical(E_fs)) = K(~logical(E_fs));


y_feedback = zeros(2*n_areas,1);
for k = 1:length(t_L) 
    
    % y_feedback(1:2:end) = y_increment(1:3:end);
    % y_feedback(2:2:end) = -y_increment(3:3:end);

    if rem(k,1000) == 0
        k
    end

    if s_h >= t_sh || k == 1        
        %dist = -K_neighbour*y_feedback;
        %dist = -K_neighbour*y_feedback;
        dist = -K_neighbour*(x_nL_d(:,k)-x0);
        s_h = 0.0;
    end
    s_h = s_h + h;



    %delta_u_nld(:,k) = -K_local*y_feedback+dist;

    delta_u_nld(:,k) = -K_local*(x_nL_d(:,k)-x0)+dist;
    
    %delta_u_nld(:,k) = min(max(delta_u_nld(:,k),-1),1);
    delta_u_nld(:,k) = min(max(delta_u_nld(:,k),-0.1),0.1);

    if isempty(P_res)
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,K,network,bus_ss,x0,u0,P_load(:,k),P_res,Pt0,u_index,delta_u_nld(:,k),debug),[0 h],x_nL_d(:,k),opts);
    else
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,K,network,bus_ss,x0,u0,P_load(:,k),P_res(:,k),Pt0,u_index,delta_u_nld(:,k),debug),[0 h],x_nL_d(:,k),opts);
    end

    x_nL_d(:,k+1) = x(end,:)';
    y_nL_d(:,k) = C*(x_nL_d(:,k));
    y_increment = C*(x_nL_d(:,k+1)-x0);

    % q(angle_index) = 0;
    % q(freq_index) = freq_value;
    % meas(1,k) = x_nL_d(:,k)'*diag(q)*x_nL_d(:,k);
    % q(angle_index) = int_value;
    % q(freq_index) = 0;
    % meas(2,k) = x_nL_d(:,k)'*diag(q)*x_nL_d(:,k);
    % meas(3,k) = delta_u_nld(:,k)'*R_*eye(size(B,2))*delta_u_nld(:,k);
end

% plotting(j,1) = R_;
% plotting(j,2) = freq_value;
% plotting(j,3) = int_value;
% plotting(j,4) = max(delta_u_nld,[],'all');



freq_limit = 0.05/50;
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(1:3:end,:)','LineWidth',1.5);
yline(1+freq_limit,'--');
yline(1-freq_limit,'--');
title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%savefig(sprintf('./fig/omega_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.fig',R_,int_value,t_sh,freq_value));
set(gcf,'renderer','Painters');
%saveas(gca,sprintf('./fig/omega_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.png',R_,int_value,t_sh,freq_value),'png');
hold off

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(2:3:end,:)','LineWidth',1.5);

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(3:3:end,:)','LineWidth',1.5);




figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,delta_u_nld','LineWidth',1.5);
ylabel('$\Delta u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
title(sprintf(' t_{sh} = %.2f',t_sh),'Interpreter','tex')
%savefig(sprintf('./fig/delta_u_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.fig',R_,int_value,t_sh,freq_value));
set(gcf,'renderer','Painters');
%saveas(gca,sprintf('./fig/delta_u_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.png',R_,int_value,t_sh,freq_value),'png');
hold off





toc



%% Linear simulation


% x_L = zeros(size(A,1),size(t_L,2));
% delta_u = zeros(size(B,2),size(t_L,2));
% y_L = zeros(size(C,1),size(t_L,2));
% x_L(:,1) = zeros(size(A,1),1);
% y_L(:,1) = C*x0;
% 
% %K  = LQROneStepLTI(A,B,diag(q),0.1*eye(size(B,2)),E);
% %K = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
% 
% for k = 1:length(t_L)-1
% 
% 
% 
% 
%     delta_u(:,k) = -K*x_L(:,k);
%     delta_u(:,k) = min(max(delta_u(:,k),-0.1),0.1);
% 
% 
%     x_L(:,k+1) = A*x_L(:,k) + B*delta_u(:,k)+ W*w(:,k);
%     y_L(:,k+1) = C*(x_L(:,k+1)+x0);
% 
% 
% end
% P_mech_L = (C_mech*(x_L+x0));

freq_limit = 0.05/50;

%%

%load('fig/Nonlinear_discrete/sims')

% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_nL_d(1:3:end,:)','LineWidth',1.5);
% yline(1+freq_limit,'--');
% yline(1-freq_limit,'--');
% title(sprintf('h = %.1f t_{sh} = %.1f',h,t_sh),'Interpreter','tex')
% ylabel('$\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% fname=sprintf('./fig/f_%f.png',t_sh);
% %saveas(gca,fname,'png');
% hold off

% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_nL_d(2:3:end,:)','LineWidth',1.5);
% ylabel('$P_m$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% % title=sprintf('./fig/f_%f.png',t_sh);
% % %saveas(gca,title,'png');
% hold off
% 
% %%
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_nL_d(3:3:end,:)','LineWidth',1.5);
% ylabel('$\delta$ (rad)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% % title=sprintf('./fig/f_%f.png',t_sh);
% % %saveas(gca,title,'png');
% hold off



%%
% meas = zeros(3,length(t_L));
% for k = 1:length(t_L)
%     q(angle_index) = 0;
%     q(freq_index) = 0.04;
%     meas(1,k) = x_nL_d(:,k)'*diag(q)*x_nL_d(:,k);
%     q(angle_index) = int_value;
%     q(freq_index) = 0;
%     meas(2,k) = x_nL_d(:,k)'*diag(q)*x_nL_d(:,k);
%     meas(3,k) = delta_u_nld(:,k)'*R_*eye(size(B,2))*delta_u_nld(:,k);
% end
% 
% figure
% hold on
% 
% plot(meas(1,:))
% plot(meas(2,:))
% plot(meas(3,:))
% legend({'Freq','Angle','u'})
% 

%%


figure
% mask = plotting(1:end-1,2) == 40;
% to_plot_40 = plotting(mask,:);
% to_plot_4 = plotting(~mask,:);
% to_plot_40(:,2) = [];
% to_plot_4(:,2) = [];


hold on
% scatter3(to_plot_40(:,1),to_plot_40(:,2),to_plot_40(:,3))
% scatter3(to_plot_4(:,1),to_plot_4(:,2),to_plot_4(:,3))
scatter3(plotting(:,1),plotting(:,3),plotting(:,4))
legend({"Q_{freq} = 4"},"Interpreter","tex")
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$R $','Interpreter','latex');
ylabel('$Q_{ang} $','Interpreter','latex');
zlabel('$max(u) $','Interpreter','latex');
