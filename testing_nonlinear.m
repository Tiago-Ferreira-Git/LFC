clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    


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

simulation_hours = 1;
simulation_seconds = 0 + 3600*simulation_hours;

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
% w_load = zeros(size(w_load));
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

q(1,freq_index) = 4;
q(1,angle_index) =  1;
R = 10000;



meas = zeros(3,length(t_L));

x_nL_d(:,1) = x0;
s_h = 0.0;

y_increment = C*(x_nL_d(:,1)-x0);

tic



K = dlqr(A,B,diag(q),R*eye(size(B,2)));


[K,E_fs] = slow_ss(mpc,debug,network,h);


% K = zeros(size(K));


K_local = zeros(size(K));
K_local(logical(E_fs)) = K(logical(E_fs));
K_neighbour = zeros(size(K));
K_neighbour(~logical(E_fs)) = K(~logical(E_fs));


y_feedback = zeros(2*n_areas,1);

flag_update = false;
for k = 1:length(t_L) 
    
    y_feedback(1:2:end) = y_increment(1:3:end);
    y_feedback(2:2:end) = -y_increment(3:3:end);

    if rem(k,1000) == 0
        k
    end


    % if rem(k,1800/h) == 0 && ~flag_update
    %     k
    %     Pmech = C_mech*x_nL_d(:,k);
    %     for i = 1:length(network)
    %         mpc.gen(network(i).mac_nr,2) = Pmech(network(i).mac_nr)*100;
    % 
    %         mask = ismember( mpc.bus(:,1),network(i).bus);
    % 
    %         mpc.bus(mask,3) = mpc.bus(mask,3) + w_load(i,k-100)./length(network(1).bus);
    % 
    %         if(y_increment(3*i) > 0)
    %             mpc.bus(mask,9) = mpc.bus(mask,9) + rem(rad2deg(y_increment(3*i)),180);
    %         else
    %             mpc.bus(mask,9) = mpc.bus(mask,9) + rem(rad2deg(y_increment(3*i)),-180);
    %         end
    % 
    %     end
    % 
    %     mpc = runopf(mpc,opt);
    % 
    %     [A_c,B_c,~,~,~,~,~,~,~,~,~,~,ren_ss,~] = get_global_ss(mpc,n_areas,flag_ren,debug,network);
    %     [x0,~,Pt0,~,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
    % 
    %     u0 = u0 + delta_u_nld(:,k);
    %     x0 = x_nL_d(:,k);
    %     x0(angle_index) = 0;
    % 
    %     PL0 = P_load(:,k-100);
    % 
    %     [A,B,~] = discrete_dynamics(A_c,B_c,W_c,h);
    % 
    %     [K,E_fs] = slow_ss(mpc,debug,network,h);
    % 
    % 
    %     K_local = zeros(size(K));
    %     K_local(logical(E_fs)) = K(logical(E_fs));
    %     K_neighbour = zeros(size(K));
    %     K_neighbour(~logical(E_fs)) = K(~logical(E_fs));
    %     flag_update = true;
    % 
    % 
    % end
    
    if s_h >= t_sh || k == 1
        if size(K,2) == size(A,1)
            dist = -K_neighbour*(x_nL_d(:,k)-x0);
        else
            dist = -K_neighbour*y_feedback;
        end
        %
        s_h = 0.0;
    end
    s_h = s_h + h;
    
    if size(K,2) == size(A,1)
        delta_u_nld(:,k) = -K_local*(x_nL_d(:,k)-x0)+dist;
    else
        delta_u_nld(:,k) = -K_local*y_feedback+dist;
    end
    
    delta_u_nld(:,k) = min(max(delta_u_nld(:,k),-0.1),0.1);
    
    if isempty(P_res)
        %nonlinear_model(t,x,K,network,bus_ss,x0,u0,PL,Pres,Pt0,u_index,delta_u,debug)
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
savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
hold off




figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,delta_u_nld','LineWidth',1.5);
ylabel('$\Delta u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
title(sprintf(' t_{sh} = %.2f',t_sh),'Interpreter','tex')
savefig(sprintf('./fig/delta_u_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
saveas(gca,sprintf('./fig/delta_u_t_{sh}_%.2f.png',t_sh),'png');
hold off


toc
%%


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(3:3:end,:)','LineWidth',1.5);
title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\Delta \delta$ (rad)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
% set(gcf,'renderer','Painters');
% saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
% hold off






