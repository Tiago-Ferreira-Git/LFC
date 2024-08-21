clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    
debug = 2;


flag_ren = 0;
flag_plot_metrics = 0;


% [mpc,n_res,idx_initial] = get_g('case118',flag_ren);
% idx = idx_initial;
% mpc = runopf(mpc,opt);
% clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  simulation_hours simulation_seconds h debug opt
% 
% n_areas = 30;
% [A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss,E_fs] = get_global_ss(mpc,n_areas,flag_ren,debug);
% max(A_c,[],'all')
% network_initial = network;
%%%%%%%save('data/sim_118_30')


% clear all
% clearvars -except path; close all; clc;

load('data/sim_118_30')
% load('data/sim_14_3')
% disp('Sampling time should be:')
% pi/abs(min(real(eig(A_c))))
% disp('Sampling time should be:')
% pi/max(abs(eig(A_c)))


h = 0.1;

simulation_hours = 0;
simulation_seconds = 100 + 3600*simulation_hours;

[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);



% Controller gain synthesis 
q = zeros(1,size(A,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
angle_index = cumsum(bus_ss(:,2));
q(1,freq_index) = 40;    
q(1,angle_index) = 1;
if debug == 1

    q(1,freq_index) = 40;    
    q(1,angle_index) = 1;
else
    
    q(1,freq_index) = 40;
    q(1,angle_index) = 1;

end
% q(1,freq_index) = 40;    
% q(1,angle_index) = 0.00005;


%% Initial conditions

[x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:3:end,1);

 teste = Pgen0 - (PL0 + Pt0  ); 

%% Simulation Setup



t_L = 0:h:simulation_seconds;
            
w = zeros(size(W,2),size(t_L,2));

[w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);

P_res = x0(ren_ss,1) + w_ren;
P_load = PL0 + w_load;


P_res = P_res(:,1:3600/h:end);
P_load = P_load(:,1:3600/h:end);



%% Nonlinear Simulation



tspan = [0 simulation_seconds];


u_index = zeros(length(network)+1,1);
u_index(1) = 1;
for i = 1:length(network) 
    u_index(i+1) = u_index(i) + network(i).machines;
end
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
%K = lqr(A_c,B_c,diag(q),0.1*eye(size(B,2)));
% 
% 
% [t_nL,x_nL] = ode45(@(t,x_nL) nonlinear_model(t,x_nL,K,network,bus_ss,x0,u0,P_load,P_res,Pt0,u_index), tspan,x0,opts);
% 
% 
% y_nL = C*x_nL';
% P_mech_nL = (C_mech*x_nL');


%% Nonlinear Discrete simulation

t_sh = h*50;
t_L_obs = 0:t_sh:simulation_seconds;   

x_nL_d = zeros(size(A,1),size(t_L,2));
x_nL_d_hat = zeros(size(A,1),size(t_L_obs,2));

delta_u_nld = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));


tic

R_= 0.1;
%K = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);



%load(sprintf('data/K_%.3f_%.2f.mat',R_,h));
%R_ = 10;
K = dlqr(A,B,diag(q),R_*eye(size(B,2)));
K_local = zeros(size(K));
K_local(logical(E_fs)) = K(logical(E_fs));
K_neighbour = zeros(size(K));
K_neighbour(~logical(E_fs)) = K(~logical(E_fs));
%K = zeros(size(B,2),size(A,1));
x_nL_d(:,1) = x0;
x_nL_d_hat(:,1) = x_nL_d(:,1);
s_h = 0.0;

%%


%Observer
% 
% [A_obs,B_obs,~] = discrete_dynamics(A_c,B_c,W_c,t_sh);
% 
% G = eye(size(A_c));
% 
% Q = 0.000000000000000000000000000000001*eye(size(A_c));
% R = eye(size(C,1));
% 
% L = dlqe(A_obs,G,C,Q,R);

for k = 1:length(t_L) 

    if rem(k,1000) == 0
        k
    end

    % if k == 1
    %     dist = -K_neighbour*(x_nL_d(:,k)-x0);
    % end

    if s_h >= t_sh || k == 1

        dist = -K_neighbour*(x_nL_d(:,k)-x0);


        s_h = 0.0;
    end
    s_h = s_h + h;
    


    delta_u_nld(:,k) = -K_local*(x_nL_d(:,k)-x0)+dist;
    delta_u_nld(:,k) = min(max(delta_u_nld(:,k),-0.1),0.1);

    [~,x] = ode45(@(t,x) nonlinear_model(t,x,K,network,bus_ss,x0,u0,P_load,P_res,Pt0,u_index,delta_u_nld(:,k),debug,mpc.bus(:,8:9)),[0 h],x_nL_d(:,k),opts);

    x_nL_d(:,k+1) = x(end,:)';
    y_nL_d(:,k) = C*(x_nL_d(:,k));

    
    %x_nL_d_hat(:,k) = x_nL_d_hat(k) + L*( awgn(y_nL_d(:,k),0.000000001,'measured') - y_nL_d(:,k)  );
    %x_nL_d_hat(:,k+1) = A_obs*x_nL_d_hat(:,k) + B_obs*delta_u_nld(:,k);
    

    

end





%% Linear simulation


% x_L = zeros(size(A,1),size(t_L,2));
% delta_u = zeros(size(B,2),size(t_L,2));
% y_L = zeros(size(C,1),size(t_L,2));
% x_L(:,1) = zeros(size(A,1),1);
% y_L(:,1) = C*x0;
% 
% %K  = LQROneStepLTI(A,B,diag(q),0.1*eye(size(B,2)),E);
% %K = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
% s_h = 0.0;
% for k = 1:length(t_L)-1
% 
% 
%     %k
%     s_h = s_h + h;
%     if s_h >= t_sh || k == 1
%         dist = -K_neighbour*(x_L(:,k));
%         s_h = 0.0;
%         %k
%     end
% 
% 
%     %delta_u(:,k) = -K*x_L(:,k);
%     delta_u(:,k) = -K_local*(x_L(:,k)) + dist;
%     delta_u(:,k) = min(max(delta_u(:,k),-6),6);
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



figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(1:3:end,:)','LineWidth',1.5);
yline(1+freq_limit,'--');
yline(1-freq_limit,'--');

ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
title=sprintf('./fig/f_%f.png',t_sh);
saveas(gca,title,'png');
savefig(sprintf('./fig/f_%f.fig',t_sh));
hold off

toc
% y_nL_d_hat = C*x_nL_d_hat;
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L_obs,y_nL_d_hat(1:3:end,1:end-1)','LineWidth',1.5);
% yline(1+freq_limit,'--');
% yline(1-freq_limit,'--');
% 
% ylabel('$\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(2:3:end,:)','LineWidth',1.5);
ylabel('$P_m$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% title=sprintf('./fig/f_%f.png',t_sh);
% saveas(gca,title,'png');
hold off


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_nL_d(3:3:end,:)','LineWidth',1.5);
ylabel('$\delta$ (rad)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% title=sprintf('./fig/f_%f.png',t_sh);
% saveas(gca,title,'png');
hold off

%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,delta_u_nld','LineWidth',1.5);
ylabel('$u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% title=sprintf('./fig/f_%f.png',t_sh);
% saveas(gca,title,'png');
hold off




% 
% 
% for i = 1:n_areas
% 
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
% 
%     stairs(t_L,y_L(1+4*(i-1),:)','LineWidth',1.5);
%     plot(t_nL,y_nL(1+4*(i-1),:),'LineWidth',1.5);
%     stairs(t_L,y_nL_d(1+4*(i-1),:)','LineWidth',1.5);
%     yline(1+freq_limit,'--');
%     yline(1-freq_limit,'--');
%     legend('$\omega_L$','$\omega_{nL}$','$\omega_{nL_d}$','Interpreter','latex','Location','best')
%     label = sprintf('$\\omega_{%d}$ (pu)',i);
%     ylabel(label,'interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%     title=sprintf('./fig/f_%d.fig',i);
%     savefig(title);
%     set(gcf,'renderer','Painters');
%     title=sprintf('./fig/f_%d.png',i);
%     saveas(gca,title,'png');
%     hold off
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
%     plot(t_nL,P_mech_nL(i,:),'LineWidth',1.5);
%     stairs(t_L,P_mech_L(i,:)','LineWidth',1.5);
%     legend('$P_{{mech}_{nL}}$','$P_{{mech}_{L}}$','Interpreter','latex','Location','best')
%     label = sprintf('$P_{%d}$ (pu)',i);
%     ylabel(label,'interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% 
% end


% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% plot(t,y(2:4:end,:),'LineWidth',1.5);
% legend('$ P_{m_1}$','$ P_{m_2}$','$ P_{m_3}$','Interpreter','latex')
% ylabel('$ P_{m}$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
% hold off




% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% plot(t,y(3:4:end,:),'LineWidth',1.5);
% legend('$ P_{{tie}_1}$','$ P_{{tie}_2}$','$ P_{{tie}_3}$','Interpreter','latex')
% ylabel('$ P_{tie}$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% plot(t_nL,y_nL(4:4:end,:),'LineWidth',1.5);
% legend('$ \delta_{1}$','$ \delta_{2}$','$ \delta_{3}$','Interpreter','latex')
% ylabel('$ \delta$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');




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


%%
% for i = 1:n_gen
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
%     plot(t_L,delta_u_nld(i,:),'LineWidth',1.5);
%     %legend('$ u_{1}$','$ u_{2}$','$ u_{3}$','$ u_{4}$','$ u_{5}$','Interpreter','latex','Location','southwest')
%     ylabel('$P_{m}$ (pu)','interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% 
% end

% Keep workspace tidy

clearvars -except A A_c B B_c C C_mech E E_fs K mpc network P_load P_res Pgen0 Pl0 Ploss Pt0 h W_c L  A_obs

% G = eye(size(A_c));
% C_area = [];
% for k = 1:size(network,2)
% 
%     C_area = [C_area zeros(size(C_area,1),1) ;  zeros(1,size(C_area,2)) 1]; 
%     for i = 1:network(k).machines
%         C_area = [C_area zeros(size(C_area,1),3) ; zeros(1,size(C_area,2)) [1 1.25 0]];
%     end
%     C_area = [C_area zeros(size(C_area,1),1) ; zeros(1,size(C_area,2)) -1];
% 
% 
% end
% 
% 
% 
% Q = 0.0000000000000000000000000001*eye(size(A_c));
% R = eye(size(C_area,1));
% 
% L = dlqe(A_obs,G,C_area,Q,R);
% 
% 
% eigen_controller = eig(A-B*K);
% eigen_controller = log(eigen_controller)/h;
% 
% eigen_observer = eig(A_obs-L*C_area);
% eigen_observer = log(eigen_observer)/h;
% 
% 
% figure
% 
% 
% 
% figure
% hold on
% plot(real(eigen_controller),imag(eigen_controller),"x")
% plot(real(eigen_observer),imag(eigen_observer),"x")
% axis equal
% grid on
% xlabel("Re(z)")
% ylabel("Im(z)")
% %plot(real(eigen_cont),imag(eigen_cont),"x")
% %title(sprintf('Centralized , R = %.3f',R_))
% %legend({'Decentralized','Centralized'},'Location','best')





