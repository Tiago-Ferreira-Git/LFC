clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    

% 
debug = 2;
flag_ren = 1;
flag_plot_metrics = 0;


%Column vector with buses
res_bus = [2];

[mpc,n_res,idx_initial] = get_g('case14',res_bus);
idx = idx_initial;
mpc = runpf(mpc,opt);
clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  simulation_hours simulation_seconds h debug opt

n_areas = 5;
A_c = 1;
while any(real(eig(A_c)) > 0)
    [A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss,E_fs,k_ties] = get_global_ss(mpc,n_areas,flag_ren,debug);
end


max(A_c,[],'all')
network_initial = network;


    




%save('data/sim_14_3')


% 
% clear all
% clearvars -except path; close all; clc;
% 
% load('data/sim_118_30')

figure
plot(real(eig(A_c)),imag(eig(A_c)),'x')





h = 0.1;

simulation_hours = 0;
simulation_seconds = 1000 + 3600*simulation_hours;



[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);

[P,~] = permute_matrix(A,ren_ss);

% Initial conditions

[x0,~,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:3:end,1);

teste = Pgen0 - (PL0 + Pt0 ); 

% Simulation Setup



t_L = 0:h:simulation_seconds;

%%
%(mpc,network,h,n_areas,simulation_seconds,PL0,n_res,res_buses)
[w,w_load,w_ren,P_load,P_res,u0] = get_disturbance_profile(mpc,network,h,n_areas,simulation_seconds,bus_ss,PL0);




% Nonlinear Simulation

tspan = [0 simulation_seconds];



opts = odeset('RelTol',1e-8,'AbsTol',1e-8);



% Nonlinear Discrete simulation

x_nL_d = zeros(size(A,1),size(t_L,2));
delta_u_nld = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));
x_nL_d(:,1) = zeros(size(A,1),1);

t_sh = 1*h;
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
q(freq_index) = 10.*k_ties;
q(angle_index) = 0.01.*k_ties;


R_ties = zeros(1,size(B,2));
j = 1;
for i = 1:size(network,2)
    R_ties(j:j-1+network(i).machines) = k_ties(i);
    j = j + network(i).machines;
end

% abs(eig(A-B*K*C_))

% if ~isempty(ren_ss)
%     A_hat = P\A*P;
%     A_hat = A_hat(1:end-n_res,1:end-n_res);
%     B_hat = P\B;
% 
%     K  = LQROneStepLTI(A,B,diag(q),diag(1.*R_ties),E,NaN);
%     j;
% 
% else
%     K  = LQROneStepLTI(A,B,diag(q),diag(1.*R_ties),E,NaN);
% end

E = E*(P\[eye(size(A,1)-n_res,size(A,1)); zeros(n_res,size(A,1))]*P);


 K  = LQROneStepLTI(A,B,diag(q),diag(1000.*R_ties),E,NaN);

%%
%K = dlqr(A,B,diag(q),R*eye(size(B,2)));


% [K,E_fs] = slow_ss(mpc,network,h,A_c);
% 
% 
% K = zeros(size(K));
% 
% 
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

    % 
    % if s_h >= t_sh || k == 1
    %     if size(K,2) == size(A,1)
    %         dist = -K_neighbour*(x_nL_d(:,k)-x0);
    %     else
    %         dist = -K_neighbour*y_feedback;
    %     end
    %     %
    %     s_h = 0.0;
    % end
    % s_h = s_h + h;
    % 
    % if size(K,2) == size(A,1)
    %     delta_u_nld(:,k) = -K_local*(x_nL_d(:,k)-x0)+dist;
    % else
    %     delta_u_nld(:,k) = -K_local*y_feedback+dist;
    % end

    % if rem(k,3600/h) == 0
    %     %update controller 
    %     hour = k*h/3600;
    %     Pmech = C_mech*x_nL_d(:,k);
    %     [mpc,A,B,C,D,~,machine_ss,~] = update_dynamics(mpc,network,flag_ren,hour,h,Pmech);
    % 
    %     [x0,~,Pt0,~,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
    % 
    %     u0 = u0 + delta_u_nld(:,k);
    %     x0 = x_nL_d(:,k);
    %     x0(angle_index) = 0;
    % 
    %     % PL0 = P_load(:,k-100);
    % 
    % 
    %     k
    % end


    % 

    if size(K,2) == size(A,1)
        delta_u_nld(:,k) = -K*(x_nL_d(:,k)-x0);
    else
         delta_u_nld(:,k) = -K*y_feedback;
    end

    %%u = max(u,0);

    %delta_u_nld(:,k) = min(max(delta_u_nld(:,k),-0.1),0.1);

    if isempty(P_res)
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0(:,k),P_load(:,k),P_res,Pt0,delta_u_nld(:,k),debug),[0 h],x_nL_d(:,k),opts);
    else
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0(:,k),P_load(:,k),P_res(:,k),Pt0,delta_u_nld(:,k),debug),[0 h],x_nL_d(:,k),opts);
    end


    x_nL_d(:,k+1) = x(end,:)';
    y_nL_d(:,k) = C*(x_nL_d(:,k));
    y_increment = y_nL_d(:,k) - y_nL_d(:,1);

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
stairs(t_L,y_nL_d(1:3:end,:)'-1,'LineWidth',1.5);
% yline(freq_limit,'--');
% yline(freq_limit,'--');
title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
hold off


%%

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



figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,(max(delta_u_nld+u0,0))','LineWidth',1.5);
ylabel('$ u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
title(sprintf(' t_{sh} = %.2f',t_sh),'Interpreter','tex')
savefig(sprintf('./fig/u_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
saveas(gca,sprintf('./fig/u_t_{sh}_%.2f.png',t_sh),'png');
hold off



toc
% %%
% 
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_nL_d(3:3:end,:)','LineWidth',1.5);
% title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% ylabel('$\Delta \delta$ (rad)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% % savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
% % set(gcf,'renderer','Painters');
% % saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
% % hold off




%% Linear simulation - Reduced Model Simulation
% [K,E_fs,A_reduced,B_reduced,W_reduced] = slow_ss(mpc,network,h,A_c);
% 
% C_reduced = eye(size(A_reduced));
% 
% x_L_reduced = zeros(size(A_reduced,1),size(t_L,2));
% delta_u_reduced  = zeros(size(B_reduced,2),size(t_L,2));
% y_L_reduced  = zeros(size(C_reduced,1),size(t_L,2));
% x_L_reduced(:,1) = zeros(size(A_reduced,1),1);
% 
% 
% for k = 1:length(t_L)-1
% 
% 
%     delta_u_reduced(:,k) = -K*x_L_reduced(:,k);
%     delta_u_reduced(:,k) = min(max(delta_u_reduced(:,k),-0.1),0.1);
% 
%     x_L_reduced(:,k+1) = A_reduced*x_L_reduced(:,k) + B_reduced*delta_u_reduced(:,k)+ W_reduced*w(:,k);
%     y_L_reduced(:,k+1) = C_reduced*(x_L_reduced(:,k+1));
% 
% end
% 
% 
% freq_limit = 0.05/50;
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_L_reduced(1:2:end,:)','LineWidth',1.5);
% title('Reduced Linearized Model','Interpreter','tex')
% ylabel('$\Delta \omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% legend({'Area 1','Area 2','Area 3'},'Location','best')
% hold off
% 
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,delta_u_reduced','LineWidth',1.5);
% ylabel('$\Delta u$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% 
% freq_limit = 0.05/50;
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_L_reduced(2:2:end,:)','LineWidth',1.5);
% title('Reduced Model','Interpreter','tex')
% ylabel('$\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% legend({'Area 1','Area 2','Area 3'},'Location','best')
% hold off
% mask = abs(delta_u_reduced(:,end)) < 0.1; 
% if sum(mask) < 10
%     1;
% end




%% Linear simulation - Complete Simulation


% x_L = zeros(size(A,1),size(t_L,2));
% delta_u = zeros(size(B,2),size(t_L,2));
% y_L = zeros(size(C,1),size(t_L,2));
% x_L(:,1) = zeros(size(A,1),1);
% %y_L(:,1) = C*x0;
% 
% 
% %Closed loop in the original system
% C_ = zeros(2*length(network),size(C,2));
% C_(1:2:end) = C(1:3:end,:);
% C_(2:2:end) = -C(3:3:end,:);
% 
% 
% %K = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
% %[K,E_fs] = slow_ss(mpc,network,h,A_c);
% 
% 
% y_feedback = zeros(2*n_areas,size(t_L,2));
% 
% for k = 1:length(t_L)-1
% 
%     y_feedback(1:2:end,k) = y_L(1:3:end,k);
%     y_feedback(2:2:end,k) = -y_L(3:3:end,k);
% 
% 
%     if size(K,2) == size(A,1)
%         delta_u(:,k) = -K*x_L(:,k);
%     else
%          delta_u(:,k) = -K*y_feedback(:,k);
%     end
% 
% 
%     delta_u(:,k) = min(max(delta_u(:,k),-0.1),0.1);
% 
% 
% 
%     x_L(:,k+1) = A*x_L(:,k) + B*delta_u(:,k)+ W*w(:,k);
%     y_L(:,k+1) = C*(x_L(:,k+1));
% 
%     %y_increment = C*x_L(:,k);
% 
% 
% end
% P_mech_L = (C_mech*(x_L+x0));
% 
% 
% 
% 
% freq_limit = 0.05/50;
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_L(1:3:end,:)','LineWidth',1.5);
% title('Full Linearized Model','Interpreter','tex')
% ylabel('$\Delta\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% legend({'Area 1','Area 2','Area 3'},'Location','best')
% hold off
% %%
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,delta_u','LineWidth',1.5);
% ylabel('$\Delta u$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% hold off
