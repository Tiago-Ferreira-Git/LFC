clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    


debug = 2;
flag_ren = 0;
flag_plot_metrics = 0;


[mpc,n_res,idx_initial] = get_g('case14',flag_ren);
idx = idx_initial;
mpc = runopf(mpc,opt);
mpc_initial = mpc;




clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  simulation_hours simulation_seconds h debug opt mpc_initial
for i = 1:10000
    n_areas = 3;
    [A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss,E_fs] = get_global_ss(mpc,n_areas,flag_ren,debug);
    max(A_c,[],'all')
    network_initial = network;
    K = lqr(A_c,B_c,eye(size(A_c,1)),0.1*eye(size(B_c,2)));

    if any(real(eig(A_c-B_c*K)) > 0)
        1;
    end

end

%%%%%%save('data/sim_118_30')


% 
% clear all
% clearvars -except path; close all; clc;
% 
% load('data/sim_118_30')
% 
% load('data/sim_118_30')
%%

h = 0.1;

simulation_hours = 1;
simulation_seconds = 0 + 3600*simulation_hours;


[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);




% Initial conditions

[x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:3:end,1);

teste = Pgen0 - (PL0 + Pt0 ); 

% Simulation Setup



t_L = 0:h:simulation_seconds;

w = zeros(size(W,2),size(t_L,2));

[w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);


for i = 1:n_areas
    mpc.bus(network(i).bus(1),3) = mpc.bus(network(i).bus(1),3) + w_load(i,1)*100;
end

mpc = runopf(mpc,opt);

if flag_ren
    data = load('data\solar.mat');
    data = data.data;
    data.data(:,2:end) = data.data(:,2:end)./100;
    if ceil(simulation_seconds/3600) ~= 1
        P_res = resample(data.data(1:ceil(simulation_seconds/3600),2:end),3600/0.1,1,'Dimension',1);
        P_res = P_res(1:size(w_ren,2),:);
    else
        P_res = data.data(1,2:end);
    end
    % 
    P_res = P_res';
    P_res = P_res + w_ren;
    
else
    P_res = [];
end


P_load = PL0 + w_load;tspan = [0 simulation_seconds];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

K = lqr(A_c,B_c,eye(size(A_c,1)),0.1*eye(size(B_c,2)));

%Continuous linearized model 
delta_u = zeros(size(B,2),1);
%%
%[t_L,x_L] = ode45(@(t,x) linearized_model(t,x,network,bus_ss,x0,u0,w_load(:,1),P_res,delta_u,mpc.bus(:,8:9)),[0 simulation_seconds],x0,opts);
[t_L,x_L] = ode45(@(t,x) linearized_model_discrete(t,x,bus_ss,A_c,B_c,W_c,x0,u0,w_load(:,1),delta_u,K),[0 simulation_seconds],x0,opt);

y_L = C*(x_L');
if any(isnan(x_L))
    return
end

% y_L(3:3:end,:) = mod(y_L(3:3:end,:) + pi, 2*pi) - pi;
%%



% Nonlinear Simulation


freq_limit = 0.05/50;

[t_nl,x_nL] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0,P_load(:,1),P_res,delta_u,mpc.bus(:,8:9),w_load(:,1),K),t_L,x0,opts);
y_nL = C*(x_nL');

% y_nL(3:3:end,:) = mod(y_nL(3:3:end,:) + pi, 2*pi) - pi;

if any(isnan(x_nL))
    return
end

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_nl,y_nL(1:3:end,:)','LineWidth',1.5);
yline(1+freq_limit,'--');
yline(1-freq_limit,'--');
% title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
%saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
hold off



figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_L,y_L(1:3:end,:)','LineWidth',1.5);
yline(1+freq_limit,'--');
yline(1-freq_limit,'--');
% title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
%saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
hold off




%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_nl,-y_nL(3:3:end,:)','LineWidth',1.5);
% title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\theta$ (rad)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
set(gcf,'renderer','Painters');
%saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
hold off

%%
% angles = -y_nL(3:3:end,:)';
% 
% ptie = zeros(n_areas,size(t_L,1));
% 
% 
% for i = 1:size(network,2)
% 
% 
%     neighbours = network(i).to_bus;
% 
% 
%     ptie(i,:) = ptie(i,:) + neighbours(:,end)'*(angles(i)-angles(neighbours(:,1)));
% 
% 
% 
% end
% %%
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_nl,ptie','LineWidth',1.5);
% % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% ylabel('$\theta$ (rad)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
% set(gcf,'renderer','Painters');
% %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
% hold off
% 
% 
% errors = zeros(size(A,1),size(t_L,1));
% %%
% for i= 1:n_areas
% 
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
%     stairs(t_L,y_L(3*(i-1) + 1,:)','LineWidth',1.5);
%     errors(3*(i-1) + 1,:) = abs(y_nL(3*(i-1) + 1,:) -y_L(3*(i-1) + 1,:))./y_nL(3*(i-1) + 1,:);
%     stairs(t_nl,y_nL(3*(i-1) + 1,:)','LineWidth',1.5);
%     legend({'Linearized dynamics','Nonlinear dynamics'},'Interpreter','latex','Location','best')
%     ylabel('$\omega$ (pu)','interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%     %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
%     set(gcf,'renderer','Painters');
%     %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
%     hold off
% 
% 
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
%     stairs(t_L,y_L(3*(i-1) + 2,:)','LineWidth',1.5);
%     stairs(t_nl,y_nL(3*(i-1) + 2,:)','LineWidth',1.5);
%     errors(3*(i-1) + 2,:) = abs(y_nL(3*(i-1) + 2,:) -y_L(3*(i-1) + 2,:))./y_nL(3*(i-1) + 2,:);
%     legend({'Linearized dynamics','Nonlinear dynamics'},'Interpreter','latex','Location','best')
%     % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
%     ylabel('$p_m$ (pu)','interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%     %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
%     set(gcf,'renderer','Painters');
%     %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
%     hold off
% 
% 
%     figure
%     set(gca,'TickLabelInterpreter','latex') % Latex style axis
%     hold on
%     grid on
%     box on;
% 
% 
%     stairs(t_L,-y_L(3*(i-1) + 3,:)','LineWidth',1.5);
%     stairs(t_nl,-y_nL(3*(i-1) + 3,:)','LineWidth',1.5);
%     errors(3*(i-1) + 3,:) = abs(y_nL(3*(i-1) + 3,:) -y_L(3*(i-1) + 3,:))./y_nL(3*(i-1) + 3,:);
%     legend({'Linearized dynamics','Nonlinear dynamics'},'Interpreter','latex','Location','best')
%     % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
%     ylabel('$\theta$ (rad)','interpreter','latex');
%     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
%     %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
%     set(gcf,'renderer','Painters');
%     %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
%     hold off
% 
% end
% 
% 
%%
figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex')

error = (x_nL' - x_L')./(x_nL');
%error([5,10,15,20,25],:) = 0;
for i = 1:size(A,1)
    stairs(t_L,error(i,:)','LineWidth',1.5);
end
ylabel('$\epsilon_{model}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
set(gcf,'renderer','Painters');
hold off



%%





% Nonlinear Discrete simulation

x_nL_d = zeros(size(A,1),size(t_L,2));
delta_u_nld = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));
x_nL_d(:,1) = zeros(size(A,1),1);


flag_update = false;
for k = 1:length(t_L) 


    if rem(k,1000) == 0
        k
    end


    delta_u_nld(:,k) = min(max(delta_u_nld(:,k),-0.1),0.1);

    if isempty(P_res)
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0,P_load(:,k),P_res,delta_u_nld(:,k),mpc.bus(:,8:9)),[0 h],x_nL_d(:,k),opts);
    else
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0,P_load(:,k),P_res(:,k),delta_u_nld(:,k),mpc.bus(:,8:9)),[0 h],x_nL_d(:,k),opts);
    end


    x_nL_d(:,k+1) = x(end,:)';
    y_nL_d(:,k) = C*(x_nL_d(:,k));

end




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




% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,delta_u_nld','LineWidth',1.5);
% ylabel('$\Delta u$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% title(sprintf(' t_{sh} = %.2f',t_sh),'Interpreter','tex')
% savefig(sprintf('./fig/delta_u_t_{sh}_%.2f.fig',t_sh));
% set(gcf,'renderer','Painters');
% saveas(gca,sprintf('./fig/delta_u_t_{sh}_%.2f.png',t_sh),'png');
% hold off
% 
% 
% toc
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
% 
% % % freq_limit = 0.05/50;
% % % figure
% % % set(gca,'TickLabelInterpreter','latex') % Latex style axis
% % % hold on
% % % grid on
% % % box on;
% % % stairs(t_L,y_L_reduced(2:2:end,:)','LineWidth',1.5);
% % % title('Reduced Model','Interpreter','tex')
% % % ylabel('$\omega$ (pu)','interpreter','latex');
% % % xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% % % legend({'Area 1','Area 2','Area 3'},'Location','best')
% % % hold off
% 
% 
% 
% 
% %% Linear simulation - Complete Simulation
% 
% 
% x_L = zeros(size(A,1),size(t_L,2));
% delta_u = zeros(size(B,2),size(t_L,2));
% y_L = zeros(size(C,1),size(t_L,2));
% x_L(:,1) = zeros(size(A,1),1);
% %y_L(:,1) = C*x0;
% 
% % K  = LQROneStepLTI(A,B,diag(q),0.1*eye(size(B,2)),E,NaN);
% % K = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
% [K,E_fs] = slow_ss(mpc,network,h,A_c);
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
% 
% % figure
% % set(gca,'TickLabelInterpreter','latex') % Latex style axis
% % hold on
% % grid on
% % box on;
% % stairs(t_L,delta_u','LineWidth',1.5);
% hold off
