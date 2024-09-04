clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    

% 
% 
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
% save('data/sim_118_30')
% 

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
simulation_seconds = 400 + 3600*simulation_hours;

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
    freq_value = 40;
    q(1,freq_index) = 40; 
    int_value = 5e-05;
    q(1,angle_index) = int_value;

end
% q(1,freq_index) = 40;    
% q(1,angle_index) = 1;


%% Initial conditions

[x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:3:end,1);

 teste = Pgen0 - (PL0 + Pt0 - Ploss); 

%% Simulation Setup



t_L = 0:h:simulation_seconds;

w = zeros(size(W,2),size(t_L,2));

[w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);

P_res = x0(ren_ss,1) + w_ren;
P_load = PL0 + w_load;

% 
% P_res = P_res(:,1:3600/h:end);
% P_load = P_load(:,1:3600/h:end);



%% Nonlinear Simulation




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


%% Nonlinear Discrete simulation

x_nL_d = zeros(size(A,1),size(t_L,2));
delta_u_nld = zeros(size(B,2),size(t_L,2));
y_nL_d = zeros(size(C,1),size(t_L,2));
x_nL_d(:,1) = zeros(size(A,1),1);

t_sh = 1*h;
tic


E_so = zeros(size(E));

index_ss = cumsum(bus_ss(:,2))+1;
index_ss = [1 ; index_ss];
tg_size = cumsum([ 1 ; bus_ss(:,3)]);

for i = 1:size(network,2)


    neighbours = network(i).to_bus;

    E_so(tg_size(i):tg_size(i+1)-1,index_ss(i):index_ss(i+1)-1) = ones(bus_ss(i,3),bus_ss(i,2));
    
    to_bus = [];
    unique_areas = unique(network(i).to_bus(:,1),'rows');


    for j = unique_areas'
        neighbours = [neighbours ; network(j).to_bus];
    end

    for j = 1:size(neighbours,1)
        E_so(tg_size(i):tg_size(i+1)-1,index_ss(neighbours(j,1)):index_ss(neighbours(j,1)+1)-1) = ones(bus_ss(i,3),bus_ss(neighbours(j,1),2));
    end


end



E_so;

%R_= 1000*10000;
myDir = pwd; %gets directory
myFiles = dir(fullfile(myDir,'data/first_order','*.mat')); %gets all wav files in struct

plotting = zeros(length(myFiles),4);

for j = 1:length(myFiles)
    baseFileName = myFiles(j).name;
    %baseFileName = 'K_1.000_0.10.mat';
    %baseFileName = 'K_10.000_0.10_angle_0.100000.mat'
    %baseFileName = 'K_1000000.000_0.10.mat';

    %Este é o lento
    baseFileName = 'K_1000000.000_0.10.mat';  %omega_R_1e+06_Angle_5e-05_t_{sh}_0.10_freq_40.00

   
    %baseFileName = 'K_10000000.000_0.10_q_0.00050.mat';


    %baseFileName = 'K_100000.000_0.10_angle_100.000000.mat';
    %baseFileName = 'K_10.000_0.10_angle_0.100000.mat';

    nome = split(baseFileName,'_');
    if length(nome) == 5
        freq_value = 4;    
        R_ = str2double(nome{2});
        int_value = nome{5};
        int_value = str2double(int_value(1:end-4));
    elseif length(nome) == 3
        freq_value = 40;    
        R_ = str2double(nome{2});
        int_value = 0.00005;
    else
        continue
    end




    meas = zeros(3,length(t_L));

    x_nL_d = zeros(size(A,1),size(t_L,2));
    
    delta_u_nld = zeros(size(B,2),size(t_L,2));
    y_nL_d = zeros(size(C,1),size(t_L,2));
    
    x_nL_d(:,1) = x0;
    x_nL_d_hat(:,1) = x_nL_d(:,1);
    s_h = 0.0;
    

    %load(sprintf('data/K_%.3f_%.2f.mat',R_,h));
    load(fullfile('data/first_order',baseFileName));

    %K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
    %K = dlqr(A,B,diag(q),10*eye(size(B,2)));
    
    if any(isnan(K),"all") || size(K,2) ~= 222
        continue
    end
    K_local = zeros(size(K));
    K_local(logical(E_fs)) = K(logical(E_fs));
    K_neighbour = zeros(size(K));
    K_neighbour(~logical(E_fs)) = K(~logical(E_fs));


    for k = 1:length(t_L) 
    
        if rem(k,1000) == 0
            k
        end
    
        if s_h >= t_sh || k == 1        
            dist = -K_neighbour*(x_nL_d(:,k)-x0);
            s_h = 0.0;
        end
        s_h = s_h + h;
        
    
    
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


        q(angle_index) = 0;
        q(freq_index) = freq_value;
        meas(1,k) = x_nL_d(:,k)'*diag(q)*x_nL_d(:,k);
        q(angle_index) = int_value;
        q(freq_index) = 0;
        meas(2,k) = x_nL_d(:,k)'*diag(q)*x_nL_d(:,k);
        meas(3,k) = delta_u_nld(:,k)'*R_*eye(size(B,2))*delta_u_nld(:,k);
    end

    plotting(j,1) = R_;
    plotting(j,2) = freq_value;
    plotting(j,3) = int_value;
    plotting(j,4) = max(delta_u_nld,[],'all');


    
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
    savefig(sprintf('./fig/omega_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.fig',R_,int_value,t_sh,freq_value));
    set(gcf,'renderer','Painters');
    saveas(gca,sprintf('./fig/omega_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.png',R_,int_value,t_sh,freq_value),'png');
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
    savefig(sprintf('./fig/delta_u_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.fig',R_,int_value,t_sh,freq_value));
    set(gcf,'renderer','Painters');
    saveas(gca,sprintf('./fig/delta_u_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.png',R_,int_value,t_sh,freq_value),'png');
    hold off




    figure
    hold on
    
    plot(t_L,meas(1,:))
    plot(t_L,meas(2,:))
    plot(t_L,meas(3,:))
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    ylabel('$J$ ','interpreter','latex');
    title(sprintf(' t_{sh} = %.2f',t_sh),'Interpreter','tex')
    legend({'Freq','Angle','u'})
    savefig(sprintf('./fig/meas_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.fig',R_,int_value,t_sh,freq_value));
    set(gcf,'renderer','Painters');
    saveas(gca,sprintf('./fig/meas_R_%0.1g_Angle_%0.1g_t_{sh}_%.2f_freq_%.2f.png',R_,int_value,t_sh,freq_value),'png');
    hold off

end


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
% saveas(gca,fname,'png');
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
% % saveas(gca,title,'png');
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
% % saveas(gca,title,'png');
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
mask = plotting(1:end-1,2) == 40;
to_plot_40 = plotting(mask,:);
to_plot_4 = plotting(~mask,:);
to_plot_40(:,2) = [];
to_plot_4(:,2) = [];


hold on
scatter3(to_plot_40(:,1),to_plot_40(:,2),to_plot_40(:,3))
scatter3(to_plot_4(:,1),to_plot_4(:,2),to_plot_4(:,3))
legend({"Q_{freq} = 40","Q_{freq} = 4"},"Interpreter","tex")
set(gca,'xscale','log')
set(gca,'yscale','log')
xlabel('$R $','Interpreter','latex');
ylabel('$Q_{ang} $','Interpreter','latex');
zlabel('$max(u) $','Interpreter','latex');
