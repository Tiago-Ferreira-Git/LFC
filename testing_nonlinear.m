clearvars -except path; close all; clc;
% %run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


flag_ren = 0;
flag_plot_metrics = 0;


%Column vector with buses

res_bus = [];
[mpc,n_res,idx_initial] = get_g('case14',res_bus,flag_ren);
idx = idx_initial;
mpc = runopf(mpc,opt);
clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  debug opt 

n_areas = 3;
A_c = 1;
while any(real(eig(A_c)) > 0)
    [A_c,B_c,C,~,W_c,~,~,~,E,~,network,bus_ss,ren_ss,E_fs,k_ties] = get_global_ss(mpc,n_areas,flag_ren);
    % eig(A_c)
    1;
end


%% Preload partitioned cases (ignores the code above )
% clearvars -except path; close all; clc;
% 
% load('data/sim_118_30')
% load('data/sim_57_3')

%%

h = 0.01;

simulation_hours = 0;
simulation_seconds = 200 + 3600*simulation_hours;

%%

[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);

[P,~] = permute_matrix(A,ren_ss);

% Initial conditions
[w,w_load,w_ren,P_load,P_res,u0,P_forecasted] = get_disturbance_profile(mpc,network,h,n_areas,simulation_seconds,bus_ss);




[x0,~,Pt0,~,~]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);


% Simulation Setup



t = 0:h:simulation_seconds; %time span

sampling_measurements = h*1000;
t_second = 0:sampling_measurements:simulation_seconds; % keep the results for higher sampling just to plot lighter arrays


%% Controller gain synthesis 
q = ones(1,size(A,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
angle_index = cumsum(bus_ss(:,2));


q(1,freq_index) = 4;
q(1,angle_index) =  20;


R_ = 1e3;


R = R_*eye(size(B,2));





E = E*(P\[eye(size(A,1)-n_res,size(A,1)); zeros(n_res,size(A,1))]*P);


opts_controller.verbose = true;
opts_controller.maxIT = 1e5;

[K_decentralized,~]  = LQROneStepLTI(A,B,diag(q),R,E,opts_controller);
[K_centralized,S,~]  = dlqr(A,B,diag(q),R);

%% Nonlinear Discrete simulation (centralized and decentralized)

opts_sim = odeset('RelTol',1e-8,'AbsTol',1e-8);


x_decentralized = zeros(size(A,1),size(t,2));
delta_u_decentralized= zeros(size(B,2),size(t,2));
y_decentralized= zeros(size(C,1),size(t,2));
y_d = zeros(size(C,1),size(t_second,2));
x_decentralized(:,1) = x0;

x_centralized = zeros(size(A,1),size(t,2));
delta_u_centralized = zeros(size(B,2),size(t,2));
y_centralized  = zeros(size(C,1),size(t,2));
y_c = zeros(size(C,1),size(t_second,2));
x_centralized(:,1) = x0;



tic

j = 1;

for k = 1:length(t) 


    if rem(k,1000) == 0
        k
    end


    if rem(k,3600/h) == 0
        %update controller 
        hour = k*h/3600+1;
        K_decentralized = update_dynamics(mpc,network,flag_ren,h,n_areas,simulation_seconds,hour,diag(q),R,E,0);

        K_centralized = update_dynamics(mpc,network,flag_ren,h,n_areas,simulation_seconds,hour,diag(q),R,E,1);

    end




    delta_u_decentralized(:,k) = -K_decentralized*(x_decentralized(:,k)-x0);
    delta_u_centralized(:,k) = -K_centralized*(x_centralized(:,k)-x0);



    if isempty(P_res)
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0(:,k),P_load(:,k),P_res,Pt0,delta_u_decentralized(:,k)),[0 h],x_decentralized(:,k),opts_sim);
    else
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0(:,k),P_load(:,k),P_res(:,k),Pt0,delta_u_decentralized(:,k)),[0 h],x_decentralized(:,k),opts_sim);
    end


    x_decentralized(:,k+1) = x(end,:)';
    y_decentralized(:,k) = C*(x_decentralized(:,k));





    if isempty(P_res)
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0(:,k),P_load(:,k),P_res,Pt0,delta_u_centralized(:,k)),[0 h],x_centralized(:,k),opts_sim);
    else
        [~,x] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0(:,k),P_load(:,k),P_res(:,k),Pt0,delta_u_centralized(:,k)),[0 h],x_centralized(:,k),opts_sim);
    end


    x_centralized(:,k+1) = x(end,:)';
    y_centralized(:,k) = C*(x_centralized(:,k));


    if rem(k,1/h) == 0
        y_d(:,j) = y_decentralized(:,k);
        y_c(:,j) = y_centralized(:,k);
       j = j+1;
    end



end


%%

[J_c,J_d] = plot_metrics(n_areas,size(mpc.gen,1),simulation_seconds,sampling_measurements/h,network,t,u0,delta_u_decentralized,x_decentralized,y_decentralized,delta_u_centralized,x_centralized,y_centralized,h,freq_index,angle_index,q,R,flag_plot_metrics);



%% Linear simulation 

K_cent =  dlqr(A,B,diag(q),R);

x_L = zeros(size(A,1),size(t,2));
delta_u = zeros(size(B,2),size(t,2));
y_L = zeros(size(C,1),size(t,2));
x_L(:,1) = zeros(size(A,1),1);




y_feedback = zeros(n_areas*2,size(t,2));

for k = 1:length(t)-1

    delta_u(:,k) = -K_cent*x_L(:,k);

    x_L(:,k+1) = A*x_L(:,k) + B*delta_u(:,k) + W*w(:,k);
    y_L(:,k+1) = C*(x_L(:,k+1));

end


%% Keep workspace tidy

clearvars angle_index ans bus_ss E_fs flag_ren flag_plot_metrics freq_index idx idx_initial j J_c J_d k k_ties n_res opt opts_controller opts_sim ren_ss S 