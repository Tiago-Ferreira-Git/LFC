
clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


flag_ren = 0;
flag_plot_metrics = 0;


[mpc,n_res,idx_initial] = get_g('case118',flag_ren);
idx = idx_initial;
mpc = runopf(mpc,opt);
clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial

mpc_initial = mpc;

mac_remove =  [12 40];
%mac_remove =  [3];

n_gen = size(mpc.gen,1);

n_areas = 250;
n_areas = 30;

[A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss] = get_global_ss(mpc,n_areas,flag_ren);
network_initial = network;
%%
h = 0.2;
[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);
save('data/ss.mat',"A","B","C","W")
bus_remove = [];
% for i = 1:5:n_areas
%     if(size(network(i).to_bus,1) > 1)
%         bus_remove = [bus_remove ; network(i).to_bus(1,2:3)];
%     end
% end


if isempty(ren_ss)
   W_ = eye(size(A,1));
else 
   W_ = permute_matrix(A,ren_ss);
end



%number of controlable nodes
n_C = size(A,1) - size(ren_ss,2);

%plot_network(areas,line,n_areas);

%%
tic


simulation_hours = 2;
simulation_seconds = 1000 + 3600*simulation_hours;
simulation_hours = ceil(simulation_seconds/3600);
t_fault = 5;
% t_fault = 3600*(t_fault)/h;

t = 0:h:simulation_seconds;




[nominal,nominal_fault] = nominal_profiles(simulation_hours,n_gen,n_res,n_areas,mpc,network,idx_initial,t_fault,mac_remove,bus_remove,flag_ren,h,C_mech,bus_ss);

w = zeros(size(W,2),size(t,2));
[w,~,~] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);

R_ = 0.1;

q = zeros(1,size(A,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
q(1,freq_index) = 40;
angle_index = cumsum(bus_ss(:,2));
q(1,angle_index) = 1;

%perfomance metrics
time_settling = zeros(1,simulation_hours);
frequency_error_cost = zeros(2,n_areas);
disp_cost_machine = zeros(2,size(B,2));
disp_cost_area = zeros(2,n_areas);
time_settling_cost = zeros(2,simulation_hours);
to_plot = zeros(simulation_hours,n_areas);

%Active generation before/after fault
P_gen = zeros(2,size(mpc.gen,1));
Q_gen = zeros(2,size(mpc.gen,1));





%Frequency output limit according to [1] 
freq_limit = 0.05/50;




%%


for control_type = 1:1
    load('data/ss')
    % Controller gain synthesis 
    if(control_type==2)
        E = ones(size(E));
    end

    if n_C == size(A,1)
        %K = get_gain(A,B,E,R_,q);
        K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);

    else
        n = size(A_c,1);
        E_ = E*inv(W_)*[eye(n_C,n); zeros(n-n_C,n) ]*W_  ;

        if(control_type==2)
            K = dlqr(A,B,diag(q),R_*eye(size(B,2)));
        else
            K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E_);
        end
    end
    if isnan(K)
        toc
        error 'Could not compute Gains'
    end

    hour = 1;
    day = 0;
    flag = 1;

    mpc = mpc_initial;
    network = network_initial;


    [x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
    k_ = 1+3600/h:3600/h:3600*simulation_hours/h + 1;
    x = zeros(size(A,1),size(t,2));
    delta_u = zeros(size(B,2),size(t,2));
    u = zeros(size(B,2),size(t,2));
    u_area = zeros(n_areas,size(t,2));
    y = zeros(size(C,1),size(t,2));
    t_settling = 0;
    x(:,1) = x0;
    x0_ = x0;
    u0_ = u0;
    y(:,1) = C*x0;
    for k = 1:length(t)-1
        
        
        delta_u(:,k) = -K*(x(:,k)-x0) ;
        delta_u(:,k) = min(max(delta_u(:,k),-0.1),0.1);
        
        x(:,k+1) = x0 + A*(x(:,k)-x0) + B*delta_u(:,k)  + B*(u0-u0_) + W*w(:,k);
        y(:,k+1) = C*(x(:,k+1));

        y(1:4:end,k) = min(max(y(1:4:end,k),1-freq_limit),1+freq_limit);
        
        %Controller performance metric 
        if all(abs(y(1:4:end,k+1)) < 1e-9) && flag
            t_settling = ((k-k_(hour))*h );
            time_settling(1,hour) =  t_settling;
            flag = 0;

        end

        %Get the control action per area
        for i=1:n_areas
            u_area(i,k) = sum(u(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+bus_ss(i,3),3));
        end

        u(:,k) = delta_u(:,k) + u0;

        if(k == k_(hour) && k ~= 1)

            [mpc,A,B,C,D,W,~,k_tie] = update_dynamics(mpc,network,flag_ren,hour,mpc_initial,idx,h);

            flag = 1;

            if(hour == t_fault)

                idx = idx_initial - size(mac_remove,2);
                [A,B,C,W,machine_ss,E,network,mpc,bus_ss] = fault(mpc,network,flag_ren,mac_remove,bus_remove,h);


            end 
            hour = hour + 1;
            [x0_,u0_,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),n_gen,bus_ss(:,2),network,mpc);
        end


    end

    disp_cost_area(control_type,:) = sum(u_area,2);
    disp_cost_machine(control_type,:) = sum(u,2);


    frequency_error_cost(control_type,:) = sum(abs(y(1:4:end,:)),2)';
    
    time_settling_cost(control_type,:) = time_settling;
    y = y';



    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,y(:,1:4:end),'LineWidth',1.5);
    %ylim([min(y(:,1:4:end),[],'all')*1.3,max(y(:,1:4:end),[],'all')*1.3])
    line(1+ freq_limit,'--');
    yline(1 -freq_limit,'--');
    %xline((t_fault)*3600,'LineWidth',0.5,'LineStyle','--','Color','black')
    legend('$\Delta\omega_1$','$\Delta\omega_2$','$\Delta\omega_3$','Interpreter','latex','Location','best')
    ylabel('$\Delta\omega$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    xticks(0:3600:simulation_hours*3600)
    xticklabels(sprintfc('%d', 0:simulation_hours))
    hold off
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/delta_f_%d.png',control_type);
    saveas(gca,title,'png');

   
end




if flag_plot_metrics
    plot_metrics(n_areas,size(B,2),simulation_hours,frequency_error_cost,disp_cost_area,disp_cost_machine,time_settling_cost,R_,to_plot)
end



%%
P_mech = (C_mech*x);

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,P_mech','LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$ P_{m_1}$','$ P_{m_2}$','$ P_{m_3}$','$ P_{m_4}$','$ P_{m_5}$','Interpreter','latex','Location','best')
ylabel('$P_{m}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_m.png';
saveas(gca,title,'png');

%%

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,delta_u','LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta u_1$','$\Delta u_2$','$\Delta u_3$','Interpreter','latex','Location','best')
ylabel('$\Delta u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_u.png';
saveas(gca,title,'png');

%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,u','LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$u_1$','$u_2$','$u_3$','Interpreter','latex')
ylabel('$u$ (MW)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/u.png';
saveas(gca,title,'png');






figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,y(:,2:4:end),'LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta P_{m_1}$','$\Delta P_{m_2}$','$\Delta P_{m_3}$','Interpreter','latex')
ylabel('$\Delta P_{m}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_m_area.png';
saveas(gca,title,'png');




figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,y(:,3:4:end),'LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta P_{{tie}_1}$','$\Delta P_{{tie}_2}$','$\Delta P_{{tie}_3}$','Interpreter','latex')
ylabel('$\Delta P_{tie}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_tie.png';
saveas(gca,title,'png');

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,y(:,4:4:end),'LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta \delta_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{3}$','Interpreter','latex')
ylabel('$\Delta \delta$ (rads)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta.png';
saveas(gca,title,'png');
