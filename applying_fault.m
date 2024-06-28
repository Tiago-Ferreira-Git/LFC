
clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


flag_ren = 0;
flag_plot_metrics = 0;


[mpc,n_res,idx_initial] = get_g('case14',flag_ren);
idx = idx_initial;
mpc = runopf(mpc,opt);
clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial

mpc_initial = mpc;

mac_remove =  [3 12 40];
%mac_remove = [];

n_gen = size(mpc.gen,1);

n_areas = 250;
n_areas = 30;

[A,B,C,D,W,machine_ss,C_mech,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(mpc,n_areas,flag_ren);
network_initial = network;
%%
h = 2.5;
[A,B,W] = discrete_dynamics(A,B,W,h);
save('data/ss.mat',"A","B","C","D","W")
bus_remove = [];
% for i = 1:5:n_areas
%     if(size(network(i).to_bus,1) > 1)
%         bus_remove = [bus_remove ; network(i).to_bus(1,2:3)];
%     end
% end

W_ = permute_matrix(A,ren_ss);
if isempty(ren_ss)
   W_ = eye(size(A,1));
end



%number of controlable nodes
n_C = size(A,1) - size(ren_ss,2);

%plot_network(areas,line,n_areas);

%%
tic


simulation_hours = 24;
t_fault = 13;
% t_fault = 3600*(t_fault)/h;

t = 0:h:3600*simulation_hours;
%7*24*3600;

w = zeros(size(W,2),size(t,2));
x0 = zeros(size(A,1),1);

mask = t > 30;

w = get_disturbance_profile(w,h,n_areas,simulation_hours,bus_ss);

R_ = 0.01;

q = zeros(1,size(A,1));
q(1) = 1;
q(1,cumsum(bus_ss(1:end-1,2))+1) = 1;
q(1,cumsum(bus_ss(:,2))) = 1000;

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


nominal = zeros(simulation_hours,n_gen-n_res);
nominal_fault = zeros(simulation_hours,n_gen-n_res);



mpc_ = mpc;
mpc_fault = mpc;
mask = true(1,n_gen-n_res);

mac_idx = zeros(n_gen-n_res,1);

j = 1;
for i = 1:n_areas
    mac_idx(j:j+network(i).machines-1,1) =  network(i).mac_nr;
    j = j + network(i).machines;
end

for i = 1:simulation_hours
    [mpc_,~,~,~,~,~,~,~] = update_dynamics(mpc_,network,flag_ren,i,mpc_initial,idx,h);
    nominal(i,:) = mpc_.gen(1:end-n_res,2);


    if i == t_fault
        mask = true(1,n_gen-n_res);
        mask(mac_remove) = 0;
        idx = idx_initial - size(mac_remove,2);
        [~,~,~,~,~,~,~,mpc_fault,~] = fault(mpc_fault,network,flag_ren,mac_remove,bus_remove,h);
    end
    [mpc_fault,~,~,~,~,~,~,~] = update_dynamics(mpc_fault,network,flag_ren,i,mpc_initial,idx,h);
    nominal_fault(i,mask) = mpc_fault.gen(1:end-n_res,2);
    nominal_fault(i,~mask) = 0;
end

figure
hold on
p3 = plot(1:simulation_hours,nominal./mpc.baseMVA,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.4);
p1 = plot(1:simulation_hours,nominal_fault./mpc.baseMVA,'Color',[0 0.4470 0.7410],'LineWidth',1.4);
ylabel('$P_G$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
legend([p3(1) p1(1)],{'Nominal','Nominal Fault'})
hold off
set(gcf,'renderer','Painters');
title='./fig/nominal.png';
saveas(gca,title,'png');


nominal_fault = nominal_fault(:,mac_idx);
nominal = nominal(:,mac_idx);



%Frequency output limit according to [1] 
freq_limit = 0.05/50;

w_m = ( nominal - nominal_fault)'./mpc.baseMVA;
k_ = 1:3600/h:size(w,2);
for hour=1:simulation_hours
    w_mech(:,k_(hour):k_(hour+1)) = repmat(w_m(:,hour),1,3600/h +1);
end



%%


for control_type = 1:2
    load('data/ss')
    % Controller gain synthesis 
    if(control_type==2)
        E = ones(size(E));
    end
    if n_C == size(A,1)
        %K = get_gain(A,B,E,R_,q);
        K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
    else
        K = get_gain(A,B,E,R_,q,W_,n_C);
    end
    if isnan(K)
        toc
        error 'Could not compute Gains'
    end

    hour = 1;
    day = 0;
    flag = 1;
    t_settling = 0;

    mpc = mpc_initial;
    network = network_initial;



    k_ = 1+3600/h:3600/h:length(t);
    x = zeros(size(A,1),size(t,2));
    delta_u = zeros(size(B,2),size(t,2));
    u = zeros(size(B,2),size(t,2));
    u_area = zeros(n_areas,size(t,2));
    y = zeros(size(C,1),size(t,2));
    t_settling = 0;
    x(:,1) = x0;
    for k = 1:length(t)-1
        
        if(k == k_(hour) && k ~= 1)
           
            [mpc,A,B,C,D,W,~,k_tie] = update_dynamics(mpc,network,flag_ren,hour,mpc_initial,idx,h);
            %to_plot(hour,:) = k_tie;
            
            flag = 1;
            hour = hour + 1;

            if(hour == t_fault)

                idx = idx_initial - size(mac_remove,2);
                [A,B,C,W,machine_ss,E,network,mpc,bus_ss] = fault(mpc,network,flag_ren,mac_remove,bus_remove,h);
                
            
            end 
        end
        

        
            
        
        delta_u(:,k) = -K*x(:,k);
        %delta_u(:,k) = min(max(delta_u(:,k),-0.1),0.1);

        x(:,k+1) = A*x(:,k) + B*delta_u(:,k) + B*w_m(:,hour) + W*w(:,k);
        y(:,k+1) = C*x(:,k+1);

        y(1:4:end,k) = min(max(y(1:4:end,k),-freq_limit),freq_limit);
        
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

        
        
        %delta_u(:,k)
        u(:,k) = delta_u(:,k)*100 + nominal_fault(hour,:)';

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
    stairs(t,y(:,1:4:end)*1e6,'LineWidth',1.5);
    %xlim([t_fault*h - 3600  t_fault*h + 3600])
    %ylim([min(y(:,1:4:end),[],'all')*1.3*1e6,max(y(:,1:4:end),[],'all')*1.3*1e6])
    ylim([-210 20])
    yline(freq_limit*1e6,'--');
    yline(-freq_limit*1e6,'--');
    xline((t_fault-1)*3600,'LineWidth',0.5,'LineStyle','--','Color','black')
    legend('$\Delta\omega_1$','$\Delta\omega_2$','$\Delta\omega_3$','Interpreter','latex','Location','best')
    ylabel('$\Delta\omega$ ($\mu$pu)','interpreter','latex');
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




P_mech = (C_mech*x);


%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,P_mech','LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta P_{m_1}$','$\Delta P_{m_2}$','$\Delta P_{m_3}$','$\Delta P_{m_4}$','$\Delta P_{m_5}$','Interpreter','latex','Location','southwest')
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
stairs(t,w_mech','LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$ w_{m_1}$','$w_{m_2}$','$w_{m_3}$','$w_{m_4}$','$w_{m_5}$','Interpreter','latex','Location','best')
ylabel('Disturbance (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/disturbance_mech.png';
saveas(gca,title,'png');


%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,delta_u','LineWidth',1.5);
xline((t_fault-1)*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta u_1$','$\Delta u_2$','$\Delta u_3$','Interpreter','latex','Location','northwest')
ylabel('$\Delta u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_u.png';
saveas(gca,title,'png');


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
ylabel('$\Delta \delta$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta.png';
saveas(gca,title,'png');
