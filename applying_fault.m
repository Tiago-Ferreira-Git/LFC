
clearvars -except path; close all; clc;
run('../plot_options');
addpath('../../pstess/')

%[g,bus,line] = get_g('data3');
[g,bus,line] = get_g('IEEE_bus_118_res');

% g.line(1,:) = [];
tol = 1e-8;                      % tolerance for convergence
itermax = 50;                    % maximum number of iterations
acc = 1.0;                       % acceleration factor




[bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);

bus_initial = bus;
clearvars -except g bus t line bus_initial

flag_integrator = 1;
flag_ren = 1;
flag_plot_metrics = 1;


n_areas = 250;
n_areas = 30;
[A,B,C,D,W,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,n_areas,flag_ren);



%%
h = 2.5;

[A,B,W] = discrete_dynamics(A,B,W,h);

bus_remove = [];
for i = 1:4:size(network,2)
    if(size(network(i).to_bus,1) > 2)
        bus_remove = [bus_remove ; network(i).to_bus(1:2,2:3)];
    end
end

W_ = permute_matrix(A,ren_ss);
if isempty(ren_ss)
   W_ = eye(size(A,1));
end



%number of controlable nodes
n_C = size(A,1) - size(ren_ss,2);

%plot_network(areas,line,n_areas);

%%
tic


simulation_hours = 2;
t_fault = 1;
t_fault = 3600*t_fault/h;

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
q(1,cumsum(bus_ss(:,2))) = 100;

%perfomance metrics
time_settling = zeros(1,simulation_hours);
frequency_error_cost = zeros(2,n_areas);
disp_cost_machine = zeros(2,size(B,2));
disp_cost_area = zeros(2,n_areas);
time_settling_cost = zeros(2,simulation_hours);
to_plot = zeros(simulation_hours,n_areas);

%Frequency output limit according to [1] 
freq_limit = 0.05/50;



for control_type = 1:2
    load('ss')
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
    bus = bus_initial;
    k_ = 1:3600/h:length(t);
    x = zeros(size(A,1),size(t,2));
    u = zeros(size(B,2),size(t,2));
    u_area = zeros(n_areas,size(t,2));
    y = zeros(size(C,1),size(t,2));
    t_settling = 0;
    x(:,1) = x0;
    for k = 1:length(t)-1
        
        if(k == k_(hour+1) )
           hour = hour + 1;
           flag = 1;
        end

        u(:,k) = -K*x(:,k);
        u(:,k) = min(max(u(:,k),-0.1),0.1);
        x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);

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

        if(k == t_fault)
            [A,B,C,W,E,~,~,~] = fault(g,bus,network,flag_ren,[],bus_remove,h);
            %K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
        end
        
        k;

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
    ylim([min(y(:,1:4:end),[],'all')*1.3*1e6,max(y(:,1:4:end),[],'all')*1.3*1e6])
    yline(freq_limit*1e6,'--');
    yline(-freq_limit*1e6,'--');
    xline(t_fault*h,'LineWidth',0.5,'LineStyle','--','Color','black')
    legend('$\Delta\omega_1$','$\Delta\omega_2$','$\Delta\omega_3$','Interpreter','latex')
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


%%
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,u','LineWidth',1.5);
xline(t_fault*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$u_1$','$u_2$','$u_3$','Interpreter','latex')
ylabel('$u$ (pu)','interpreter','latex');
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
stairs(t,y(:,2:4:end),'LineWidth',1.5);
xline(t_fault*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta P_{m_1}$','$\Delta P_{m_2}$','$\Delta P_{m_3}$','Interpreter','latex')
ylabel('$\Delta P_{m}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_m.png';
saveas(gca,title,'png');




figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,y(:,3:4:end),'LineWidth',1.5);
xline(t_fault*h,'LineWidth',0.5,'LineStyle','--','Color','black')
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
xline(t_fault*h,'LineWidth',0.5,'LineStyle','--','Color','black')
legend('$\Delta \delta_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{3}$','Interpreter','latex')
ylabel('$\Delta \delta$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
xticks(0:3600:simulation_hours*3600)
xticklabels(sprintfc('%d', 0:simulation_hours))
hold off
set(gcf,'renderer','Painters');
title='./fig/delta.png';
saveas(gca,title,'png');