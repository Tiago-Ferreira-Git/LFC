
clearvars -except path; close all; clc;
%run('../plot_options');
addpath('../../pstess/')

[g,bus,line] = get_g('IEEE_bus_118_res');
tol = 1e-8;                      % tolerance for convergence
itermax = 50;                    % maximum number of iterations
acc = 1.0;                       % acceleration factor

buses_with_generation = bitor(bus(:,10)==2 , bus(:,10)==1 );
bus(buses_with_generation,11) = 9999;
bus(buses_with_generation,12) = -9999;
bus(119:end,11) = 0.001;
bus(119:end,12) = -0.001;
bus(:,6) = 0;
bus(:,7) = 0;


bus_initial = bus;

flag_integrator = 1;
flag_ren = 1;
flag_plot_metrics = 0;


n_areas = 250;
n_areas = 30;
[A,B,C,D,W,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,n_areas,flag_ren,flag_integrator);
cond(A)
A_c = A;

%%
h = 2.5;

[A,B,W] = discrete_dynamics(A,B,h,W);

rank(ctrb(A,B))

W_ = permute_matrix(A,ren_ss);
if isempty(ren_ss)
   W_ = eye(size(A,1));
end



%number of controlable nodes
n = size(A,1);
n_C = n - size(ren_ss,2);

%plot_network(areas,line,n_areas);

%%
tic


simulation_hours = 24;
t = 0:h:3600*simulation_hours;
%7*24*3600;

w = zeros(size(W,2),size(t,2));
x0 = zeros(size(A,1),1);

mask = t > 30;

w = get_disturbance_profile(w,h,n_areas,simulation_hours,bus_ss);

R_ = 0.1;

q = zeros(1,size(A,1));
q(1) = 1;
q(1,cumsum(bus_ss(1:end-1,2))+1) = 1;


if flag_integrator
    % weight of the integrator 
    q(1,cumsum(bus_ss(:,2))) = 10;
end

% q(1,cumsum(bus_ss(1:end-1,2))) = 100;

%perfomance metrics
time_settling = zeros(1,simulation_hours);
frequency_error_cost = zeros(2,n_areas);
disp_cost_machine = zeros(2,size(B,2));
disp_cost_area = zeros(2,n_areas);
time_settling_cost = zeros(2,simulation_hours);
to_plot = zeros(simulation_hours,n_areas);

%Frequency output limit according to [1] 
freq_limit = 0.05/50;

for control_type = 1:1
    % Controller gain synthesis 
    if(control_type==1)
        decentralized = true;
    else
        decentralized = false;
    end
    if ~decentralized  E = ones(size(E)) ; end
    %E = ones(size(E));
    % if n_C == size(A,1)
        %K = get_gain(A,B,E,R_,q);
        
    %E = E*inv(W_)*[eye(n_C,n); zeros(n-n_C,n) ] *W_  ;
    %K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
    % else
    K = get_gain(A,B,E,R_,q,W_,n_C);
    % end
    if isnan(K)
        toc
        error 'Could not compute Gains'
    end

    hour = 1;
    day = 0;
    flag = 1;
    % t_settling = 0;
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

           if(24 <= hour -day*24 )
                 day = day + 1
           end
        end

        u(:,k) = -K*x(:,k);
        u(:,k) = min(max(u(:,k),-0.1),0.1);
        x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);

        y(:,k+1) = C*x(:,k+1);



        if flag_integrator
            y(1:4:end,k) = min(max(y(1:4:end,k),-freq_limit),freq_limit);
        else
            y(1:3:end,k) = min(max(y(1:3:end,k),-freq_limit),freq_limit);
        end
        %Controller performance metric
        if flag_integrator
            if all(abs(y(1:4:end-n_areas,k+1)) < 1e-9) && flag
                t_settling = ((k-k_(hour))*h );
                time_settling(1,hour) =  t_settling;
                flag = 0;

            end
        else
            if all(abs(y(1:3:end-n_areas,k+1)) < 1e-9) && flag
                t_settling = ((k-k_(hour))*h );
                time_settling(1,hour) =  t_settling;
                flag = 0;

            end
        end


        %Get the control action per area
        for i=1:n_areas
            u_area(i,k) = sum(u(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+bus_ss(i,3),3));
        end

    end

    disp_cost_area(control_type,:) = sum(u_area(1:end,:),2);
    disp_cost_machine(control_type,:) = sum(u(1:end,:),2);
    if flag_integrator
        frequency_error_cost(control_type,:) = sum(abs(y(1:4:end,:)),2)';
    else
        frequency_error_cost(control_type,:) = sum(abs(y(1:3:end,:)),2)';
    end
    time_settling_cost(control_type,:) = time_settling;
    y = y';


   
end

%%


if flag_plot_metrics
    plot_metrics(n_areas,size(B,2),simulation_hours,frequency_error_cost,disp_cost_area,disp_cost_machine,time_settling_cost,R_,to_plot)
end



figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,u','LineWidth',1.5);
legend('$u_1$','$u_2$','$u_3$','Interpreter','latex')
ylabel('$u$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_u.png';
saveas(gca,title,'png');

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;

if flag_integrator 
    stairs(t,y(:,1:4:end),'LineWidth',1.5);
    ylim([min(min(y(:,1:4:end)))*1.3,max(max(y(:,1:4:end)))*1.3])
else
    stairs(t,y(:,1:3:end),'LineWidth',1.5);
    ylim([min(min(y(:,1:3:end)))*1.3,max(max(y(:,1:3:end)))*1.3])
end
yline(freq_limit,'--');
yline(-freq_limit,'--');
legend('$\Delta\omega_1$','$\Delta\omega_2$','$\Delta\omega_3$','Interpreter','latex')
ylabel('$\Delta\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_f.png';
saveas(gca,title,'png');


%%

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
if flag_integrator 
    stairs(t,y(:,2:4:end),'LineWidth',1.5);
else
    stairs(t,y(:,2:3:end),'LineWidth',1.5);
end
legend('$\Delta P_{m_1}$','$\Delta P_{m_2}$','$\Delta P_{m_3}$','Interpreter','latex')
ylabel('$\Delta P_{m}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_m.png';
saveas(gca,title,'png');




figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
if flag_integrator 
    stairs(t,y(:,3:4:end),'LineWidth',1.5);
else
    stairs(t,y(:,3:3:end),'LineWidth',1.5);
end
legend('$\Delta P_{{tie}_1}$','$\Delta P_{{tie}_2}$','$\Delta P_{{tie}_3}$','Interpreter','latex')
ylabel('$\Delta P_{tie}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_tie.png';
saveas(gca,title,'png');

if flag_integrator
    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,y(:,4:4:end),'LineWidth',1.5);
    legend('$\Delta \delta_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{3}$','Interpreter','latex')
    ylabel('$\Delta \delta$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/delta.png';
    saveas(gca,title,'png');
end