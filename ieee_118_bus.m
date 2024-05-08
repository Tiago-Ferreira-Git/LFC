
clearvars -except path; close all; clc;
run('plot_options');

% %% Run nonlinear simulation and store results
% cd('../../pstess/');
% %set_path('IEEE_bus_14');
% set_path('IEEE_bus_118.m');  % set the working directory in Sget_path
% %no_res_matpower
% 
% run('s_simu');  % run PSTess
% cd('../analysis/area partition/');


[g,bus,line] = get_g('IEEE_bus_118_res');

tol = 1e-8;                      % tolerance for convergence
itermax = 50;                    % maximum number of iterations
acc = 1.0;                       % acceleration factor

[bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);

bus_initial = bus;
clearvars -except g bus t line bus_initial

%%
n_areas = 30;
%n_areas = 3;
[A_global,B_global,C_global,D_global,W_global,~,E_global,L,areas,network,bus_ss] = get_global_ss(g,bus,n_areas,0);



h = 2.5;

[A,B,C,W,~,E] = discrete_dynamics(A_global,B_global,C_global,D_global,W_global,E_global,bus_ss,h,n_areas,L);
C_angle = zeros(n_areas,size(C,2));

C_angle(1:n_areas,end-n_areas+1:end) = -eye(n_areas);
C = [C ; C_angle];
%ones(size(E_global,1),n_areas);


%%
tic


simulation_hours = 168;
t = 0:h:3600*simulation_hours;
%7*24*3600;

w = zeros(size(W_global,2),size(t,2));
x0 = zeros(size(A,1),1);

mask = t > 30;

w = get_disturbance_profile(w,h,n_areas);

R_ = 0.01;

q = zeros(1,size(A,1));
q(1,1) = 1;
for i=1:n_areas
    
    q(1,sum(bus_ss(1:i,2))+1) = 1;
end

% weight of the integrator 


q(1,size(A_global,1)+1:end) = 10;



%perfomance metrics

time_settling = zeros(1,simulation_hours);

frequency_error_cost = zeros(2,n_areas);
disp_cost_machine = zeros(2,size(B,2));
disp_cost_area = zeros(2,n_areas);
time_settling_cost = zeros(2,simulation_hours);

to_plot = zeros(simulation_hours,n_areas);


for control_type = 1:2
    % Controller gain synthesis 
    if(control_type==1)
        decentralized = true;
    else
        decentralized = false;
    end
    if ~decentralized  E = ones(size(E)) ; end
    K = get_gain(A,B,E,R_,q);
    if isnan(K)
        toc
        error 'Could not compute Gains'
    end
    
    hour = 1;
    day = 0;
    flag = 1;
    % t_settling = 0;
    bus = bus_initial;
    k_ = 0;
    x = zeros(size(A,1),size(t,2));
    u = zeros(size(B,2),size(t,2));
    u_area = zeros(n_areas,size(t,2));
    y = zeros(size(C,1),size(t,2));
    x(:,1) = x0;
    for k = 1:length(t)-1
        
        if(3600 < (k-k_)*h )
           hour = hour + 1;
           flag = 1;
           k_ = k;

           if(24 <= hour -day*24 )
                 day = day + 1
                 hour-day*24+1;
                 %K = get_gain(A,B,E,R_,q);
           end
           [A_global,~,k_tie] = update_dynamics(areas,bus,bus_initial,line,n_areas,A_global,bus_ss,g,network,hour-day*24+1);
           [A,B,~,W,~,~] = discrete_dynamics(A_global,B_global,C_global,D_global,W_global,E_global,bus_ss,h,n_areas,L);
            

           
            to_plot(hour,:) = k_tie;
        end

        u(:,k) = -K*x(:,k);
        u(:,k) = min(max(u(:,k),-0.1),0.1);
        x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);

        y(:,k+1) = C*x(:,k+1);



        y(1:3:end-n_areas,k) = min(max(y(1:3:end-n_areas,k),-0.04*2*pi),0.04*2*pi);
        %Frequency output limit according to [1] 
        freq_limit = 0.05/50;
        y(1:3:end-n_areas,k) = min(max(y(1:3:end-n_areas,k),-freq_limit),freq_limit);

        %Controller performance metric
        t_settling = ((k-k_)*h );
        if all(abs(y(1:3:end-n_areas,k+1)) < 1e-9) && flag
            time_settling(1,hour) =  t_settling;
            flag = 0;

        end


        %Get the control action per area
        for i=1:n_areas
            u_area(i,k) = sum(u(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+bus_ss(i,3),3));
        end

    end

    disp_cost_area(control_type,:) = sum(u_area(1:end,:),2);
    disp_cost_machine(control_type,:) = sum(u(1:end,:),2);
    frequency_error_cost(control_type,:) = sum(abs(y(1:3:end-n_areas,:)),2)';
    time_settling_cost(control_type,:) = time_settling;
    y = y';

    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,y(:,1:3:end-n_areas),'LineWidth',1.5);
    yline(freq_limit,'--');
    yline(-freq_limit,'--');
    ylim([min(min(y(:,1:3:end-n_areas)))*1.3,max(max(y(:,1:3:end-n_areas)))*1.3])
    legend('$\Delta\omega_1$','$\Delta\omega_2$','$\Delta\omega_3$','Interpreter','latex')
    ylabel('$\Delta\omega$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/delta_f.png';
    saveas(gca,title,'png');
   
end


figure; 
hold on;
grid on;
box on;
set(gca,'TickLabelInterpreter','latex') % Latex style axis
%%%%
plot(1:n_areas,frequency_error_cost(1,:),1:n_areas,frequency_error_cost(2,:),'LineWidth',1.5)
%%%%
legend({'Decentralized','Centralized'},...
	'Location','best','Interpreter','latex');
ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
xlabel('$Area $','Interpreter','latex');
hold off;
set(gcf,'renderer','Painters');
title=sprintf('./fig/error_R%f.png',R_);
saveas(gca,title,'png');


figure
hold on;
grid on;
box on;
set(gca,'TickLabelInterpreter','latex') % Latex style axis
%%%%
plot(1:n_areas,disp_cost_area(1,:),1:n_areas,disp_cost_area(2,:),'LineWidth',1.5)
%%%%
legend({'Decentralized','Centralized'},...
	'Location','best','Interpreter','latex');
ylabel('$\sum_{k=1}^{k_{end}} \Delta u_i(k) (\mathrm{pu})$','Interpreter','latex');
xlabel('$Area $','Interpreter','latex');
hold off;
set(gcf,'renderer','Painters');
title=sprintf('./fig/u_area_R%f.png',R_);
saveas(gca,title,'png');




figure
hold on;
grid on;
box on;
set(gca,'TickLabelInterpreter','latex') % Latex style axis
%%%%
plot(1:size(B,2),disp_cost_machine(1,:),1:size(B,2),disp_cost_machine(2,:),'LineWidth',1.5)
%%%%
legend({'Decentralized','Centralized'},...
	'Location','best','Interpreter','latex');
ylabel('$\sum_{k=1}^{k_{end}} \Delta u_i(k) (\mathrm{pu})$','Interpreter','latex');
xlabel('$Machine $','Interpreter','latex');
hold off;
set(gcf,'renderer','Painters');
title=sprintf('./fig/u_machine_R%f.png',R_);
saveas(gca,title,'png');




figure; 
hold on;
grid on;
box on;
set(gca,'TickLabelInterpreter','latex') % Latex style axis
%%%%
plot(1:simulation_hours,time_settling_cost(1,:),1:simulation_hours,time_settling_cost(2,:),'LineWidth',1.5)
%%%%
legend({'Decentralized','Centralized'},...
	'Location','best','Interpreter','latex');
ylabel('Settling time $(\mathrm{s})$','Interpreter','latex');
xlabel('Hours $(\mathrm{h})$','Interpreter','latex');
hold off;
set(gcf,'renderer','Painters');
title=sprintf('./fig/time_settling_R%f.png',R_);
saveas(gca,title,'png');






%%
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
%%

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,y(:,2:3:end-n_areas),'LineWidth',1.5);
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
stairs(t,y(:,3:3:end-n_areas),'LineWidth',1.5);
legend('$\Delta P_{{tie}_1}$','$\Delta P_{{tie}_2}$','$\Delta P_{{tie}_3}$','Interpreter','latex')
ylabel('$\Delta P_{tie}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_tie.png';
saveas(gca,title,'png');


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,y(:,end-n_areas:end),'LineWidth',1.5);
legend('$\Delta \delta_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{3}$','Interpreter','latex')
ylabel('$\Delta \delta$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta.png';
saveas(gca,title,'png');


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,w(1:2:end,:)','LineWidth',1.5);
legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
ylabel('$\Delta P_{load}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/delta_p_load.png';
saveas(gca,title,'png');



figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t,w(2:2:end,:)','LineWidth',1.5);
legend('$\Delta P_{{ren}_1}$','$\Delta P_{{ren}_2}$','$\Delta P_{{ren}_3}$','Interpreter','latex')
ylabel('$\Delta P_{ren}$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/p_ren.png';
saveas(gca,title,'png');



covariances = zeros(n_areas,1);
means = zeros(n_areas,1);
for i =1:n_areas

    covariances(i) = cov(to_plot(2:end-3,i));
    means(i) = mean(to_plot(2:end-3,i));

    % title = sprintf('./fig/K_tie_%d_%d.png',i,simulation_hours);
    % figure
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % hold on
    % grid on
    % box on;
    % plot(1:simulation_hours-5,to_plot(2:end-4,i)-to_plot(2,i),'LineWidth',1.5);
    % ylabel('$T_{{tie}_{i}} - T_{{tie}_{i,0}}$ ','interpreter','latex');
    % xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    % hold off
    % set(gcf,'renderer','Painters');
    % saveas(gca,title,'png');


end

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(sqrt(covariances),'LineWidth',1.5);
%legend('$\sigma_1$','$\sigma_2}$','$\sigma_3}$','Interpreter','latex')
%y_label = sprintf('$T_{{tie}_{%d,k}} - T_{{tie}_{%d,0}}$',i,i);
ylabel('$\sigma$','interpreter','latex');
xlabel('$area$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = sprintf('./fig/K_tie_%d_%d.png',10000,simulation_hours);
saveas(gca,title,'png');


figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(means,'LineWidth',1.5);
%legend('$\sigma_1$','$\sigma_2}$','$\sigma_3}$','Interpreter','latex')
%y_label = sprintf('$T_{{tie}_{%d,k}} - T_{{tie}_{%d,0}}$',i,i);
ylabel('$\mu$','interpreter','latex');
xlabel('$area$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = sprintf('./fig/K_tie_%d_%d.png',10123,simulation_hours);
saveas(gca,title,'png');

toc


% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% plot(1:simulation_hours-1,to_plot(2:end,:),'LineWidth',1.5);
% %legend('$\Delta T_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{3}$','Interpreter','latex')
% ylabel('$ T_{{tie}_{i}}$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% hold off
% set(gcf,'renderer','Painters');
% title='./fig/K_tie.png';
% saveas(gca,title,'png');



figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;

c_ = zeros(n_areas,210);
c_(1:n_areas,end-n_areas+1:end) = -eye(n_areas);

stairs(t,(c_*x)','LineWidth',1.5);
legend('$\Delta \delta_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{2}$','Interpreter','latex')
ylabel('$\Delta \delta$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title='./fig/p_ren.png';
saveas(gca,title,'png');



%[1] - Frequency Control Concerns In The North American Electric Power System 
%%

% 
% x = 0:1:24;
% %gaussmf(x,[10*3600 12*3600])
% nominal_load_increase_profile = gaussmf(x,[4 12])*0.2 + 1;
% 
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% plot(x,nominal_load_increase_profile,'LineWidth',1.5);
% ylabel('$ P_{L}$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
% hold off
% set(gcf,'renderer','Painters');
% title='./fig/load.png';
% saveas(gca,title,'png');