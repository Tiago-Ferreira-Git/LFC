
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

[bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);

bus_initial = bus;
clearvars -except g bus t line bus_initial


n_areas = 30;
%n_areas = 3;
[A_global,B_global,C_global,D_global,W_global,~,E_global,L,areas,network,bus_ss] = get_global_ss(g,bus,n_areas,0,1);



h = 2.5;

[A,B,C,W,~,E] = discrete_dynamics(A_global,B_global,C_global,D_global,W_global,E_global,bus_ss,h,n_areas,L);
C_angle = zeros(n_areas,size(C,2));

C_angle(1:n_areas,end-n_areas+1:end) = -eye(n_areas);
C = [C ; C_angle];
%%
%plot_network(areas,line,n_areas);

%%
tic


simulation_hours = 24;
t = 0:h:3600*simulation_hours;
%7*24*3600;

w = zeros(size(W_global,2),size(t,2));
x0 = zeros(size(A,1),1);

mask = t > 30;

w = get_disturbance_profile(w,h,n_areas,simulation_hours);

R_ = 0.01;

q = zeros(1,size(A,1));
q(1,1) = 1;
for i=1:n_areas
    
    q(1,sum(bus_ss(1:i,2))+1) = 1;
end

% weight of the integrator 
q(1,size(A_global,1)+1:end) = 10;



%perfomance metrics
% time_settling = zeros(1,simulation_hours);
% frequency_error_cost = zeros(2,n_areas);
% disp_cost_machine = zeros(2,size(B,2));
% disp_cost_area = zeros(2,n_areas);
% time_settling_cost = zeros(2,simulation_hours);
to_plot = zeros(simulation_hours,n_areas);
bus_reactive_power_ren = zeros(132,simulation_hours);
bus_active_power_ren = zeros(132,simulation_hours);

bus_active_power_consumption = zeros(132,simulation_hours);
bus_reactive_power_consumption = zeros(132,simulation_hours);


%Frequency output limit according to [1] 
freq_limit = 0.05/50;

for control_type = 1:1
    % Controller gain synthesis 
    % if(control_type==1)
    %     decentralized = true;
    % else
    %     decentralized = false;
    % end
    % if ~decentralized  E = ones(size(E)) ; end
    %K = get_gain(A,B,E,R_,q);
    % if isnan(K)
    %     toc
    %     error 'Could not compute Gains'
    % end
    
    hour = 1;
    day = 0;
    flag = 1;
    % t_settling = 0;
    bus = bus_initial;
    k_ = 1:1440:length(t);
    % x = zeros(size(A,1),size(t,2));
    % u = zeros(size(B,2),size(t,2));
    % u_area = zeros(n_areas,size(t,2));
    % y = zeros(size(C,1),size(t,2));
    bus_reactive_power_ren(:,hour) = bus(:,5);
    t_settling = 0;
    x(:,1) = x0;
    for k = 1:length(t)-1
        
        if(k == k_(hour+1) )
           
           if hour == 12
                mask = ones(132,1);
                %bitand(bus(:,7)<0.5 , 0.3<bus(:,7));
           end
           
           hour = hour + 1;

           bus_active_power_ren(:,hour) = bus(:,4);
           bus_reactive_power_ren(:,hour) = bus(:,5);
           
            
           bus_active_power_consumption(:,hour) = bus(:,6);
           bus_reactive_power_consumption(:,hour) = bus(:,7);
           
           
           
           
           flag = 1;

           if(24 <= hour -day*24 )
                 day = day + 1
                 %K = get_gain(A,B,E,R_,q);
           end
           
           [A_global,bus,k_tie] = update_dynamics(areas,bus,bus_initial,line,n_areas,A_global,bus_ss,g,network,hour-day*24+1);
           [A,B,~,W,~,~] = discrete_dynamics(A_global,B_global,C_global,D_global,W_global,E_global,bus_ss,h,n_areas,L);
            

           
            to_plot(hour,:) = k_tie;
        end

        % u(:,k) = -K*x(:,k);
        % u(:,k) = min(max(u(:,k),-0.1),0.1);
        % x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);
        % 
        % y(:,k+1) = C*x(:,k+1);



        
        % y(1:3:end-n_areas,k) = min(max(y(1:3:end-n_areas,k),-freq_limit),freq_limit);
        % 
        % %Controller performance metric
        % %all(abs(y(1:3:end-n_areas,k+1)) < 1e-9)
        % if all(abs(y(1:3:end-n_areas,k+1)) < 1e-9) && flag
        %     t_settling = ((k-k_(hour))*h );
        %     time_settling(1,hour) =  t_settling;
        %     flag = 0;
        % 
        % end
        % 
        % 
        % %Get the control action per area
        % for i=1:n_areas
        %     u_area(i,k) = sum(u(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+bus_ss(i,3),3));
        % end

    end

    % disp_cost_area(control_type,:) = sum(u_area(1:end,:),2);
    % disp_cost_machine(control_type,:) = sum(u(1:end,:),2);
    % frequency_error_cost(control_type,:) = sum(abs(y(1:3:end-n_areas,:)),2)';
    % time_settling_cost(control_type,:) = time_settling;
    % y = y';


   
end
%%
covariances = zeros(n_areas,1);
means = zeros(n_areas,1);
for i =1:n_areas

    covariances(i) = cov(to_plot(2:end-3,i));
    means(i) = mean(to_plot(2:end-3,i));

end
mask = any((to_plot(2:end,:)-to_plot(2,:)) > 200);
title = sprintf('./fig/delta_K_tie_%d.png',i);
figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours-1,to_plot(2:end,mask)-to_plot(2,mask),'LineWidth',1.5);
ylabel('$T_{{tie}_{i}} - T_{{tie}_{i,0}}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
saveas(gca,title,'png');


title = sprintf('./fig/K_ties_%d.png',i);
figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours-1,to_plot(2:end,mask),'LineWidth',1.5);
ylabel('$T_{{tie}_{i}} $ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
saveas(gca,title,'png');

%%
figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(sqrt(covariances),'LineWidth',1.5);
%legend('$\sigma_1$','$\sigma_2}$','$\sigma_3}$','Interpreter','latex')
%y_label = sprintf('$T_{{tie}_{%d,k}} - T_{{tie}_{%d,0}}$',i,i);
ylabel('$\sigma$','interpreter','latex');
xlabel('$area$','Interpreter','latex');
xlim([0 35.2])
hold off
set(gcf,'renderer','Painters');
title = sprintf('./fig/cov_%d.eps',simulation_hours);
saveas(gca,title,'epsc');
title = sprintf('./fig/cov_%d.png',simulation_hours);
saveas(gca,title,'png');

 

figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
xlim([0 35.2])
plot(means,'LineWidth',1.5);
%legend('$\sigma_1$','$\sigma_2}$','$\sigma_3}$','Interpreter','latex')
%y_label = sprintf('$T_{{tie}_{%d,k}} - T_{{tie}_{%d,0}}$',i,i);
ylabel('$\mu$','interpreter','latex');
xlabel('$area$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = sprintf('./fig/mean_%d.eps',simulation_hours);
saveas(gca,title,'epsc');
title = sprintf('./fig/mean_%d.png',simulation_hours);
saveas(gca,title,'png');


figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,bus_reactive_power_ren(119:end,:),'LineWidth',1.5);
ylabel('$Q_G$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = './fig/Q_gen.png';
saveas(gca,title,'png');


figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,bus_active_power_ren(119:end,:),'LineWidth',1.5);
ylabel('$P_G$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = './fig/P_gen.png';
saveas(gca,title,'png');



%%
figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,bus_reactive_power_consumption,'LineWidth',1.5);
ylabel('$Q_L$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = './fig/Q_con.png';
saveas(gca,title,'png');



figure
set(gca,'FontSize',17);
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours,bus_active_power_consumption,'LineWidth',1.5);
ylabel('$P_L$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
hold off
set(gcf,'renderer','Painters');
title = './fig/P_con.png';
saveas(gca,title,'png');


return


%%
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
%%

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



%%
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