
clearvars -except path; close all; clc;
run('../plot_options');
addpath('../../pstess/')

[g,bus,line] = get_g('data3');
<<<<<<< HEAD
%[g,bus,line] = get_g('IEEE_bus_118');
%[g,bus,line] = get_g('AutoSynGrid_3000');
=======
%[g,bus,line] = get_g('AutoSynGrid_3000');
[g,bus,line] = get_g('IEEE_bus_118');

>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
% bus(bus(:,10)==2,11) = 0;
% bus(bus(:,10)==2,12) = 0;

tol = 1e-8;                      % tolerance for convergence
itermax = 50;                    % maximum number of iterations
acc = 1.0;                       % acceleration factor


[bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);

bus_initial = bus;
clearvars -except g bus t line bus_initial

flag_integrator = 0;
flag_ren = 0;
flag_plot_metrics = 0;


n_areas = 250;
<<<<<<< HEAD
n_areas = 3;
[A,B,C,D,W,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,n_areas,0,flag_ren,flag_integrator);


cond(A)
A_c = A;
%%
h = 0.5;

[A,B,W] = discrete_dynamics(A,B,W,h);
rank(ctrb(A,B))

% G = ss(A,B,C,D,h);
% Gn = ssbal(G);
% 
% A = Gn.A;
% B = Gn.B;
% C = Gn.C;
% D = Gn.D;

W_ = permute_matrix(A,ren_ss);
if isempty(ren_ss)
   W_ = eye(size(A,1));
end



%number of controlable nodes
n_C = size(A,1) - size(ren_ss,2);
=======
%n_areas = 20;
n_areas = 30;

flag_integrator = 0;


[A_global,B_global,C,D,W_global,E,~,g,areas,network,bus_ss,rows_NC] = get_global_ss(g,bus,n_areas,0,1,flag_integrator);
n_ren = size(rows_NC,2);
B = B_global;
check_rank = [B A_global*B (A_global^2)*B (A_global^3)*B (A_global^4)*B_global (A_global^5)*B_global (A_global^6)*B_global (A_global^7)*B_global (A_global^7)*B_global (A_global^8)*B_global (A_global^9)*B_global (A_global^10)*B_global] ;

%(A_global^4)*B_global (A_global^5)*B_global (A_global^6)*B_global (A_global^7)*B_global (A_global^7)*B_global (A_global^8)*B_global (A_global^9)*B_global (A_global^10)*B_global


n = sum(sum(check_rank==0,2)==size(check_rank,2))
n_ren
h = 2.5;

[A,B,W] = discrete_dynamics(A_global,B_global,W_global,h);
W_ = permute_matrix(A,rows_NC);
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16

%plot_network(areas,line,n_areas);

%%
tic


simulation_hours = 5;
t = 0:h:3600*simulation_hours;
%7*24*3600;

w = zeros(size(W,2),size(t,2));
x0 = zeros(size(A,1),1);

mask = t > 30;

<<<<<<< HEAD
w = get_disturbance_profile(w,h,n_areas,simulation_hours,bus_ss);
=======
w = get_disturbance_profile(w,h,n_areas,simulation_hours,g.n_ren);
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16

R_ = 0.1;

<<<<<<< HEAD
q = zeros(1,size(A,1));
q(1) = 1;
q(1,cumsum(bus_ss(1:end-1,2))+1) = 1;


if flag_integrator
    % weight of the integrator 
    q(1,cumsum(bus_ss(:,2))-1) = 10;
end

q(1,cumsum(bus_ss(1:end-1,2))) = 100;
=======
if flag_integrator
    q = sum(C(1:4:end,:),1);
else
    q = sum(C(1:3:end,:),1);
end
if flag_integrator
    q = q + -10*sum(C(3:4:end,:),1);
end


>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16

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
<<<<<<< HEAD
    K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
    %K = get_gain(A,B,E,R_,q);
    % if n_C == size(A,1)
    %     %K = get_gain(A,B,E,R_,q);
    %     K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
    % else
    %     K = get_gain(A,B,E,R_,q,W_,n_C);
    % end
=======
    K = get_gain(A,B,E,R_,q,W_,size(rows_NC,2));
    %K = dlqr(A,B,diag(q),R_*eye(size(B,2)));
    %K = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
    %K = zeros(size(B,2),size(A,1));
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
    if isnan(K)
        toc
        error 'Could not compute Gains'
    end
    
    hour = 1;
    day = 0;
    flag = 1;
    % t_settling = 0;
    bus = bus_initial;
    k_ = 1:1440:length(t);
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
                 %K = get_gain(A,B,E,R_,q);
           end
           
           %[A_global,bus,k_tie] = update_dynamics(areas,bus,bus_initial,line,n_areas,A_global,bus_ss,g,network,hour-day*24+1);
           %[A,B,~,W,~,~] = discrete_dynamics(A_global,B_global,C_global,D_global,W_global,E_global,bus_ss,h,n_areas,L);
            

           
            %to_plot(hour,:) = k_tie;
        end

        u(:,k) = -K*x(:,k);
        u(:,k) = min(max(u(:,k),-0.1),0.1);
        x(:,k+1) = A*x(:,k) + B*u(:,k) + W*w(:,k);

        y(:,k+1) = C*x(:,k+1);



<<<<<<< HEAD
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
=======
        
        % y(1:4:end-n_areas,k) = min(max(y(1:4:end-n_areas,k),-freq_limit),freq_limit);
        % 
        % %Controller performance metric
        % %all(abs(y(1:3:end-n_areas,k+1)) < 1e-9)
        % if all(abs(y(1:3:end-n_areas,k+1)) < 1e-9) && flag
        %     t_settling = ((k-k_(hour))*h );
        %     time_settling(1,hour) =  t_settling;
        %     flag = 0;
        % 
        % end
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16


        %Get the control action per area
        % for i=1:n_areas
        %     u_area(i,k) = sum(u(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+bus_ss(i,3),3));
        % end

    end

<<<<<<< HEAD
    disp_cost_area(control_type,:) = sum(u_area(1:end,:),2);
    disp_cost_machine(control_type,:) = sum(u(1:end,:),2);
    if flag_integrator
        frequency_error_cost(control_type,:) = sum(abs(y(1:4:end,:)),2)';
    else
        frequency_error_cost(control_type,:) = sum(abs(y(1:3:end,:)),2)';
    end
    time_settling_cost(control_type,:) = time_settling;
=======
    % disp_cost_area(control_type,:) = sum(u_area(1:end,:),2);
    % disp_cost_machine(control_type,:) = sum(u(1:end,:),2);
    % frequency_error_cost(control_type,:) = sum(abs(y(1:3:end-n_areas,:)),2)';
    % time_settling_cost(control_type,:) = time_settling;
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
    y = y';


   
end


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

<<<<<<< HEAD
if flag_integrator 
=======
if flag_integrator
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
    stairs(t,y(:,1:4:end),'LineWidth',1.5);
    ylim([min(min(y(:,1:4:end)))*1.3,max(max(y(:,1:4:end)))*1.3])
else
    stairs(t,y(:,1:3:end),'LineWidth',1.5);
    ylim([min(min(y(:,1:3:end)))*1.3,max(max(y(:,1:3:end)))*1.3])
end
<<<<<<< HEAD
yline(freq_limit,'--');
yline(-freq_limit,'--');
=======

yline(freq_limit,'--');
yline(-freq_limit,'--');

>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
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
<<<<<<< HEAD
if flag_integrator 
=======
if flag_integrator
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
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




if flag_integrator
    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,y(:,3:4:end),'LineWidth',1.5);
    
    legend('$\Delta \delta_{1}$','$\Delta \delta_{2}$','$\Delta \delta_{3}$','Interpreter','latex')
    ylabel('$\Delta \delta$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/delta.png';
    saveas(gca,title,'png');
end

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
<<<<<<< HEAD
if flag_integrator 
    stairs(t,y(:,3:4:end),'LineWidth',1.5);
=======

if flag_integrator
    stairs(t,y(:,4:4:end),'LineWidth',1.5);
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
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

<<<<<<< HEAD
cond(A);
cond(A_c);

=======





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
>>>>>>> 7ed45c8baf499c9855c492ee0bd319b4b80e9d16
%%

% 
% % Define A1 and A_hat
% 
% 
% teste = [A_hat A_c];
% 
% [RB,p] = rref(teste);
% 
% P_1 = RB(1:end,size(A_hat,1)+1:end);
% 
% A_hat_calculated = P_1 * A_c / P_1;
% 
% 
% B = [-3 2; 4 -2];
% 
% B_ = [-1,2; 2 -2];
% 
% teste = [B_ B];
% 
% [RB,p] = rref(teste);


