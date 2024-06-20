
clearvars -except path; close all; clc;
run('../plot_options');
addpath('../../pstess/')

[g,bus,line] = get_g('data2m');
%[g,bus,line] = get_g('IEEE_bus_118_res');
%[g,bus,line] = get_g('AutoSynGrid_3000');
% bus(bus(:,10)==2,11) = 0;
% bus(bus(:,10)==2,12) = 0;

tol = 1e-8;                      % tolerance for convergence
itermax = 50;                    % maximum number of iterations
acc = 1.0;                       % acceleration factor


[bus,line,line_flw] = loadflow(bus,line,tol,itermax,acc,'n',2);

bus_initial = bus;
clearvars -except g bus t line bus_initial

flag_integrator = 1;
flag_ren = 1;
flag_plot_metrics = 0;


n_areas = 250;
n_areas = 1;
[A,B,C,D,W,~,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,n_areas,flag_ren,flag_integrator);

%% Testing controlability

A(end,1) = 1e3;
A(1,end) = W(1);

A12 = zeros(size(A));

A12(end,1) = -1e3;
A21 = A12;

A_ =  [A A12; A21 A];
B_= [B zeros(size(B));zeros(size(B)) B ];
C_ = zeros(1,size(A_,1));
C_(1,1) = 1;
%C_(2,6) = 1;

s = tf('s');

global_ = C_*inv(s*eye(size(A_)) - A_)*B_;

w_tie = [1 0 0 0 0]*inv(s*eye(size(A)) - A)*B;
% roots(w_tie.Numerator{1})
% roots(w_tie.Denominator{1})



wo_tie = [1 0 0 0] *inv(s*eye(size(A(1:4,1:4))) - A(1:4,1:4))*B(1:4);
roots(global_.Numerator{1,1});
roots(global_.Denominator{1,1});

roots(wo_tie.Numerator{1});
roots(wo_tie.Denominator{1});


P1_tie = -A(end,1)/(s*(s*network(1).inertia + network(1).damping));
P2_tie = P1_tie;

P1 = wo_tie;
P2 = P1;

%G1 = wo_tie*(s-P2_tie)/(s-P1_tie-P2_tie)
%G1 = wo_tie*(s)/(s-P1_tie)

G1 = -(P1*(1-P2_tie))/(1-P1_tie-P2_tie);

disp('Poles Calulated by Matlab of G1(s)')
roots(G1.Denominator{1})

disp('Poles Calulated by Matlab of P1(s)')
roots(P1.Denominator{1})
disp('Supposedly added Poles from tie-lines Calulated by Matlab')
roots([network(1).inertia network(1).damping 2e3])


disp('Poles Calulated by Matlab from state-space')
roots(global_.Denominator{1,1})

disp('Poles of ss (should be equal to G1(s))')
roots(global_.Denominator{1})


G2 = -(P1_tie*(P2))/(1-P1_tie-P2_tie);
roots(G2.Denominator{1})
roots(global_.Denominator{1,2})



-A(end,1)*A(1,end)
A_c = A;
%% Test

A_mech = A(2:4,2:4);
B_mech = B(2:4);
C_mech = [1 network(1).tg_con(1,8)*network(1).tg_con(1,9) 0];

Gmech = C_mech*inv(s*eye(size(A_mech)) - A_mech)*B_mech;
eig(A_mech)
eig(A(1:4,1:4))

Gteste = Gmech/(s*network(1).inertia+network(1).damping);
roots(Gteste.Denominator{1})
 


%%
h = 2.5;

[A,B,W] = discrete_dynamics(A,B,h,W);

rank(ctrb(A,B))

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

           if(24 <= hour -day*24 )
                 day = day + 1
                 %K = get_gain(A,B,E,R_,q);
           end
           
            %[~,bus,k_tie] = update_dynamics(areas,bus,bus_initial,line,n_areas,A,bus_ss,g,network,hour-day*24+1);
            %to_plot(hour,:) = k_tie;
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


title = sprintf('./fig/K_tie_%d_%d.eps',i,simulation_hours);
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(1:simulation_hours-5,to_plot(2:end-4,i),'LineWidth',1.5);
ylabel('$T_{{tie}_{i}}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
savefig('./fig/filename.fig');
set(gcf,'renderer','Painters');
saveas(gca,title,'epsc');




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

cond(A);
cond(A_c);

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


