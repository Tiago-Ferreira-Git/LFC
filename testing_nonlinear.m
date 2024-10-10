clearvars -except path; close all; clc;
run('../plot_options');
run("D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Simulations\Power System Toolbox\matpower7.1\startup.m")
opt = mpoption('verbose',0,'out.all',0);


    

% 
% debug = 2;
% flag_ren = 0;
% flag_plot_metrics = 0;
% 
% 
% [mpc,n_res,idx_initial] = get_g('case14',flag_ren);
% idx = idx_initial;
% mpc = runopf(mpc,opt);
% mpc_initial = mpc;
% 
% 
% 
% clearvars -except mpc flag_plot_metrics flag_ren n_res idx idx_initial  simulation_hours simulation_seconds h debug opt mpc_initial
% 
% n_areas = 3;
% [A_c,B_c,C,~,W_c,~,C_mech,~,E,~,network,bus_ss,ren_ss,E_fs] = get_global_ss(mpc,n_areas,flag_ren,debug);
% if real(eig(A_c)) > 0
%     error 'Unstable Model'
% end

% figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
% hold on;
% grid on;
% box on;
% set(gca,'FontSize',20);
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% %%%%
% plot(real(eig(A_c)),imag(eig(A_c)),'x')
% %%%%
% % legend({'Temperature of agent $\nu_1$'},...
% % 	'Location','best','Interpreter','latex');
% ylabel('Imaginary Axis','Interpreter','latex');
% xlabel('Real Axis','Interpreter','latex');
% hold off;
% % Save figure to .fig and .eps formats
% savefig('./fig/poles118.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'./fig/poles11.eps','epsc');



% 
% max(A_c,[],'all')
% network_initial = network;
% save('data/sim_14_3')

%%%%%%save('data/sim_118_30')



% clear all
% clearvars -except path; close all; clc;
% 
% %load('data/sim_118_30')
% 
load('data/sim_14_3')
%%

h = 0.1;

simulation_hours = 1;
simulation_seconds = 0 + 3600*simulation_hours;


[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);


angle_index = cumsum(bus_ss(:,2));
freq_index = [1 ; angle_index+1];

% Initial conditions

[x0,u0,Pt0,PL0,Ploss]  = initial_conditions(size(A_c,1),size(B,2),bus_ss(:,2),network,mpc);
Pgen0 = C*x0;
Pgen0 = Pgen0(2:3:end,1);

teste = Pgen0 - (PL0 + Pt0 ); 

% Simulation Setup



t_L = 0:h:simulation_seconds;

w = zeros(size(W,2),size(t_L,2));

[w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss);


% for i = 1:n_areas
%     mpc.bus(network(i).bus(1),3) = mpc.bus(network(i).bus(1),3) + w_load(i,1)*100;
% end

% mpc = runopf(mpc,opt);

if flag_ren
    data = load('data\solar.mat');
    data = data.data;
    data.data(:,2:end) = data.data(:,2:end)./100;
    if ceil(simulation_seconds/3600) ~= 1
        P_res = resample(data.data(1:ceil(simulation_seconds/3600),2:end),3600/0.1,1,'Dimension',1);
        P_res = P_res(1:size(w_ren,2),:);
    else
        P_res = data.data(1,2:end);
    end
    % 
    P_res = P_res';
    P_res = P_res + w_ren;
    
else
    P_res = [];
end

% w_load = zeros(size(w_load));
P_load = PL0 + w_load;

tspan = [0 simulation_seconds];
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);

K = lqr(A_c,B_c,eye(size(A_c,1)),0.1*eye(size(B_c,2)));
K = zeros(size(K));


%Continuous linearized model 
delta_u = zeros(size(B,2),1);
%%
%[t_L,x_L] = ode45(@(t,x) linearized_model(t,x,network,bus_ss,x0,u0,w_load(:,1),P_res,delta_u,mpc.bus(:,8:9)),[0 simulation_seconds],x0,opts);

[t_L,x_L] = ode23(@(t,x) linearized_model_ss(t,x,bus_ss,A_c,B_c,W_c,x0,u0,w_load(:,1),delta_u,K),[0 simulation_seconds],x0,opt);



y_L = C*(x_L');
if any(isnan(x_L))
    return
end

% y_L(3:3:end,:) = mod(y_L(3:3:end,:) + pi, 2*pi) - pi;
%%



% Nonlinear Simulation
opts = odeset('RelTol',1e-8,'AbsTol',1e-8);





freq_limit = 0.05/50;
[t_nl,x_nL] = ode45(@(t,x) nonlinear_model(t,x,network,bus_ss,x0,u0,P_load(:,1),P_res,delta_u,mpc.bus(:,8:9),w_load(:,1),K),0:h:simulation_seconds,x0,opts);
y_nL = C*(x_nL');

% y_nL(3:3:end,:) = mod(y_nL(3:3:end,:) + pi, 2*pi) - pi;

if any(isnan(x_nL))
    return
end

figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_nl,y_nL(1:3:end,:)','LineWidth',1.5);
% yline(1+freq_limit,'--');
% yline(1-freq_limit,'--');
% title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
savefig('fig/Nonlinear.fig');
set(gcf,'renderer','Painters');
saveas(gca,'fig/Nonlinear.png','png');
hold off


%%

% 
% x_t = zeros(size(x_nL));
% 
% t_nl;
%     % 
%     % func = @(sigma) expm(A_c(area_index,area_index).*sigma);
%     % int = integral(func,0,h,'ArrayValued', true); 
% 
% for k = length(t_nl)
% 
%     func = @(s) expm(A_c*(t_nl(k)-s))*W_c*w_load(:,1);
%     x_t(k,:) = x0 + expm(A_c*t_nl(k))*(x0-x0) + integral(func,0,t_nl(k),'ArrayValued', true); 
% 
% 
% end
% 
% 
% 
% y_formula = C*(x_t');
% 
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_nl,y_formula(1:3:end,:)','LineWidth',1.5);
% % yline(1+freq_limit,'--');
% % yline(1-freq_limit,'--');
% % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% ylabel('$\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% savefig('fig/Nonlinear.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'fig/Nonlinear.png','png');
% hold off


%%

% W_c_ = W_c;
% 
% W_c_(1,1) = W_c(1,1).*w_load(1,1);
% W_c_(12,2) = W_c(12,2).*w_load(2,1);
% W_c_(17,3) = W_c(17,3).*w_load(3,1);
% 
% W_c_ = sum(W_c_,2);
% 
% sys = ss(A_c,W_c_,C(1:3:end,:),zeros(size(C(1:3:end,:),1),1));
% 
% step(sys)

%%

% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_L(1:3:end,:)','LineWidth',1.5);
% % yline(1+freq_limit,'--');
% % yline(1-freq_limit,'--');
% % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% ylabel('$\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% savefig('fig/ode45.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'fig/ode45.png','png');
% hold off
% 
% 
% 
% %%
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_nl,y_nL(3:3:end,:)','LineWidth',1.5);
% % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% ylabel('$\theta$ (rad)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
% set(gcf,'renderer','Painters');
% %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
% hold off
% 
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t_L,y_L(3:3:end,:)','LineWidth',1.5);
% % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% ylabel('$\theta$ (rad)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
% set(gcf,'renderer','Painters');
% %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
% hold off
% 
% 
% 
% 
% %%
% % figure
% % set(gca,'TickLabelInterpreter','latex') % Latex style axis
% % hold on
% % grid on
% % box on;
% % stairs(t_nl,-y_nL(3:3:end,:)','LineWidth',1.5);
% % % title(sprintf('t_{sh} = %.2f s',t_sh),'Interpreter','tex')
% % ylabel('$\theta$ (rad)','interpreter','latex');
% % xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% % %savefig(sprintf('./fig/f_t_{sh}_%.2f.fig',t_sh));
% % set(gcf,'renderer','Painters');
% % %saveas(gca,sprintf('./fig/f_t_{sh}_%.2f.png',t_sh),'png');
% % hold off
% 
% 
% %%
% figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
% hold on;
% grid on;
% box on;
% set(gca,'FontSize',20);
% set(gca,'TickLabelInterpreter','latex')
% 
% errors = (x_nL' - x_L');
% %error([5,10,15,20,25],:) = 0;
% for i = 1:size(A_c)
%     stairs(t_L,errors(i,:)','LineWidth',1.5);
% end
% ylabel('$\epsilon_{model}$ ','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% set(gcf,'renderer','Painters');
% savefig('fig/errors_integrator.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'fig/errors_integrator.png','png');
% saveas(gca,'./fig/errors_integrator.eps','epsc');
% 
% hold off
% 
% %%
% 
% figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
% hold on;
% grid on;
% box on;
% set(gca,'FontSize',20);
% set(gca,'TickLabelInterpreter','latex')
% errors(angle_index,:) = 0;
% 
% 
% for i = 1:size(A,1)
%     stairs(t_L,errors(i,:)','LineWidth',1.5);
% end
% ylabel('$\epsilon_{model}$ ','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% set(gcf,'renderer','Painters');
% savefig('fig/errors_wo_integrator.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'fig/errors_wo_integrator.png','png');
% saveas(gca,'./fig/errors_wo_integrator.eps','epsc');
% hold off




%% Discretization


% 
% A_d = zeros(size(A));
% B_d = zeros(size(B));
% W_d = zeros(size(W));
% 
% for i = 1:n_areas
% 
%     area_index = [network(i).mec_index(1)-1 network(i).mec_index network(i).mec_index(end)+1 ];
% 
% 
%     [A_d(area_index,area_index),B_d(area_index,network(i).u_index),W_d(area_index,i)] = discrete_dynamics(A_c(area_index,area_index),B_c(area_index,network(i).u_index),W_c(area_index,i),h);
% 
%     %teste = expm(A_c(area_index,area_index).*h);
%     func = @(sigma) expm(A_c(area_index,area_index).*sigma);
%     int = integral(func,0,h,'ArrayValued', true); 
% 
%     for j = 1:n_areas
%         if i ~= j
%             freq_index(j):angle_index(j);
%             A_neighbour = int*A_c(area_index,freq_index(j):angle_index(j));
% 
%             A_d(area_index,freq_index(j):angle_index(j)) = A_neighbour;
%         end
%     end
%     i;
% end
% 
% [~,I] = max(abs(A-A_d),[],'all')
% [~,I] = max(abs(B-B_d),[],'all')
% max(abs(W-W_d),[],'all')
% 
% 
% 
% 
% 
% 
% % Nonlinear Discrete simulation Method 1 - decentralized framework first
% 
% 
% t = 0:h:simulation_seconds; 
% 
% x_method1 = zeros(size(A,1),size(t,2));
% delta_u_method1  = zeros(size(B,2),size(t,2));
% y_method1  = zeros(size(C,1),size(t,2));
% x_method1(:,1) = zeros(size(A,1),1);
% 
% K = zeros(size(K));
% 
% 
% for k = 1:length(t)-1
% 
%     delta_u_method1(:,k) = -K*x_method1(:,k);
%     delta_u_method1(:,k) = min(max(delta_u_method1(:,k),-0.1),0.1);
% 
%     x_method1(:,k+1) = A_d*x_method1(:,k) + B_d*delta_u_method1(:,k)+ W_d*w(:,k);
%     y_method1(:,k+1) = C*(x_method1(:,k+1));
% 
% end


h = 0.1;

simulation_hours = 1;
simulation_seconds = 0 + 3600*simulation_hours;


[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);

x_method2 = zeros(size(A,1),length(t_nl));
delta_u_method2  = zeros(size(B,2),length(t_nl));
y_method2  = zeros(size(C,1),length(t_nl));
x_method2(:,1) = x0;
y_method2(:,1) = C*x0;
p_tie = zeros(2*n_areas,length(t_nl));


for k = 1:length(t_nl)-1

    delta_u_method2(:,k) = -K*x_method2(:,k);
    delta_u_method2(:,k) = min(max(delta_u_method2(:,k),-0.1),0.1);
    % 
    % for i = 1:n_areas
    %     p_tie(2*(i-1) + 1,k) = network(i).to_bus(:,end)'*(x_method2(angle_index(i),k) - x_method2(angle_index(network(i).to_bus(:,1)),k));
    % 
    % 
    %     p_tie(2*(i-1) + 2,k) = network(i).to_bus(:,end)'*sin(x_method2(angle_index(i),k) - x_method2(angle_index(network(i).to_bus(:,1)),k));
    % 
    % end
    
    x_method2(:,k+1) = x0 + A*(x_method2(:,k)-x0) + B*delta_u_method2(:,k)+ W*w(:,k);
    y_method2(:,k+1) = C*(x_method2(:,k+1));

end



% for i = 1:n_areas
%     figure
%     hold on
%     stairs(t,p_tie(2*(i-1)+1,:)' - p_tie(2*(i-1)+2,:)')
%     legend({'Approx','Real'})
% end


% freq_limit = 0.05/50;
% figure
% set(gca,'TickLabelInterpreter','latex') % Latex style axis
% hold on
% grid on
% box on;
% stairs(t,y_method1(1:3:end,:)','LineWidth',1.5);
% yline(1+freq_limit,'--');
% yline(1-freq_limit,'--');
% ylabel('$\omega$ (pu)','interpreter','latex');
% xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
% hold off
% set(gcf,'renderer','Painters');
% savefig('fig/discrete_method1.fig');
% set(gcf,'renderer','Painters');
% saveas(gca,'fig/discrete_method1.png','png');


freq_limit = 0.05/50;
figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
stairs(t_nl,y_method2(1:3:end,:)','LineWidth',1.5);
% yline(1+freq_limit,'--');
% yline(1-freq_limit,'--');
ylabel('$\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
set(gcf,'renderer','Painters');
savefig('fig/discrete_method2.fig');
set(gcf,'renderer','Painters');
saveas(gca,'fig/discrete_method2.png','png');


hold off






figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
hold on;
grid on;
box on;
set(gca,'FontSize',20);
set(gca,'TickLabelInterpreter','latex')

mask = ones(1,size(A_c,1));
mask(freq_index(1:end-1)) = 1;
mask(angle_index) = 0;
mask = logical(mask);


errors = (x_nL' - x_method2);
errors(~mask,:) = 0;
for i = 1:size(A_c)
    stairs(t_nl,errors(i,:)','LineWidth',1.5);
end
ylabel('$\epsilon_{discrete}$ ','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
set(gcf,'renderer','Painters');
savefig('fig/errors_discrete.fig');
set(gcf,'renderer','Painters');
saveas(gca,'fig/errors_discrete.png','png');
saveas(gca,'./fig/errors_discrete.eps','epsc');

hold off


