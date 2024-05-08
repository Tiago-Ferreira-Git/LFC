%this version uses k mean clustering to determine the area partitioning

% Data file = d2m_pwrmod1.m
clearvars -except path; close all; clc;
run('plot_options');

%% Run nonlinear simulation and store results
%cd('../pstess/');
addpath('../../analysis/')
addpath('../../pstess/')
%set_path('data16m.m');
set_path('IEEE_bus_14.m');
%set_path('IEEE_bus_118_no_res');  % set the working directory in Sget_path
%set_path('d2asb.m');

run('s_simu');  % run PSTess

clearvars -except g bus t S1

%%

[A_global,B_global,C_global,D_global,u] = get_global_ss(g,bus,3,1);

sys = ss(A_global,B_global,C_global,D_global);
%sys = c2d(sys,0.01)
g.freq.bus_freq = g.freq.bus_freq - 1;
g.mac.pelect = g.mac.pelect-g.mac.pelect(:,1);
g.mac.pmech = g.mac.pmech-g.mac.pmech(:,1);

%%

mask = (t < 10000);
%u(2,:) = abs(u(2,:));
[y,~] = lsim(sys,u(:,mask),t(mask));

t = t(mask);


figure
hold on
plot(t,u)
set(gca,'TickLabelInterpreter','latex') % Latex style axis
legend({'$u_1$','$u_2$','$u_3$','$u_4$'},'interpreter','latex')
ylabel('Input Signal (pu)','interpreter','latex');
xlabel('Time (s)','interpreter','latex')
hold off
set(gcf,'renderer','Painters');
title_=sprintf('P_c.png');
saveas(gca,title_,'png');

%%
figure
hold on
plot(t,g.freq.bus_freq(1:end-1,mask))
set(gca,'TickLabelInterpreter','latex') % Latex style axis
legend({'$\Delta\omega_1$','$\Delta\omega _2$','$\Delta\omega _3$'},'interpreter','latex')
ylabel('$\Delta\omega$ (pu)','interpreter','latex');
xlabel('Time (s)','interpreter','latex')
hold off
title_=sprintf('freq.png');
saveas(gca,title_,'png');





figure
set(gca,'TickLabelInterpreter','latex') % Latex style axis
hold on
grid on
box on;
plot(t,y(:,1:2:end));
legend('$\Delta\omega_1$','$\Delta\omega_2$','$\Delta\omega_3$','Interpreter','latex')
ylabel('$\Delta\omega$ (pu)','interpreter','latex');
xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
hold off
title_=sprintf('freq_sim.png');
saveas(gca,title_,'png');


