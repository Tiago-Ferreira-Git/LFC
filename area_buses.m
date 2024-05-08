% Data file = d2m_pwrmod1.m
clearvars -except path; close all; clc;
run('plot_options');

%% Run nonlinear simulation and store results
%cd('../pstess/');
addpath('../../analysis/')
addpath('../../pstess/')
set_path('data2m.m');
%set_path('IEEE_bus_118_no_res');  % set the working directory in Sget_path
%set_path('d2asb.m');

run('s_simu');  % run PSTess

clearvars -except g bus t S1

base_mva = 100;
%%

% area_1.bus = [ 1 2 117 3 12 4 11 13 5 7 6 16 14 33 15 30 8 9 10 38];
% 
% area_2.bus = [40 41 39 36 35 37 34 43 44 45 46 48 49 47 42 53 52 50 51 67 66 54 57 56 55 63 58 64 61 60 62 65 59];
% 
% area_3.bus = [17 18 19 20 21 22 23 24 25 26 27 28 29 31 32 113 114 115 72 ];
% 
% area_4.bus = [71 73 69 68 81 116 78 79 74 75 118 76 77 80 98 82 96 97 70];
% 
% area_5.bus = [85 84 83 95 88 89 92 94 93 86 87 90 91 102 101 ];
% 
% area_6.bus = [99 100 107 105 103 106 111 110 109 112 108 104];

area_1.bus = [ 1 10 2 20 3 4 101];

area_2.bus = [ 13 14 120 11 110 12];



%mac_areas = [12 4 10 40 42 61 59 49 66 65 46 31 113 72 24 26 27 73 69 116 89 90 87 91 100 107 103 111 112]'
area_1.inertia = 0;
area_2.inertia = 0;
% area_3.inertia = 0;
% area_4.inertia = 0;
% area_5.inertia = 0;
% area_6.inertia = 0;
% 
area_1.damping = 0;
area_2.damping = 0;
% area_3.damping = 0;
% area_4.damping = 0;
% area_5.damping = 0;
% area_6.damping = 0;
% 
area_1.machines = 0;
area_2.machines = 0;
% area_3.machines = 0;
% area_4.machines = 0;
% area_5.machines = 0;
% area_6.machines = 0;
% 
% 
area_1.tg_con = [];
area_2.tg_con = [];
% area_3.tg_con = [];
% area_4.tg_con = [];
% area_5.tg_con = [];
% area_6.tg_con = [];


area_1.tg_sig = [];
area_2.tg_sig = [];
% area_3.tg_sig = [];
% area_4.tg_sig = [];
% area_5.tg_sig = [];
% area_6.tg_sig = [];


area_1.mac_base = [];
area_2.mac_base = [];
% area_3.mac_base = [];
% area_4.mac_base = [];
% area_5.mac_base = [];
% area_6.mac_base = [];

g.mac.mac_con(:,16) = g.mac.mac_con(:,16).*g.mac.mac_con(:,3)/base_mva;
g.mac.mac_con(:,17) = g.mac.mac_con(:,17).*g.mac.mac_con(:,3)/base_mva;


%% obtaining area inertia and damping
for i=1:118
    
    mask = g.mac.mac_con(:,2) == i;
    if any(mask)
        
        if ismember(i, area_1.bus)
    
            area_1.inertia = area_1.inertia + g.mac.mac_con(mask,16)*g.mac.mac_con(mask,3)/base_mva;
            area_1.damping = area_1.damping + g.mac.mac_con(mask,17);
            area_1.machines = area_1.machines + 1;
            area_1.tg_con = [area_1.tg_con g.tg.tg_con(g.mac.mac_con(mask,1),:)];
            area_1.tg_sig = [area_1.tg_sig g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
            area_1.mac_base = [area_1.mac_base g.mac.mac_con(mask,3)];

        elseif ismember(i, area_2.bus)

            area_2.inertia = area_2.inertia + g.mac.mac_con(mask,16);
            area_2.damping = area_2.damping + g.mac.mac_con(mask,17);
            area_2.machines = area_2.machines + 1;
            area_2.tg_con = [area_2.tg_con g.tg.tg_con(g.mac.mac_con(mask,1),:)];
            area_2.tg_sig = [area_2.tg_sig g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
            area_2.mac_base = [area_2.mac_base g.mac.mac_con(mask,3)];

        elseif ismember(i, area_3.bus)

            area_3.inertia = area_3.inertia + g.mac.mac_con(mask,16);
            area_3.damping = area_3.damping + g.mac.mac_con(mask,17);
            area_3.machines = area_3.machines + 1;
            area_3.tg_con = [area_3.tg_con g.tg.tg_con(g.mac.mac_con(mask,1),:)];
            area_3.tg_sig = [area_3.tg_sig g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
            area_3.mac_base = [area_3.mac_base g.mac.mac_con(mask,3)];
    
        elseif ismember(i, area_4.bus)

            area_4.inertia = area_4.inertia + g.mac.mac_con(mask,16);
            area_4.damping = area_4.damping + g.mac.mac_con(mask,17);
            area_4.machines = area_4.machines + 1;
            area_4.tg_con = [area_4.tg_con g.tg.tg_con(g.mac.mac_con(mask,1),:)];
            area_4.tg_sig = [area_4.tg_sig g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
            area_4.mac_base = [area_4.mac_base g.mac.mac_con(mask,3)];
    
        elseif ismember(i, area_5.bus)

            area_5.inertia = area_5.inertia + g.mac.mac_con(mask,16);
            area_5.damping = area_5.damping + g.mac.mac_con(mask,17);
            area_5.machines = area_5.machines + 1;
            area_5.tg_con = [area_5.tg_con g.tg.tg_con(g.mac.mac_con(mask,1),:)];
            area_5.tg_sig = [area_5.tg_sig g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
            area_5.mac_base = [area_5.mac_base g.mac.mac_con(mask,3)];
    
        elseif ismember(i, area_6.bus)

            area_6.inertia = area_6.inertia + g.mac.mac_con(mask,16);
            area_6.damping = area_6.damping + g.mac.mac_con(mask,17);
            area_6.machines = area_6.machines + 1;
            area_6.tg_con = [area_6.tg_con g.tg.tg_con(g.mac.mac_con(mask,1),:)];
            area_6.tg_sig = [area_6.tg_sig g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
            area_6.mac_base = [area_6.mac_base g.mac.mac_con(mask,3)];

        end
    end

end

area_1.damping = area_1.damping/area_1.machines;
area_2.damping = area_2.damping/area_2.machines;
% area_3.damping = area_3.damping/area_3.machines;
% area_4.damping = area_4.damping/area_4.machines;
% area_5.damping = area_5.damping/area_5.machines;
% area_6.damping = area_6.damping/area_6.machines;



%% Retrieving lmod index
[~,~,index_lmod] = intersect(area_1.bus,g.lmod.lmod_con(:,2),'stable');
area_1.lmod = index_lmod;

[~,~,index_lmod] = intersect(area_2.bus,g.lmod.lmod_con(:,2),'stable');
area_2.lmod = index_lmod;

% [~,~,index_lmod] = intersect(area_3.bus,g.lmod.lmod_con(:,2),'stable');
% area_3.lmod = index_lmod;
% 
% [~,~,index_lmod] = intersect(area_4.bus,g.lmod.lmod_con(:,2),'stable');
% area_4.lmod = index_lmod;
% 
% [~,~,index_lmod] = intersect(area_5.bus,g.lmod.lmod_con(:,2),'stable');
% area_5.lmod = index_lmod;
% 
% [~,~,index_lmod] = intersect(area_6.bus,g.lmod.lmod_con(:,2),'stable');
% area_6.lmod = index_lmod;


%% Getting lines that connect to other areas

area_1.to_bus = [];
area_2.to_bus = [];
% area_3.to_bus = [];
% area_4.to_bus = [];
% area_5.to_bus = [];
% area_6.to_bus = [];

areas = [area_1; area_2; ];
    %area_3; area_4; area_5; area_6];

neighbour = 1:2;

for j = 1:size(areas)
    mask = neighbour == j;
    to_areas = neighbour(~mask);
    for i=1:size(areas(j).bus,2)
        
        mask = g.line(:,1) == areas(j).bus(i);

        if any(ismember(g.line(mask,2), areas(j).bus) == 0 ) 

            for k = 1:length(to_areas)
                mask_ismember  = ismember(g.line(mask,2), areas(to_areas(k)).bus);
                if any(mask_ismember)
                    to_line = g.line(mask,:);
                    to_line = to_line(mask_ismember,:);
                    areas(j).to_bus = [areas(j).to_bus; ones(size(to_line,1),1)*to_areas(k) to_line];
                end
                
            end
        end
        
        mask = g.line(:,2) == areas(j).bus(i);

        if any(ismember(g.line(mask,1), areas(j).bus) == 0 ) 

            for k = 1:length(to_areas)
                mask_ismember  = ismember(g.line(mask,1), areas(to_areas(k)).bus);
                if any(mask_ismember)
                    to_line = g.line(mask,:);
                    to_line = to_line(mask_ismember,:);
                    areas(j).to_bus = [areas(j).to_bus; ones(size(to_line,1),1)*to_areas(k) to_line];

                    areas(j).to_bus(:,[2,3]) = areas(j).to_bus(:,[3,2]);
                end
                
            end
        end
   
    end
end



%% Generate the ss for each area
bus_ss = [];
A_global = [];
B_global = [];
C_global = [];

A_area = [];
B_area = [];
C_area = [];
base_mva = 100;
u = [];
s = tf('s');

for i=1:length(areas)
    u_area = [];
    
    if( 0 < size(areas(i).lmod,1))
        u_area = [u_area ; -sum(g.lmod.lmod_sig(areas(i).lmod,:),1)];
    else
        u_area = [u_area ; zeros(1,size(g.lmod.lmod_sig(1,:),2))];
    end


    for j = 1:size(areas(i).tg_con,1)

        %p_mech ss
        
        Ts = areas(i).tg_con(j,6); Tc = areas(i).tg_con(j,7); T3 = areas(i).tg_con(j,8); T4 = areas(i).tg_con(j,9); T5 = areas(i).tg_con(j,10);
    
        
        sys = (1/(1 + s*Ts)) * ((1+s*T3)/(1+s*Tc)) * ((1+s*T4)/(1+s*T5));
        %sys = sys*areas(i).mac_base(j)/base_mva;
        [A_mech,B_mech,C_mech,D_mech] = tf2ss(sys.Numerator{1},sys.Denominator{1});
        
        %Control signal input
        u_area = [u_area ;  areas(i).tg_sig(j,:)];

       
        A_area = [A_area zeros(size(A_area,1),size(A_mech,2)) ; zeros(size(A_mech,1),size(A_area,2)) A_mech];
        B_area = [B_area zeros(size(B_area,1),size(B_mech,2)) ; zeros(size(B_area,1),size(B_mech,2)) B_mech];
        C_area = [C_area C_mech];

    end


    intertia_sys = 1/(s*areas(i).inertia + areas(i).damping);
    %/base_mva
    [A_freq,B_freq,C_freq,D_freq] = tf2ss(intertia_sys.Numerator{1},intertia_sys.Denominator{1});
    
    %we can concactenate it like this because B is one dimensional to 1
    A = [A_freq B_freq.*C_mech ; -areas(i).tg_con(:,4).*B_mech A_mech];
    

    %inputs of this system will be p_load and p_c - tg_sig | the last
    %row is to null the effect of input to the p_tie state
    B = [B_freq 0 ; zeros(size(B_mech,1),1)  B_mech ; 0 0 ];

    %Add tie-line state
    A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
    A(1,end) = -B_freq;
    

    C = zeros(2,size(A,1));
    C(1,1) = 1;
    C(2,end) = 1;

    u = [u;u_area];
    
    bus_ss = [bus_ss; g.bus.bus_int(i) size(A,1)];
    A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
    B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
    C_global = [C_global zeros(size(C_global,1), size(C,2)); zeros(size(C,1), size(C_global,2)) C];
end



tol = 1e-8;   % tolerance for convergence
iter_max = 30; % maximum number of iterations
acc = 1.0;   % acceleration factor
[bus_sol,line,line_flw] = loadflow(bus,g.line,tol,iter_max,acc,'n',2);
bus_sol(:,3) = deg2rad(bus(:,3));


%Add tie-lines

for i = 1:size(areas)
    

    neighbours = areas(i).to_bus;

    T_i = 0;
    for j = 1:size(neighbours)
        T_ji = 377*cos(bus_sol(g.bus.bus_int(neighbours(j,3)),3) - bus_sol(g.bus.bus_int(neighbours(j,2)),3))/(neighbours(j,5));
        %377*
        %T_ji = T_ji * gain_tie_lines(i);
        A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:neighbours(j,1)-1,2))+1 ) = -T_ji;
        %A_global(sum(bus_ss(1:i,2))+1 ,sum(bus_ss(1:i,2)) ) = T_ji;


        T_i =  T_i + T_ji;
    end

    A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:i-1,2))+1 ) = T_i;

end

D_global = zeros(size(C_global,1),size(B_global,2));
sys = ss(A_global,B_global,C_global,D_global);
% sys = c2d(sys,0.01);
%eig(A_global);

%%
V1 = g.bus.bus_v(g.bus.bus_int(g.line(:,1)),:);
V2 = g.bus.bus_v(g.bus.bus_int(g.line(:,2)),:);
% V1 = ones(size(V1)).*exp(1j.*angle(V1));
% V2 = ones(size(V2)).*exp(1j.*angle(V2));
R = g.line(:,3);
X = g.line(:,4);
B = g.line(:,5);
tap = g.line(:,6);
phi = g.line(:,7);
[S1,S2] = line_pq(V1,V2,R,X,B,tap,phi);


S1 = real(S1) - real(S1(:,1));
S2 = real(S2) - real(S2(:,1));
g.freq.bus_freq = g.freq.bus_freq - 1;
g.mac.pelect = g.mac.pelect-g.mac.pelect(:,1);
g.mac.pmech = g.mac.pmech-g.mac.pmech(:,1);
%%

% %plot real tie-lines in comparison with aproximation
tie_lines = zeros(size(S1(1,:),2),length(g.bus.bus_int));
tie_lines_approx = zeros(size(S1(1,:),2),length(g.bus.bus_int));


for i=1:length(g.bus.bus_int)
    pf = zeros(size(S1(1,:),2),1);
    approx = zeros(size(S1(1,:),2),1);
    
    %Get power income
    mask = g.line(:,2) ==  g.bus.bus_int(i,1);
    
    if(any(mask))
        if(size(S1(mask,:)',2) ~= 1)
            pf = pf - sum(S1(mask,:))';
        else
            pf = pf - (S1(mask,:))';
        end
    end
    
    
    %Subtract power flow
    
    mask = g.line(:,1) ==  g.bus.bus_int(i,1);

    if(any(mask))
        if(size(S1(mask,:)',2) ~= 1)
            pf = pf + sum(S1(mask,:))';
        else
            pf = pf + (S1(mask,:))';
        end
    end
    if i == 1 pf_ = pf;end
    tie_lines(:,i) = pf;
end
%If tie-lines are wrong, experiment with Vbuses at 1pu



mask = (t < 200);
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




figure
hold on;
box on;
set(gca,'TickLabelInterpreter','latex') % Latex style axis
plot(t,tie_lines(:,g.bus.bus_int(101)),t,tie_lines(:,g.bus.bus_int(13)));
legend('$\Delta P_{tie_1}$','$\Delta P_{tie_2}$','$\Delta P_{tie_3}$','$\Delta P_{tie_4}$')
xlabel('Time - [s]')
ylabel('power flow per unit');
hold off;



figure
hold on
box on;
set(gca,'TickLabelInterpreter','latex') % Latex style axis
grid on;
plot(t,y(:,2:2:end))
set(gca,'TickLabelInterpreter','latex') % Latex style axis
legend('$\Delta P_{tie_1}$','$\Delta P_{tie_2}$','$\Delta P_{tie_3}$','$\Delta P_{tie_4}$')
xlabel('Time - [s]')
ylabel('power flow per unit');
hold off



%%
figure
hold on
plot(t,g.lmod.lmod_sig)
xlabel('time in seconds')
ylabel('power in pu on system base')
legend({'$P_{L_1}$','$P_{L_2}$','$P_{L_3}$','$P_{L_4}$','$P_{L_5}$'},'Interpreter','latex')
title('load electrical power')




figure
hold on
plot(t,g.mac.pmech*991/100)
xlabel('time in seconds')
ylabel('power in pu on system base')
legend({'$Pm_1$','$Pm_2$','$Pm_3$','$Pm_4$','$Pm_5$'},'Interpreter','latex')
title('turbine power')


figure
hold on
plot(t,g.mac.pmech*991/100)
xlabel('time in seconds')
ylabel('power in pu on system base')
legend({'$Pm_1$','$Pm_2$','$Pm_3$','$Pm_4$','$Pm_5$'},'Interpreter','latex')
title('turbine power')



%%
figure
hold on
plot(t,g.mac.pelect)

xlabel('time in seconds')
ylabel('power in pu on system base')
legend({'$Pe_1$','$Pe_2$','$Pe_3$','$Pe_4$','$Pe_5$'},'Interpreter','latex')
title('generator electrical power')



figure
hold on
plot(t,g.tg.tg_sig)

xlabel('time in seconds')
ylabel('power in pu on system base')
legend({'$Pc_1$','$Pc_2$','$Pc_3$','$Pc_4$','$Pc_5$'},'Interpreter','latex')
title('Control signal')



