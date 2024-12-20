function [w,w_load,w_res,P_load,P_res,u0,P_load_forecasted,mpc_] = get_disturbance_profile(mpc,network,h,n_areas,simulation_seconds,bus_ss)
%get_disturbance_profile   Generate disturbance profiles for each area using real data.
%
%   Inputs:                          
%       
%       mpc - The MATPOWER case struct after running OPF.
%       network - Vector of area objects returned from the get_global_ss function.
%       h - Sampling time in seconds.
%       n_areas - Number of areas in the network.
%       simulation_seconds - Total simulation duration in seconds.
%       bus_ss - A matrix that contains useful information about the state
%           sapce representation.
%
%
%   Outputs:
%
%       w - Concactenated vector of w_load and w_res. 
%       w_load - Incremental disturbance vector associated with loads.
%       w_res - Incremental disturbance vector associated with RESs.
%       P_load - Absolute real load power.
%       P_res - Absolute res power.
%       u0 - Initial control action (deprecated).
%       P_load_forecasted - Absolute forecasted load power.
%       mpc_ - The  MATPOWER case after running opf with initial changes in load
%       power / RESs.
    

    load('data\load_profile.mat');

    if (size(mpc.bus,1) == 118)
        load('data\res_profile_118.mat');
    else
        load('data\res_profile.mat');
        max_res = 200;
        res.forecast = res.forecast*max_res;
        res.measured = res.measured*max_res;
        
    end

    


    n_res = sum(mpc.isolar_mask);
    legend_text = cell(n_res,1);


    % Repeat the forecasted data to match the number of buses of the network or if t is specified to be greater than an year
    if size(load_.forecast,1) < length(mpc.bus) 
        load_.forecast = [repmat(load_.forecast, floor(length(mpc.bus)/size(load_.forecast,1)), 1); ...
                        load_.forecast(1:mod(length(mpc.bus), size(load_.forecast,1)),:)];
        load_.measured = [repmat(load_.measured, floor(length(mpc.bus)/size(load_.measured,1)), 1); ...
                load_.measured(1:mod(length(mpc.bus), size(load_.measured,1)),:)];
    else 
        load_.forecast = load_.forecast(1:length(mpc.bus),:);
        load_.measured = load_.measured(1:length(mpc.bus),:);
    end

    if size(res.forecast,1) < n_res 
        res.forecast = [repmat(res.forecast, floor(n_res/size(res.forecast,1)), 1); ...
                        res.forecast(1:mod(n_res, size(res.forecast,1)),:)];
        res.measured = [repmat(res.measured, floor(n_res/size(res.measured,1)), 1); ...
                res.measured(1:mod(n_res, size(res.measured,1)),:)];
    else 
        res.forecast = res.forecast(1:n_res,:);
        res.measured = res.measured(1:n_res,:);
    end



    simulation_hours = ceil(simulation_seconds/3600);


    n_machines = sum(mpc.isgen);
    
    plots_angles = zeros(size(mpc.bus,1),simulation_hours);
    plots_p = zeros(n_machines,simulation_hours);
    plots_p_load = zeros(size(mpc.bus,1),simulation_hours);
    area_load = zeros(n_areas,simulation_hours);
    area_load_measured = zeros(n_areas,simulation_hours);
    plots_p_res = zeros(n_res,simulation_hours);
    plots_tielines = zeros(n_areas,simulation_hours);

    u0 = zeros(n_machines,simulation_hours);

    P_load = zeros(n_areas,simulation_hours);
    P_load_forecasted = zeros(n_areas,simulation_hours);
    P_res = zeros(n_res,simulation_hours);
    w_load_hour = zeros(n_areas,simulation_hours);
    w_res_hour = zeros(n_res,simulation_hours);

    w = zeros(n_areas+n_res,simulation_hours);
    

    %Test if all load and renewables converge
    mpc_initial = mpc;
    j = 1;
    for k = 1:simulation_hours+1
        load_measured = zeros(size(mpc.bus,1),1);
        for i = 1:n_areas

            %Change loads
            %mpc_initial.bus(network(i).bus,3).*load_.forecast(network(i).bus,k);
            mpc.bus(network(i).bus,3) = mpc_initial.bus(network(i).bus,3).*load_.forecast(network(i).bus,k);
            load_measured(network(i).bus) = mpc_initial.bus(network(i).bus,3).*load_.measured(network(i).bus,k);
           
            
            area_load_measured(i,k) = sum(load_measured(network(i).bus));
            area_load(i,k) = sum( mpc.bus(network(i).bus,3));

            w_load_hour(i,k) = sum(load_measured(network(i).bus) - mpc.bus(network(i).bus,3))./mpc.baseMVA;

            
            P_load(i,k) = sum(load_measured(network(i).bus))./mpc.baseMVA;


            P_load_forecasted(i,k) = sum(mpc.bus(network(i).bus,3))./mpc.baseMVA;


            
            % Change RES
            idx = network(i).res_nr;

            if j ~= n_res+1
                legend_text(j) = cellstr(sprintf('Area %d',i));
                j = j+1;
            end

            mpc.gen(idx,2) = res.forecast(idx-n_machines,k);


            P_res(idx-n_machines,k) = sum(res.forecast(idx-n_machines,k))./mpc.baseMVA;

            w_res_hour(idx-n_machines,k)= sum(res.measured(idx-n_machines,k) - res.forecast(idx-n_machines,k))./mpc.baseMVA;


            %Active power limits
            mpc.gen(idx,9) = mpc.gen(idx,2) + 0.000001 ;
            mpc.gen(idx,10) = mpc.gen(idx,2) - 0.000001;
        
        
            %Reactive Power Limitis
            mpc.gen(idx,4) = min(1,res.forecast(idx-n_machines,k)./10)+0.001;
            mpc.gen(idx,4) = mpc.gen(idx,4);
        
            mpc.gen(idx,5) = 0;
            mpc.gen(idx,5) = mpc.gen(idx,5);
            
        end     
        [mpc,flag] = runopf(mpc,mpoption('verbose',0,'out.all',0));
        if ~flag
            error 'Load or Renewables profile made runopf not converge!'
        end

        [~,u0(:,k),~,~,Ploss]  = initial_conditions(sum(bus_ss(:,2)),n_machines,bus_ss(:,2),network,mpc);
        w_load_hour(:,k) = w_load_hour(:,k) + Ploss;
        P_load(:,k) = P_load(:,k) + Ploss;
        P_load_forecasted(:,k) = P_load_forecasted(:,k);
        
        plots_angles(:,k) = mpc.bus(:,9);
        plots_p(:,k) = mpc.gen(mpc.isgen,2);
        plots_p_load(:,k) = mpc.bus(:,3);
        plots_p_res(:,k) = mpc.gen(mpc.isolar_mask,2);
        plots_tielines(:,k) = compute_tielines(mpc,network);

        if k == 1
            mpc_ = mpc;
        end

    end

    %Plots
    if simulation_hours > 1
        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,deg2rad(plots_angles(:,1:end-1))','LineWidth',1.5);
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        ylabel('$\mathrm{Bus\;angles\;} (\mathrm{rad})$','Interpreter','latex');
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
        hold off;
        % Save figure to .fig and .eps formats
        savefig('./fig/angles.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/angles.eps','epsc');
        saveas(gca,'./fig/angles.png','png');




        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,plots_p(:,1:end-1)','LineWidth',1.5);
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        ylabel('$\mathrm{Generators\;Power\;} (\mathrm{MW})$','Interpreter','latex');
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
        hold off;
        % Save figure to .fig and .eps formats
        savefig('./fig/generator_power.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/generator_power.eps','epsc');
        saveas(gca,'./fig/generator_power.png','png');



        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,plots_p_res(:,1:end-1)','LineWidth',1.5);
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        ylabel('$\mathrm{RES\;Power\;} (\mathrm{MW})$','Interpreter','latex');
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
        hold off;
        % Save figure to .fig and .eps formats
        savefig('./fig/res_power.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/res_power.eps','epsc');
        saveas(gca,'./fig/res_power.png','png');
        



        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,plots_p_load(:,1:end-1)','LineWidth',1.5);
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        ylabel('$\mathrm{Forecasted\;Load\;Power\;} (\mathrm{MW})$','Interpreter','latex');
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
        hold off;
        %Save figure to .fig and .eps formats
        savefig('./fig/load.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/load.eps','epsc');
        saveas(gca,'./fig/load.png','png');




        % figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        % hold on;
        % grid on;
        % box on;
        % set(gca,'FontSize',20);
        % set(gca,'TickLabelInterpreter','latex') % Latex style axis
        % stairs(1:simulation_hours,plots_p_load_measured(:,1:end-1)','LineWidth',1.5);
        % % legend({'Temperature of agent $\nu_1$'},...
        % %     'Location','best','Interpreter','latex');
        % ylabel('$\mathrm{Measured\;Load\;Power\;} (\mathrm{MW})$','Interpreter','latex');
        % xlabel('$t (\mathrm{h})$','Interpreter','latex');
        % hold off;
        % %Save figure to .fig and .eps formats
        % savefig('./fig/measured_load.fig');
        % set(gcf,'renderer','Painters');
        % saveas(gca,'./fig/measured_load.eps','epsc');


        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,plots_tielines(:,1:end-1)','LineWidth',1.5);
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        %legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','northeast','Interpreter','latex');
        ylabel('$T_{i}$','Interpreter','latex');
        xlabel('$t (\mathrm{h})$','Interpreter','latex');
        hold off;
        savefig('./fig/tielines.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/tielines.eps','epsc');
        saveas(gca,'./fig/tielines.png','png');


        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,w_load_hour(:,1:end-1)'.*mpc.baseMVA,'LineWidth',1.5);
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        ylabel('$\mathrm{Load\;Power\;Difference\;} (\mathrm{MW})$','Interpreter','latex');
        %legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','best','Interpreter','latex');
        xlabel('$t (\mathrm{h})$','Interpreter','latex');
        hold off;
        % Save figure to .fig and .eps formats
        savefig('./fig/w_load.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/w_load.eps','epsc');
        saveas(gca,'./fig/w_load.png','png');



        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,w_res_hour(:,1:end-1)'.*mpc.baseMVA,'LineWidth',1.5);
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        %legend(legend_text,'Location','best','Interpreter','latex');
        ylabel('$\mathrm{RES\;Power\;Difference\;} (\mathrm{MW})$','Interpreter','latex');
        xlabel('$t (\mathrm{h})$','Interpreter','latex');
        hold off;
        %Save figure to .fig and .eps formats
        savefig('./fig/w_res.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/w_res.eps','epsc');
        saveas(gca,'./fig/w_res.png','png');
        



        %  figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        % hold on;
        % grid on;
        % box on;
        % set(gca,'FontSize',20);
        % set(gca,'TickLabelInterpreter','latex') % Latex style axis
        % stairs(0:simulation_hours-1,area_load(:,1:end-1)','LineWidth',1.5);
        % xticks(0:1:(simulation_hours-1))
        % xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        % xlim([0 (simulation_hours-1)])
        % %legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','northeast','Interpreter','latex');
        % xlim([0 (simulation_hours-1)*3600])
        % ylabel('$\mathrm{Forecasted\;Load\;Power\;} (\mathrm{MW})$','Interpreter','latex');
        % xlabel('$t (\mathrm{h})$','Interpreter','latex');
        % hold off;
        % %Save figure to .fig and .eps formats
        % savefig('./fig/load.fig');
        % set(gcf,'renderer','Painters');
        % saveas(gca,'./fig/load.eps','epsc');




        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        stairs(0:simulation_hours-1,area_load_measured(:,1:end-1)','LineWidth',1.5);
        %legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','northeast','Interpreter','latex');
        xticks(0:1:(simulation_hours-1))
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)])
        ylabel('$\mathrm{Measured\;Load\;Power\;} (\mathrm{MW})$','Interpreter','latex');
        xlabel('$t (\mathrm{h})$','Interpreter','latex');
        hold off;
        %Save figure to .fig and .eps formats
        savefig('./fig/area_load.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/area_load.eps','epsc');
        saveas(gca,'./fig/area_load.png','png');

    end
   
    
    %% Resample to match the sampling time h
    load_mask = zeros(1,size(w,1));
    load_mask(1) = 1;
    load_mask(cumsum(1+bus_ss(1:end-1,4))+1) = 1;
    load_mask = logical(load_mask);

    res_mask = ~load_mask;


        

    u0 = interp1(0:3600:(simulation_hours)*3600,u0',0:h:simulation_seconds,'previous');
    u0 = u0';



    P_load = interp1(0:3600:(simulation_hours)*3600,P_load',0:h:simulation_seconds,'previous');
    P_load = P_load';
    if n_res ~= 0
        
        P_res = interp1(0:3600:(simulation_hours)*3600,P_res',0:h:simulation_seconds,'previous');
        if (size(P_res,1) ~=1)
            P_res = P_res';
        end
    end


    w_load = interp1(0:3600:(simulation_hours)*3600,w_load_hour',0:h:simulation_seconds,'previous');
    w_load = w_load';

    if n_res ~= 0
        w_res = interp1(0:3600:(simulation_hours)*3600,w_res_hour',0:h:simulation_seconds,'previous');
        if (size(w_res,1) ~=1)
            w_res = w_res';
        end
    end

    
    if n_res ~= 0
        w = zeros(n_areas+n_res,size(w_load,2));
        w(res_mask,:) = w_res;
    else
        w = zeros(n_areas,size(w_load,2));
        w_res = [];
    end
    w(load_mask,:) = w_load;
    
    % for i = 1:size(w,1)
    % 
    %     %Low Pass filter
    %     % freq_pole = 0.003;
    %     % a = [1 -exp(-h*freq_pole)];
    %     % b = freq_pole*h;
    %     % %-exp(-h*freq_pole);
    %     % w(i,:) = filter(b,a,w(i,:));
    % 
    % 
    % 
    %     windowSize = 3600/(h*2); 
    %     b = (1/windowSize)*ones(1,windowSize);
    %     a = 1;
    %     w(i,:) = filter(b,a,w(i,:));
    % 
    % 
    % 
    % end



end

