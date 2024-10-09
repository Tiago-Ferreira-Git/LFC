function [w,w_load_,w_res,P_load_,P_res] = get_disturbance_profile(mpc,network,h,n_areas,simulation_seconds,PL0,n_res,res_buses)
    % w_load__hour = load_('data\w_load_.mat');
    % w_load__hour = w_load__hour.w_load_(:,1:24);
    % 
    % w_res_hour = load_('data\w_res.mat');
    % w_res_hour = w_res_hour.w_res(:,1:24);



    load('data\load_profile.mat');
    load('data\res_profile.mat');






    % n_res = sum(bus_ss(:,4));
    %if the number of areas in the data available is lower than the areas 
    if size(load_.forecast,1) < n_areas 
        load_.forecast = [repmat(load_.forecast, floor(n_areas/size(load_.forecast,1)), 1); ...
                        load_.forecast(1:mod(n_areas, size(load_.forecast,1)),:)];
        load_.measured = [repmat(load_.measured, floor(n_areas/size(load_.measured,1)), 1); ...
                load_.measured(1:mod(n_areas, size(load_.measured,1)),:)];
    else 
        load_.forecast = load_.forecast(1:n_areas,:);
        load_.measured = load_.measured(1:n_areas,:);
    end

    if size(res.forecast,1) < n_areas 
        res.forecast = [repmat(res.forecast, floor(n_areas/size(res.forecast,1)), 1); ...
                        res.forecast(1:mod(n_areas, size(res.forecast,1)),:)];
        res.measured = [repmat(res.measured, floor(n_areas/size(res.measured,1)), 1); ...
                res.measured(1:mod(n_areas, size(res.measured,1)),:)];
    else 
        res.forecast = res.forecast(1:n_areas,:);
        res.measured = res.measured(1:n_areas,:);
    end

    simulation_hours = ceil(simulation_seconds/3600);



    %Test if all load and renewables converge
    mpc_initial = mpc;
    for k = 1:simulation_hours
        for i = 1:n_areas
            mpc.bus(network(i).bus,3) = mpc_initial.bus(network(i).bus,3)*load_.forecast(i,k);
    
        end     
        [~,flag] = runopf(mpc);
        if ~flag
            error 'Load or Renewables profile made runopf not converge!'
        end
    end
    
    

    load__mask = zeros(1,size(w,1));
    load__mask(1) = 1;
    load__mask(cumsum(1+bus_ss(1:end-1,4))+1) = 1;
    load__mask = logical(load__mask);

    res_mask = ~load__mask;


    k_ = 1:3600/h:simulation_hours*3600/h;

    for k=1:simulation_seconds/h + 1
        
        if ~isempty(find(k == k_, 1))
            hour = find(k == k_);
        end
        w(load__mask,k) = w_load__hour(:,hour)*100;
        w(res_mask,k) = w_res_hour(:,hour);
        

        %w(load__mask,k_(hour):k_(hour+1)) = repmat(w_load__hour(:,hour),1,3600/h +1).*100;
        %w(res_mask,k_(hour):k_(hour+1)) = repmat(w_res_hour(:,hour),1,3600/h +1);
    end
    
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

    w_load_ = w(load__mask,:);
    w_res = w(res_mask,:);
    

    t = 0:h:(size(w,2)-1)*h;
    % figure
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % hold on
    % grid on
    % box on;
    % stairs(t,w(load__mask,:)','LineWidth',1.5);
    % legend('$\Delta P_{{load_}_1}$','$\Delta P_{{load_}_2}$','$\Delta P_{{load_}_3}$','Interpreter','latex')
    % ylabel('$\Delta P_{load_}$ (pu)','interpreter','latex');
    % xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    % xticks(0:3600:simulation_hours*3600)
    % xticklabels(sprintfc('%d', 0:simulation_hours))
    % hold off
    % set(gcf,'resderer','Painters');
    % title='./fig/delta_p_load_.png';
    % saveas(gca,title,'png');
    % 
    % 
    % if(any(res_mask))
    %     figure
    %     set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %     hold on
    %     grid on
    %     box on;
    %     stairs(t,w(res_mask,:)','LineWidth',1.5);
    %     legend('$\Delta P_{{res}_1}$','$\Delta P_{{res}_2}$','$\Delta P_{{res}_3}$','Interpreter','latex')
    %     ylabel('$\Delta P_{res}$ (pu)','interpreter','latex');
    %     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    %     xticks(0:3600:simulation_hours*3600)
    %     xticklabels(sprintfc('%d', 0:simulation_hours))
    %     hold off
    %     set(gcf,'resderer','Painters');
    %     title='./fig/p_res.png';
    %     saveas(gca,title,'png');
    % end



    
    if(any(res_mask))
        data = load_('data\solar.mat');
        data = data.data;
        %data.data(:,2:end) = data.data(:,2:end)./100;
        if ceil(simulation_seconds/3600) ~= 1
            P_res = resample(data.data(1:ceil(simulation_seconds/3600),2:end),3600/h,1,'Dimension',1);
            if size(w_res,2) < size(P_res,1)
                P_res = P_res(1:size(w_res,2),:);
            else
                P_res(end+ size(w_res,2) - size(P_res,1),:) = P_res(end,:);
            end
        else
            P_res = data.data(1,2:end);
        end
        
        %Negative values might be introduced by the resampling
        mask = P_res < 0;
        P_res(mask) = 0;
    
        P_res = P_res';
        P_res = P_res + w_res;
        
    else
        P_res = [];
    end
    
    %w_load_ = zeros(size(w_load_));
    P_load_ = PL0+ w_load_;


    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,P_load_','LineWidth',1.5);
    %legend('$\Delta P_{{load_}_1}$','$\Delta P_{{load_}_2}$','$\Delta P_{{load_}_3}$','Interpreter','latex')
    ylabel('$P_{L}$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    xticks(0:3600:simulation_hours*3600)
    xticklabels(sprintfc('%d', 0:simulation_hours))
    hold off
    set(gcf,'resderer','Painters');
    title='./fig/p_load_.png';
    saveas(gca,title,'png');

    if(any(res_mask))
        figure
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        hold on
        grid on
        box on;
        stairs(t,P_res','LineWidth',1.5);
        %legend('$\Delta P_{{load_}_1}$','$\Delta P_{{load_}_2}$','$\Delta P_{{load_}_3}$','Interpreter','latex')
        ylabel('$P_{RES}$ (pu)','interpreter','latex');
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:simulation_hours*3600)
        xticklabels(sprintfc('%d', 0:simulation_hours))
        hold off
        set(gcf,'resderer','Painters');
        title='./fig/p_res.png';
        saveas(gca,title,'png');
    end

end

