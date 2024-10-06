function [w,w_load,w_ren,P_load,P_res] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss,PL0)
    w_load_hour = load('data\w_load.mat');
    w_load_hour = w_load_hour.w_load(:,1:24);
    
    w_ren_hour = load('data\w_ren.mat');
    w_ren_hour = w_ren_hour.w_ren(:,1:24);

    n_ren = sum(bus_ss(:,4));
    %if the number of areas in the data available is lower than the areas 
    if size(w_load_hour,1) < n_areas 
        w_load_hour = [repmat(w_load_hour, floor(n_areas/size(w_load_hour,1)), 1); ...
                        w_load_hour(1:mod(n_areas, size(w_load_hour,1)),:)];
    else 
        w_load_hour = w_load_hour(1:n_areas,:);
    end

    if size(w_ren_hour,1) < n_ren 
        w_ren_hour = [repmat(w_ren_hour, floor(n_ren/size(w_ren_hour,1)), 1); ...
                        w_ren_hour(1:mod(n_ren, size(w_ren_hour,1)),:)];
    else 
        w_ren_hour = w_ren_hour(1:n_ren,:);
    end

    simulation_hours = ceil(simulation_seconds/3600);
    
    %if the hours in the data available is lower than the simulation hour
    if 24 < simulation_hours 
        w_load_hour = repmat(w_load_hour,1,floor(simulation_hours/24)+1);
        w_ren_hour =  repmat(w_ren_hour,1,floor(simulation_hours/24)+1);
    else            
        w_load_hour = w_load_hour(:,1:simulation_hours);
        w_ren_hour = w_ren_hour(:,1:simulation_hours);
    end
    

    load_mask = zeros(1,size(w,1));
    load_mask(1) = 1;
    load_mask(cumsum(1+bus_ss(1:end-1,4))+1) = 1;
    load_mask = logical(load_mask);

    ren_mask = ~load_mask;


    k_ = 1:3600/h:simulation_hours*3600/h;

    for k=1:simulation_seconds/h + 1
        
        if ~isempty(find(k == k_, 1))
            hour = find(k == k_);
        end
        w(load_mask,k) = w_load_hour(:,hour)*100;
        w(ren_mask,k) = w_ren_hour(:,hour);
        

        %w(load_mask,k_(hour):k_(hour+1)) = repmat(w_load_hour(:,hour),1,3600/h +1).*100;
        %w(ren_mask,k_(hour):k_(hour+1)) = repmat(w_ren_hour(:,hour),1,3600/h +1);
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

    w_load = w(load_mask,:);
    w_ren = w(ren_mask,:);
    

    t = 0:h:(size(w,2)-1)*h;
    % figure
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % hold on
    % grid on
    % box on;
    % stairs(t,w(load_mask,:)','LineWidth',1.5);
    % legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
    % ylabel('$\Delta P_{load}$ (pu)','interpreter','latex');
    % xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    % xticks(0:3600:simulation_hours*3600)
    % xticklabels(sprintfc('%d', 0:simulation_hours))
    % hold off
    % set(gcf,'renderer','Painters');
    % title='./fig/delta_p_load.png';
    % saveas(gca,title,'png');
    % 
    % 
    % if(any(ren_mask))
    %     figure
    %     set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %     hold on
    %     grid on
    %     box on;
    %     stairs(t,w(ren_mask,:)','LineWidth',1.5);
    %     legend('$\Delta P_{{ren}_1}$','$\Delta P_{{ren}_2}$','$\Delta P_{{ren}_3}$','Interpreter','latex')
    %     ylabel('$\Delta P_{ren}$ (pu)','interpreter','latex');
    %     xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    %     xticks(0:3600:simulation_hours*3600)
    %     xticklabels(sprintfc('%d', 0:simulation_hours))
    %     hold off
    %     set(gcf,'renderer','Painters');
    %     title='./fig/p_ren.png';
    %     saveas(gca,title,'png');
    % end



    
    if(any(ren_mask))
        data = load('data\solar.mat');
        data = data.data;
        %data.data(:,2:end) = data.data(:,2:end)./100;
        if ceil(simulation_seconds/3600) ~= 1
            P_res = resample(data.data(1:ceil(simulation_seconds/3600),2:end),3600/h,1,'Dimension',1);
            if size(w_ren,2) < size(P_res,1)
                P_res = P_res(1:size(w_ren,2),:);
            else
                P_res(end+ size(w_ren,2) - size(P_res,1),:) = P_res(end,:);
            end
        else
            P_res = data.data(1,2:end);
        end
        
        %Negative values might be introduced by the resampling
        mask = P_res < 0;
        P_res(mask) = 0;
    
        P_res = P_res';
        P_res = P_res + w_ren;
        
    else
        P_res = [];
    end
    
    %w_load = zeros(size(w_load));
    P_load = PL0+ w_load;


    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,P_load','LineWidth',1.5);
    %legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
    ylabel('$P_{L}$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    xticks(0:3600:simulation_hours*3600)
    xticklabels(sprintfc('%d', 0:simulation_hours))
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/p_load.png';
    saveas(gca,title,'png');

    if(any(ren_mask))
        figure
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        hold on
        grid on
        box on;
        stairs(t,P_res','LineWidth',1.5);
        %legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
        ylabel('$P_{RES}$ (pu)','interpreter','latex');
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:simulation_hours*3600)
        xticklabels(sprintfc('%d', 0:simulation_hours))
        hold off
        set(gcf,'renderer','Painters');
        title='./fig/p_res.png';
        saveas(gca,title,'png');
    end

end

