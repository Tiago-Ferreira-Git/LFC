function [w,w_load,w_ren] = get_disturbance_profile(w,h,n_areas,simulation_seconds,bus_ss)
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
        w(load_mask,k) = w_load_hour(:,hour).*100;
        w(ren_mask,k) = w_ren_hour(:,hour);
        

        %w(load_mask,k_(hour):k_(hour+1)) = repmat(w_load_hour(:,hour),1,3600/h +1).*100;
        %w(ren_mask,k_(hour):k_(hour+1)) = repmat(w_ren_hour(:,hour),1,3600/h +1);
    end
    
    w_load = w(load_mask,:);
    w_ren = w(ren_mask,:);
    

    
    
    for i = 1:size(w,1)

        freq_pole = 0.00001;
        a = [1 -exp(-h*freq_pole)];
        b = 1-exp(-h*freq_pole);
        w(i,:) = filter(b,a,w(i,:));

    end

    t = 0:h:(size(w,2)-1)*h;
    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,w(load_mask,:)','LineWidth',1.5);
    legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
    ylabel('$\Delta P_{load}$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    xticks(0:3600:simulation_hours*3600)
    xticklabels(sprintfc('%d', 0:simulation_hours))
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/delta_p_load.png';
    saveas(gca,title,'png');
    
    
    if(any(ren_mask))
        figure
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        hold on
        grid on
        box on;
        stairs(t,w(ren_mask,:)','LineWidth',1.5);
        legend('$\Delta P_{{ren}_1}$','$\Delta P_{{ren}_2}$','$\Delta P_{{ren}_3}$','Interpreter','latex')
        ylabel('$\Delta P_{ren}$ (pu)','interpreter','latex');
        xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
        xticks(0:3600:simulation_hours*3600)
        xticklabels(sprintfc('%d', 0:simulation_hours))
        hold off
        set(gcf,'renderer','Painters');
        title='./fig/p_ren.png';
        saveas(gca,title,'png');
    end

end

