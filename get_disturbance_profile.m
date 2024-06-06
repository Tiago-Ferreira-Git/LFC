function w = get_disturbance_profile(w,h,n_areas,simulation_hours,bus_ss)
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

    
    
    %if the hours in the data available is lower than the simulation hour
    if 24 < simulation_hours 
        w_load_hour = repmat(w_load_hour,1,floor(simulation_hours/24)+1);
        w_ren_hour =  repmat(w_ren_hour,1,floor(simulation_hours/24)+1);
    end
    
    load_mask = zeros(1,size(w,1));
    load_mask(1) = 1;
    load_mask(cumsum(1+bus_ss(1:end-1,4))+1) = 1;
    load_mask = logical(load_mask);

    ren_mask = ~load_mask;

    k_ = 1:3600/h:size(w,2);
    for hour=1:simulation_hours
        w(load_mask,k_(hour):k_(hour+1)) = repmat(w_load_hour(:,hour),1,3600/h +1 )*100;
        w(ren_mask,k_(hour):k_(hour+1)) = repmat(w_ren_hour(:,hour),1,3600/h +1);
    end


    t = 0:h:3600*simulation_hours;
    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,w(load_mask,:)','LineWidth',1.5);
    legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
    ylabel('$\Delta P_{load}$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
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
        stairs(t,w(2:2:end,:)','LineWidth',1.5);
        legend('$\Delta P_{{ren}_1}$','$\Delta P_{{ren}_2}$','$\Delta P_{{ren}_3}$','Interpreter','latex')
        ylabel('$\Delta P_{ren}$ (pu)','interpreter','latex');
        xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
        hold off
        set(gcf,'renderer','Painters');
        title='./fig/p_ren.png';
        saveas(gca,title,'png');
    end

end

