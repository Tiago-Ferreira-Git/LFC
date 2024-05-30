function w = get_disturbance_profile(w,h,n_areas,simulation_hours,areas)
    w_load_hour = load('data\w_load.mat');
    w_load_hour = w_load_hour.w_load(:,1:24);
    
    w_ren_hour = load('data\w_ren.mat');
    w_ren_hour = w_ren_hour.w_ren(:,1:24);


    n_ren = sum(areas);

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
        
        


    if 24 < simulation_hours 
        w_load_hour = repmat(w_load_hour,1,floor(simulation_hours/24)+1);
        w_ren_hour =  repmat(w_ren_hour,1,floor(simulation_hours/24)+1);
    end

    load_index = zeros(1,size(w_load_hour,1));
    ren_index = zeros(1,size(w_ren_hour,1));

    load_index(1) = 1;
    k = 1;
    for i = 1:n_areas-1
       load_index(i+1) = areas(i) + load_index(i)+1;
       if areas(i) ~= 0
            ren_index(k:k+areas(i)-1) = load_index(i)+1:load_index(i+1)-1;
            k = k + areas(i);
       end
    end
    
    if areas(n_areas) ~= 0
        ren_index(k:k+areas(n_areas)-1) = load_index(n_areas)+1:load_index(n_areas)+areas(n_areas);
        k = k + areas(i);
    end
   


    k_ = 1:1440:size(w,2);
    for hour=1:simulation_hours
        
        w(load_index,k_(hour):k_(hour+1)) = repmat(w_load_hour(:,hour),1,1441)*100;
        w(ren_index,k_(hour):k_(hour+1)) = repmat(w_ren_hour(:,hour),1,1441);
        
    end

    t = 0:h:3600*simulation_hours;
    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,w(load_index,:)','LineWidth',1.5);
    legend('$\Delta P_{{load}_1}$','$\Delta P_{{load}_2}$','$\Delta P_{{load}_3}$','Interpreter','latex')
    ylabel('$\Delta P_{load}$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/delta_p_load.png';
    saveas(gca,title,'png');
    
    
    
    figure
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    hold on
    grid on
    box on;
    stairs(t,w(ren_index,:)','LineWidth',1.5);
    legend('$\Delta P_{{ren}_1}$','$\Delta P_{{ren}_2}$','$\Delta P_{{ren}_3}$','Interpreter','latex')
    ylabel('$\Delta P_{ren}$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/p_ren.png';
    saveas(gca,title,'png');

end

