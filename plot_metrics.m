function  plot_metrics(n_areas,n_machines,simulation_hours,frequency_error_cost,disp_cost_area,disp_cost_machine,time_settling_cost,R_,to_plot)


    covariances = zeros(n_areas,1);
    means = zeros(n_areas,1);
    for i =1:n_areas
    
        covariances(i) = cov(to_plot(2:end-3,i));
        means(i) = mean(to_plot(2:end-3,i));
    
        % title = sprintf('./fig/K_tie_%d_%d.png',i,simulation_hours);
        % figure
        % set(gca,'TickLabelInterpreter','latex') % Latex style axis
        % hold on
        % grid on
        % box on;
        % plot(1:simulation_hours-5,to_plot(2:end-4,i)-to_plot(2,i),'LineWidth',1.5);
        % ylabel('$T_{{tie}_{i}} - T_{{tie}_{i,0}}$ ','interpreter','latex');
        % xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
        % hold off
        % set(gcf,'renderer','Painters');
        % saveas(gca,title,'png');
    
    
    end



    %%
    figure; 
    hold on;
    grid on;
    box on;
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(1:n_areas,frequency_error_cost(1,:),1:n_areas,frequency_error_cost(2,:),'LineWidth',1.5)
    %%%%
    legend({'Decentralized','Centralized'},...
	    'Location','best','Interpreter','latex');
    ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    xlabel('$Area $','Interpreter','latex');
    hold off;
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/error_R%f.png',R_);
    saveas(gca,title,'png');
    %%
    
    figure
    hold on;
    grid on;
    box on;
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(1:n_areas,disp_cost_area(1,:),1:n_areas,disp_cost_area(2,:),'LineWidth',1.5)
    %%%%
    legend({'Decentralized','Centralized'},...
	    'Location','best','Interpreter','latex');
    ylabel('$\sum_{k=1}^{k_{end}} \Delta u_i(k) (\mathrm{pu})$','Interpreter','latex');
    xlabel('$Area $','Interpreter','latex');
    hold off;
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/u_area_R%f.png',R_);
    saveas(gca,title,'png');
    
    
    
    
    figure
    hold on;
    grid on;
    box on;
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(1:n_machines,disp_cost_machine(1,:),1:n_machines,disp_cost_machine(2,:),'LineWidth',1.5)
    %%%%
    legend({'Decentralized','Centralized'},...
	    'Location','best','Interpreter','latex');
    ylabel('$\sum_{k=1}^{k_{end}} \Delta u_i(k) (\mathrm{pu})$','Interpreter','latex');
    xlabel('$Machine $','Interpreter','latex');
    hold off;
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/u_machine_R%f.png',R_);
    saveas(gca,title,'png');
    
    
    
    %%
    figure; 
    hold on;
    grid on;
    box on;
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(1:simulation_hours,time_settling_cost(1,:),1:simulation_hours,time_settling_cost(2,:),'LineWidth',1.5)
    %%%%
    legend({'Decentralized','Centralized'},...
	    'Location','best','Interpreter','latex');
    ylabel('Settling time $(\mathrm{s})$','Interpreter','latex');
    xlabel('Hours $(\mathrm{h})$','Interpreter','latex');
    hold off;
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/time_settling_R%f.png',R_);
    saveas(gca,title,'png');

end

