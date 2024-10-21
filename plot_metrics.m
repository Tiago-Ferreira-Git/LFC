function  plot_metrics(n_areas,n_machines,simulation_hours,network,t,u,x,y,h,freq_index,angle_index,q,R)


    %Obtain ACE
    ace = zeros(n_areas,length(t));
    for i =1:n_areas       
        ace(i,:) = y((i-1)*4+4,:) + network(i).bias_factor.*(y((i-1)*4+1,:)-1); 
    end


    %frequency_error_cost,disp_cost_area,disp_cost_machine,time_settling_cost,R_,to_plot

    %Time settling - defined as the time unde 1e-5 freq

    % time_settling  = zeros(1,n_areas);
    % 
    % %Frequency limits in pu
    % freq_lim = 2e-6;
    % %time transient to ignore in seconds
    % t_transient = 100;
    % for i =1:n_areas
    %     %find the first time where delta_freq was under 1e-5pu
    % 
    %     index = find(abs((y((i-1)*4+1,:)-1)) < freq_lim);
    % 
    %     mask = find(t >= t_transient);
    %     mask = mask(1);
    %     %Ignore initial transient
    %     mask = index > mask;
    %     index(~mask) = [];
    % 
    %     for j = 1:length(index)- round(5/h)
    %         if(index(j+round(5/h)) - index(j) == round(5/h))
    %             index = index(j);
    %             break
    %         end
    %     end
    % 
    %     if length(index) ~=1
    %         error "Could not find time settling"
    %     end
    % 
    %     time_settling(i) = t(index);
    % end


    figure; 
    hold on;
    grid on;
    box on
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(t,ace','LineWidth',1.5)
    %%%%
    %legend({'Frequency','Angle','Control Action'},...
	    %'Location','best','Interpreter','latex');
    %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    ylabel('ACE ','Interpreter','latex');
    if simulation_hours > 1
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:simulation_hours*3600)
        xticklabels(sprintfc('%d', 0:simulation_hours))
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/ACE.png');
    saveas(gca,title,'png');






    %Cost function

    q_freq = zeros(size(q));
    q_freq(freq_index) = q(freq_index);

    q_angle = zeros(size(q));
    q_angle(angle_index) = q(angle_index);

    cost_freq = zeros(size(t));
    cost_angle = zeros(size(t));
    cost_u = zeros(size(t));

    for i = 1:length(t)
        cost_freq(i) = x(:,i)'*diag(q_freq)*x(:,i);
        cost_angle(i) = x(:,i)'*diag(q_angle)*x(:,i);
        cost_u(i) = u(:,i)'*R*u(:,i);
    end
    
    



    figure; 
    hold on;
    grid on;
    box on
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(t,cost_freq','LineWidth',1.5)
    plot(t,cost_angle','LineWidth',1.5)
    plot(t,cost_u','LineWidth',1.5)

    %%%%
    legend({'Frequency','Angle','Control Action'},...
	    'Location','best','Interpreter','latex');
    %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    ylabel('Cost Function ','Interpreter','latex');
    if simulation_hours > 1
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:simulation_hours*3600)
        xticklabels(sprintfc('%d', 0:simulation_hours))
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    set(gcf,'renderer','Painters');
    title=sprintf('./fig/J.png');
    saveas(gca,title,'png');


    % 
    % %%
    % figure; 
    % hold on;
    % grid on;
    % box on;
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % %%%%
    % plot(1:n_areas,frequency_error_cost(1,:),1:n_areas,frequency_error_cost(2,:),'LineWidth',1.5)
    % %%%%
    % legend({'Decentralized','Centralized'},...
	%     'Location','best','Interpreter','latex');
    % ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    % xlabel('$Area $','Interpreter','latex');
    % hold off;
    % set(gcf,'renderer','Painters');
    % title=sprintf('./fig/error_R%f.png',R_);
    % saveas(gca,title,'png');
    % %%
    % 
    % figure
    % hold on;
    % grid on;
    % box on;
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % %%%%
    % plot(1:n_areas,disp_cost_area(1,:),1:n_areas,disp_cost_area(2,:),'LineWidth',1.5)
    % %%%%
    % legend({'Decentralized','Centralized'},...
	%     'Location','best','Interpreter','latex');
    % ylabel('$\sum_{k=1}^{k_{end}} \Delta u_i(k) (\mathrm{pu})$','Interpreter','latex');
    % xlabel('$Area $','Interpreter','latex');
    % hold off;
    % set(gcf,'renderer','Painters');
    % title=sprintf('./fig/u_area_R%f.png',R_);
    % saveas(gca,title,'png');
    % 
    % 
    % 
    % 
    % figure
    % hold on;
    % grid on;
    % box on;
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % %%%%
    % plot(1:n_machines,disp_cost_machine(1,:),1:n_machines,disp_cost_machine(2,:),'LineWidth',1.5)
    % %%%%
    % legend({'Decentralized','Centralized'},...
	%     'Location','best','Interpreter','latex');
    % ylabel('$\sum_{k=1}^{k_{end}} \Delta u_i(k) (\mathrm{pu})$','Interpreter','latex');
    % xlabel('$Machine $','Interpreter','latex');
    % hold off;
    % set(gcf,'renderer','Painters');
    % title=sprintf('./fig/u_machine_R%f.png',R_);
    % saveas(gca,title,'png');
    % 
    % 
    % 
    % %%
    % figure; 
    % hold on;
    % grid on;
    % box on;
    % set(gca,'TickLabelInterpreter','latex') % Latex style axis
    % %%%%
    % plot(1:simulation_hours,time_settling_cost(1,:),1:simulation_hours,time_settling_cost(2,:),'LineWidth',1.5)
    % %%%%
    % legend({'Decentralized','Centralized'},...
	%     'Location','best','Interpreter','latex');
    % ylabel('Settling time $(\mathrm{s})$','Interpreter','latex');
    % xlabel('Hours $(\mathrm{h})$','Interpreter','latex');
    % hold off;
    % set(gcf,'renderer','Painters');
    % title=sprintf('./fig/time_settling_R%f.png',R_);
    % saveas(gca,title,'png');

end

