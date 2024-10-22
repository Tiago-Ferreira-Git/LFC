function  plot_metrics(n_areas,n_machines,simulation_seconds,network,t,u_d,x_d,y_d,u_c,x_c,y_c,h,freq_index,angle_index,q,R)

    simulation_hours = ceil(simulation_seconds/3600);

    %Obtain ACE
    ace_d = zeros(n_areas,length(t));
    for i =1:n_areas       
        ace_d(i,:) = y_d((i-1)*4+4,:) + network(i).bias_factor.*(y_d((i-1)*4+1,:)-1); 
    end

    ace_c = zeros(n_areas,length(t));
    for i =1:n_areas       
        ace_c(i,:) = y_c((i-1)*4+4,:) + network(i).bias_factor.*(y_c((i-1)*4+1,:)-1); 
    end


    %frequency_error_cost,disp_cost_area,disp_cost_machine,time_settling_cost,R_,to_plot

    %Time settling - defined as the time unde 1e-5 freq

    time_settling_c  = zeros(simulation_hours,n_areas);

    %Frequency limits in pu
    freq_lim = 2e-5;
    %time transient to ignore in seconds
    k_ = 1;
    for k = 0:3600:simulation_seconds-3600
        t_transient = k +  100;
        for i =1:n_areas
            time_settling_c(k_,i) = t_settling(t,t_transient,freq_lim,(y_c((i-1)*4+1,:)-1),h,k);
            
        end
        k_ = k_ + 1;
    end

    time_settling_d  = zeros(simulation_hours,n_areas);

    %time transient to ignore in seconds
    k_ = 1;
    for k = 0:3600:simulation_seconds-3600
        t_transient = k +  100;
        for i =1:n_areas
            time_settling_d(k_,i) = t_settling(t,t_transient,freq_lim,(y_d((i-1)*4+1,:)-1),h,k);        
        end
        k_ = k_ + 1;
    end

    figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(t,ace_c,'LineWidth',1.5)
    %%%%
    %legend({'Frequency','Angle','Control Action'},...
	    %'Location','best','Interpreter','latex');
    %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    ylabel('ACE ','Interpreter','latex');
    if simulation_hours > 1

        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    savefig('./fig/ace_c.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/ace_c.eps','epsc');


    figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(t,ace_d,'LineWidth',1.5)
    %%%%
    %legend({'Frequency','Angle','Control Action'},...
	    %'Location','best','Interpreter','latex');
    %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    ylabel('ACE ','Interpreter','latex');
    if simulation_hours > 1
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    savefig('./fig/ace_d.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/ace_d.eps','epsc');







    %Cost function

    q_freq = zeros(size(q));
    q_freq(freq_index) = q(freq_index);

    q_angle = zeros(size(q));
    q_angle(angle_index) = q(angle_index);

    cost_freq_d = zeros(size(t));
    cost_angle_d = zeros(size(t));
    cost_u_d = zeros(size(t));
    energy_d = zeros(size(t));


    
    cost_freq_c = zeros(size(t));
    cost_angle_c = zeros(size(t));
    cost_u_c = zeros(size(t));
    energy_c = zeros(size(t));

    for i = 1:length(t)
        cost_freq_d(i) = x_d(:,i)'*diag(q_freq)*x_d(:,i);
        cost_angle_d(i) = x_d(:,i)'*diag(q_angle)*x_d(:,i);
        cost_u_d(i) = u_d(:,i)'*R*u_d(:,i);
        energy_d(i) = u_d(:,i)'*u_d(:,i);
    end

    for i = 1:length(t)
        cost_freq_c(i) = x_c(:,i)'*diag(q_freq)*x_c(:,i);
        cost_angle_c(i) = x_c(:,i)'*diag(q_angle)*x_c(:,i);
        cost_u_c(i) = u_c(:,i)'*R*u_c(:,i);
        energy_c(i) = u_c(:,i)'*u_c(:,i);
    end


    



    figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(t,cost_freq_d','LineWidth',1.5)
    plot(t,cost_angle_d','LineWidth',1.5)
    plot(t,cost_u_d','LineWidth',1.5)

    %%%%
    legend({'Frequency','Angle','Control Action'},...
	    'Location','best','Interpreter','latex');
    %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    ylabel('Cost Function ','Interpreter','latex');
    
    if simulation_hours > 1
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    savefig('./fig/J_d.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/J_d.eps','epsc');
    
    figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    plot(t,cost_freq_c','LineWidth',1.5)
    plot(t,cost_angle_c','LineWidth',1.5)
    plot(t,cost_u_c','LineWidth',1.5)

    %%%%
    legend({'Frequency','Angle','Control Action'},...
	    'Location','best','Interpreter','latex');
    %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
    ylabel('Cost Function ','Interpreter','latex');
    if simulation_hours > 1
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    savefig('./fig/J_c.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/J_c.eps','epsc');



    figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    
   
    plot(t,energy_c','LineWidth',1.5)
    plot(t,energy_d','LineWidth',1.5)

    %%%%
    legend({'Centralized','Decentralized'},...
	    'Location','best','Interpreter','latex');
    
    ylabel('Energy - $|\Delta u_\mathrm{i} |^2$','Interpreter','latex');
    %ylabel('Cost Function ','Interpreter','latex');
    if simulation_hours > 1
        xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
    else
        xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
    end
    savefig('./fig/energy.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/energy.eps','epsc');



    figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    %%%%
    
   
    plot(mean(time_settling_c,2),'LineWidth',1.5)
    plot(mean(time_settling_d,2),'LineWidth',1.5)

    %%%%
    legend({'Centralized','Decentralized'},...
	    'Location','best','Interpreter','latex');
    
    ylabel('Time Settling $\;[\mathrm{s }]$','Interpreter','latex');
    %ylabel('Cost Function ','Interpreter','latex');
    xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    xticks(0:1:simulation_hours)
    savefig('./fig/time_settling.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/time_settling.eps','epsc');





end

