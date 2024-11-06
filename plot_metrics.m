function  [J_c, J_d] = plot_metrics(n_areas,n_machines,simulation_seconds,sampling_measurements,network,t,u0,u_d,x_d,y_d,u_c,x_c,y_c,h,freq_index,angle_index,q,R,flag_plot_metrics)
    

    simulation_hours = ceil(simulation_seconds/3600);
   
    if flag_plot_metrics
        %Obtain ACE
        ace = zeros(n_areas,length(t));
        for i =1:n_areas       
            ace(i,:) = ace(i,:) +  abs(y_c(i*4,:) + network(i).bias_factor.*(y_c((i-1)*4+1,:)-1)) - abs(y_d(i*4,:) + network(i).bias_factor.*(y_d((i-1)*4+1,:)-1)); 
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
        plot(t(1:sampling_measurements:end),ace(:,1:sampling_measurements:end),'LineWidth',1.5)
        %%%%
        %ylim([-0.2 0.4])
        legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','best','Interpreter','latex');
        %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
        ylabel('$|\mathbf{ACE}_\mathrm{c}| - |\mathbf{ACE}_\mathrm{d}|$ ','Interpreter','latex');
        if simulation_hours > 1
    
            xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
            xticks(0:3600:(simulation_hours-1)*3600)
            xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
            xlim([0 (simulation_hours-1)*3600])
        else
            xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
        end
        savefig('./fig/ace.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/ace.eps','epsc');
        saveas(gca,'./fig/ace.png','png');
    
    
    
    
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
        J_c = cost_freq_c + cost_angle_c + cost_u_c;
        J_d = cost_freq_d + cost_angle_d + cost_u_d;
    
        %%%%
        plot(t(1:sampling_measurements:end),J_c(1:sampling_measurements:end)','LineWidth',1.5)
        plot(t(1:sampling_measurements:end),J_d(1:sampling_measurements:end)','LineWidth',1.5)
    
        %%%%
        legend({'Centralized','Decentralized'},...
	        'Location','best','Interpreter','latex');
        %ylabel('Frequency error ($|0 - \Delta \omega |$) $(\mathrm{pu})$','Interpreter','latex');
        ylabel('Cost Function ','Interpreter','latex');
    
        if simulation_hours > 1
            xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
            xticks(0:3600:(simulation_hours-1)*3600)
            xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
            xlim([0 (simulation_hours-1)*3600])
        else
            xlabel('$t \;[\mathrm{s }]$','Interpreter','latex');
        end
        savefig('./fig/J.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/J.eps','epsc');
        saveas(gca,'./fig/J.png','png');
    
        J_c = sum(J_c);
        J_d = sum(J_d);
    
        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        %%%%
    
    
        plot(t(1:sampling_measurements:end),energy_c(1:sampling_measurements:end)'-energy_d(1:sampling_measurements:end)','LineWidth',1.5)
        %plot(t,','LineWidth',1.5)
    
        %%%%
        %legend({'Centralized','Decentralized'},...
	    %    'Location','best','Interpreter','latex');
    
        ylabel('$\mathbf{E}_\mathrm{c} - \mathbf{E}_\mathrm{d}$','Interpreter','latex');
        %ylabel('Cost Function ','Interpreter','latex');
        if simulation_hours > 1
            xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
            xticks(0:3600:(simulation_hours-1)*3600)
            xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
            xlim([0 (simulation_hours-1)*3600])
        else
            xlabel('$t \;[\mathrm{s}]$','Interpreter','latex');
        end
        savefig('./fig/energy.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/energy.eps','epsc');
        saveas(gca,'./fig/energy.png','png');
    
    
    
        figure('Position',4*[0 0 192 144]); % Nice aspect ratio for double column
        hold on;
        grid on;
        box on;
        set(gca,'FontSize',20);
        set(gca,'TickLabelInterpreter','latex') % Latex style axis
        scatter(0:simulation_hours-1,mean(time_settling_c,2),"filled")
        scatter(0:simulation_hours-1,mean(time_settling_d,2),"filled")
        legend({'Centralized','Decentralized'},...
	        'Location','best','Interpreter','latex');
        xlim([-0.5 23.5])
        ylabel('Time Settling $\;(\mathrm{s})$','Interpreter','latex');
        %ylabel('Cost Function ','Interpreter','latex');
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
        %xticks(-1:1:simulation_hours)
        savefig('./fig/time_settling.fig');
        set(gcf,'renderer','Painters');
        saveas(gca,'./fig/time_settling.eps','epsc');
        saveas(gca,'./fig/time_settling.png','png');
    end



    t_zoom = 20;
    t_zoom_2 = 2;
    figure('Position',4*[0 0 192*1.5 144*1.1]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    title('Centralized')
    stairs(t(1:sampling_measurements:end-1),y_c(1:4:end,1:sampling_measurements:end-1)','LineWidth',1.5);
    ylim([0.98 1.015])
    yline(1+4e-3,'--','LineWidth',1.5)
    yline(1-4e-3,'--','LineWidth',1.5)
    ylabel('$\omega$ (pu)','interpreter','latex');
    xlabel('$t \;(\mathrm{s})$','Interpreter','latex');
    if simulation_hours > 1
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');

        axes('position',[.55 .175 .25 .25])
        box on % put box around new pair of axes
        mask = (t < (t_zoom+1)*3600) & (t > -100 + (t_zoom)*3600); % range of t near perturbation
        stairs(t(mask),y_c(1:4:end,mask)','LineWidth',1.5) % plot on new axes
        %stairs(t_second(1:end-1),y_c(1:4:end,1:end-1)','LineWidth',1.5);
        xticks((t_zoom)*3600:3600:(t_zoom+1)*3600)
        xticklabels(sprintfc('%d', t_zoom:t_zoom+1))
        axis padded
        
        axes('position',[.275 .655 .25 .25])
        box on % put box around new pair of axes
        mask = (t < (t_zoom_2+5)*3600) & (t > -100 + (t_zoom_2)*3600); % range of t near perturbation
        stairs(t(mask),y_c(1:4:end,mask)','LineWidth',1.5) % plot on new axes
        %stairs(t_second(1:end-1),y_c(1:4:end,1:end-1)','LineWidth',1.5);
        xticks((t_zoom_2)*3600:3600:(t_zoom_2+6)*3600)
        xticklabels(sprintfc('%d', t_zoom_2:t_zoom_2+6))
        axis padded
        hold off
    else
        J_c = 0;
        J_d = 0;
    end
    
    savefig('./fig/frequency_c.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/frequency_c.eps','epsc')
    saveas(gca,'./fig/frequency_c.png','png')
    
    
    figure('Position',4*[0 0 192*1.5 144*1.1]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    title('Decentralized')
    stairs(t(1:sampling_measurements:end-1),y_d(1:4:end,1:sampling_measurements:end-1)','LineWidth',1.5);
    ylim([0.98 1.015])
    yline(1+4e-3,'--','LineWidth',1.5)
    yline(1-4e-3,'--','LineWidth',1.5)
    ylabel('$\omega$ (pu)','interpreter','latex');
    xlabel('$t \;(\mathrm{s})$','Interpreter','latex');
    if simulation_hours > 1
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
        axes('position',[.55 .175 .25 .25])
        box on % put box around new pair of axes
        mask = (t < (t_zoom+1)*3600) & (t > -100 + (t_zoom)*3600); % range of t near perturbation
        stairs(t(mask),y_d(1:4:end,mask)','LineWidth',1.5) % plot on new axes
        xticks((t_zoom)*3600:3600:(t_zoom+1)*3600)
        xticklabels(sprintfc('%d', t_zoom:t_zoom+1))
        axis padded
        yline(1+4e-3)
        yline(1+4e-3)
        axes('position',[.275 .655 .25 .25])
        box on % put box around new pair of axes
        mask = (t < (t_zoom_2+5)*3600) & (t > -100 + (t_zoom_2)*3600); % range of t near perturbation
        stairs(t(mask),y_d(1:4:end,mask)','LineWidth',1.5) % plot on new axes
        xticks((t_zoom_2)*3600:3600:(t_zoom_2+6)*3600)
        xticklabels(sprintfc('%d', t_zoom_2:t_zoom_2+6))
        axis padded
        hold off
    end
    
    
    savefig('./fig/frequency_d.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/frequency_d.eps','epsc')
    saveas(gca,'./fig/frequency_d.png','png')
    
    
    
    figure('Position',4*[0 0 192*1.5 144*1.1]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    title('Decentralized')
    stairs(t(1:sampling_measurements:end),(u0(:,1:sampling_measurements:end) + u_d(:,1:sampling_measurements:end))','LineWidth',1.5);
    %pause(5)
    legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','west','Interpreter','latex');
    xlabel('$t \;(\mathrm{s})$','Interpreter','latex');
    if simulation_hours > 1
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
    end
    ylabel('$\mathbf{u}$ (pu)','interpreter','latex');
    savefig('./fig/u_d.fig');
    set(gcf,'renderer','Painters');
    pause(5)
    saveas(gca,'./fig/u_d.eps','epsc');
    saveas(gca,'./fig/u_d.png','png');
    hold off
    
    
    
    figure('Position',4*[0 0 192*1.5 144*1.1]); % Nice aspect ratio for double column
    hold on;
    grid on;
    box on;
    set(gca,'FontSize',20);
    set(gca,'TickLabelInterpreter','latex') % Latex style axis
    title('Centralized')
    stairs(t(1:sampling_measurements:end),(u0(:,1:sampling_measurements:end) + u_c(:,1:sampling_measurements:end))','LineWidth',1.5);
    %pause(5)
    legend({'Area 1','Area 2','Area 3','Area 4','Area 5'},'Location','west','Interpreter','latex');
    xlabel('$t \;(\mathrm{s})$','Interpreter','latex');
    if simulation_hours > 1
        xticks(0:3600:(simulation_hours-1)*3600)
        xticklabels(sprintfc('%d', 0:(simulation_hours-1)))
        xlim([0 (simulation_hours-1)*3600])
        xlabel('$t \;(\mathrm{h})$','Interpreter','latex');
    end
    ylabel('$\mathbf{u}$ (pu)','interpreter','latex');
    savefig('./fig/u_c.fig');
    set(gcf,'renderer','Painters');
    pause(5)
    saveas(gca,'./fig/u_c.eps','epsc');
    saveas(gca,'./fig/u_c.png','png');
    hold off
 




end

