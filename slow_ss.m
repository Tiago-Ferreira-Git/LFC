function [K,E_fs] = slow_ss(mpc,debug,network,h)
    

    %% Generate the ss for each area
    
    bus_ss = [];
    A_global = [];
    B_global = [];

    sum_area = 0;
    for i=1:length(network)
        
        B = zeros(2,network(i).machines);

        for j = 1:size(network(i).tg_con,1)
            
            Ts = network(i).tg_con(j,6); Tc = network(i).tg_con(j,7);  T5 = network(i).tg_con(j,10);

            a = Ts + Tc + T5;

               
           sum_area = sum_area - network(i).tg_con(j,4)/(a);
           
           B(1,j) = 1/(network(i).inertia*a);

        end
    

        A = zeros(2,2);
        A(1,1) = (sum_area + network(i).damping)/network(i).inertia;

        %Add error integrator 
        A(end,1) = -2*pi*60;
        
        
        A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
        B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
        bus_ss = [bus_ss; size(network(i).tg_con,1) ];
        
    end
    
    
    
    %%

    bus_sol = mpc.bus(:,9);
    bus_sol = deg2rad(bus_sol);


    %Add tie-lines
    E = zeros(size(B_global,2),size(A_global,1));
    E_fs = zeros(size(B_global,2),size(A_global,1));


    tg_size = cumsum([ 1 ; bus_ss]);


    for i = 1:size(network,2)


        neighbours = network(i).to_bus;

        E(tg_size(i):tg_size(i+1)-1,i*2-1:i*2) = ones(network(i).machines,2);
        E_fs(tg_size(i):tg_size(i+1)-1,i*2-1:i*2) = ones(network(i).machines,2);
        T_i = 0;
        for j = 1:size(neighbours,1)
            
            if debug == 1
                T_ji = 2*pi*60*cos(bus_sol(neighbours(j,3)) - bus_sol(neighbours(j,2)))/(neighbours(j,5)); 
            else
                T_ji = cos(bus_sol(neighbours(j,3)) - bus_sol(neighbours(j,2)))/(neighbours(j,5));
            end

            neigbour_index = neighbours(j,1)*2-1:neighbours(j,1)*2;

            %Changin ptie for error integral
            A_global(2*i-1, neigbour_index(end)) =  A_global(2*i-1, neigbour_index(end)) -T_ji/network(i).inertia;
            
            
            
            E(tg_size(i):tg_size(i+1)-1,neigbour_index) = ones(network(i).machines,2);


            T_i =  T_i + T_ji;
        end
        
        A_global(2*i-1, i*2) = T_i/network(i).inertia;

    end

   
    [A,B,~] = discrete_dynamics(A_global,B_global,zeros(size(A_global,1),1),h);

    q = zeros(1,size(A,1));
    q(1:2:end) = 0;
    q(2:2:end) = 1000;


    
    % E_to = zeros(size(E));
    % 
    % 
    % for i = 1:size(network,2)
    % 
    % 
    %     neighbours = network(i).to_bus;
    % 
    %     E_to(tg_size(i):tg_size(i+1)-1,i*2-1:i*2) = ones(network(i).machines,2);
    % 
    % 
    %     unique_areas = unique(network(i).to_bus(:,1),'rows');
    % 
    % 
    %     for j = unique_areas'
    %         neighbours = [neighbours ; network(j).to_bus];
    %     end
    % 
    %     % Third order neighbours
    %     unique_areas = unique(neighbours(:,1),'rows');
    % 
    %     for j = unique_areas'
    %         neighbours = [neighbours ; network(j).to_bus];
    %     end
    % 
    % 
    %     for j = 1:size(neighbours,1)
    %         neigbour_index = neighbours(j,1)*2-1:neighbours(j,1)*2;
    %         E_to(tg_size(i):tg_size(i+1)-1,neigbour_index) = ones(network(i).machines,2);
    %     end
    % 
    % end
    % 
    % 
    % E = E_to;

    tic
    [K,~,trace_records] = LQROneStepLTI(A,B,diag(q),0.1*eye(size(B,2)),E,NaN);
    figure
    plot(trace_records(trace_records>0))
    xlabel('Iterations')
    ylabel("trace(P_{inf})",'Interpreter','tex')
    
    savefig('./fig/trace_my_method.fig');
    set(gcf,'renderer','Painters');
    saveas(gca,'./fig/trace_my_method.png','png');
    toc

end

