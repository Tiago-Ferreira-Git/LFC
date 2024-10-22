function [K,E_fs,A,B,W] = slow_ss(mpc,network,h,A_orig,Pmech0)
    

    %% Generate the ss for each area
    
    bus_ss = [];
    A_global = [];
    B_global = [];
    W_global = [];

    % for i=1:length(network)
    %     network(i).tg_con(:,4) = 0;
    % end
    
    for i=1:length(network)
        
        B = zeros(2,network(i).machines);
        sum_area = 0;
        for j = 1:size(network(i).tg_con,1)
            
            Ts = network(i).tg_con(j,6); Tc = network(i).tg_con(j,7);  T5 = network(i).tg_con(j,10);

            a = Ts + Tc + T5;

               
           sum_area = sum_area - network(i).tg_con(j,4)/(a);
           
           B(1,j) = 1/(network(i).inertia*a);

        end
    

        A = zeros(2,2);

        %- network(i).damping
        A(1,1) = (sum_area - network(i).damping)/network(i).inertia;

        %Add error integrator 
        A(end,1) = 2*pi*60;

        W = zeros(2,1);

        W(1,1) = -1/(network(i).inertia);

        
        
        A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
        B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
        W_global = [W_global zeros(size(W_global,1), size(W,2)); zeros(size(W,1), size(W_global,2)) W];
        bus_ss = [bus_ss; size(network(i).tg_con,1) ];
        
    end
    
    
    
    %%

    bus_sol = mpc.bus(:,9);
    bus_sol = deg2rad(bus_sol);


    %Add tie-lines
    E = zeros(size(B_global,2),size(A_global,1));
    E_fs = zeros(size(B_global,2),size(A_global,1));
    k_ties = zeros(1,length(network));

    tg_size = cumsum([ 1 ; bus_ss]);


    for i = 1:length(network)


        neighbours = network(i).to_bus;

        E(tg_size(i):tg_size(i+1)-1,i*2-1:i*2) = ones(network(i).machines,2);
        E_fs(tg_size(i):tg_size(i+1)-1,i*2-1:i*2) = ones(network(i).machines,2);
        T_i = 0;
        for j = 1:size(neighbours,1)
            
            T_ji = cos(bus_sol(neighbours(j,3)) - bus_sol(neighbours(j,2)))/(neighbours(j,5));

            neigbour_index = neighbours(j,1)*2-1:neighbours(j,1)*2;

            %Changin ptie for error integral
            A_global(2*i-1, neigbour_index(end)) =  A_global(2*i-1, neigbour_index(end)) +T_ji/network(i).inertia;
            
            
            
            E(tg_size(i):tg_size(i+1)-1,neigbour_index) = ones(network(i).machines,2);


            T_i =  T_i + T_ji;
        end
        
        A_global(2*i-1, i*2) = -T_i/network(i).inertia;
        k_ties(i) = T_i;
    end


    figure
    hold on
    title("Reduced vs Original Model poles")
    plot(real(eig(A_orig)),imag(eig(A_orig)),'x')
    plot(real(eig(A_global)),imag(eig(A_global)),'x')
    legend({'Original','Reduced'},'Location','best')
    hold off


    figure
    hold on
    title("Reduced vs Original Model poles zoomed")
    plot(real(eig(A_orig)),imag(eig(A_orig)),'x')
    plot(real(eig(A_global)),imag(eig(A_global)),'x')
    legend({'Original','Reduced'},'Location','best')
    xline(0)
    xlim([-0.205,0.003])
    hold off
   
    [A,B,W] = discrete_dynamics(A_global,B_global,W_global,h);

    q = zeros(1,size(A,1));

    %k_ties = A_global(1:2:end,2:2:end);

    q(1:2:end) = 10.*k_ties;
    q(2:2:end) = 0.1.*k_ties;
    % 
    q(1:2:end) = 10;
    q(2:2:end) = 100;

   


    R_ties = zeros(1,size(B,2));

    j = 1;

    for i = 1:size(network,2)

        network(i).machines*size(network(i).to_bus,1);
        R_ties(j:j-1+network(i).machines) = k_ties(i);
        j = j + network(i).machines;
    end


    
    
    % E = E_to;
    % E = ones(size(E));
    %diag(1e2.*R_)
    tic
    %K = dlqr(A,B,diag(q),diag(1e2.*R_ties));


    %K = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q,R,E,300e3,1e-8,A,B,r,Z,T);

    %diag(1e2.*R_ties)
    % Pmech0 = C_mech*x0+1e-1;
    % % R = diag(R_./Pmech0);
    [K,~,trace_records] = LQROneStepLTI(A,B,diag(q),1e5*eye(size(B,2)),E,NaN);
    figure
    plot(trace_records(trace_records>0))
    xlabel('Iterations')
    ylabel("trace(P_{inf})",'Interpreter','tex')


    % 
    % opts.verbose = true;
    % opts.maxOLIt = 10;
    % opts.W = 40;
    % [K_cent,S] = dlqr(A,B,1e2*diag(q),1*eye(size(B,2)));
    % [K,Pinf] = LQRFiniteHorizonLTI(A,B,1e2*diag(q),1*eye(size(B,2)),E,K_cent,S,opts);
    % trace(Pinf)




    % K = dlqr(A,B,diag(q),diag(1e2*R_));
    % 
    % savefig('./fig/trace_my_method.fig');
    % set(gcf,'renderer','Painters');
    % saveas(gca,'./fig/trace_my_method.png','png');
    % toc

end

