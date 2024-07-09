function dxdt = nonlinear_model(t,x,K,network,bus_ss,x0,u0,PL,Pres,u_index,C_mech,delta_u)
   

    delta_x = x - x0;

    if nargin < 12
        
        delta_u = -K*delta_x;
    end



    k = find(t >= 0:3600:24*3600);
    k = k(1);

    dxdt = zeros(size(x0));

    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];
    
    
    
    u = u0 + delta_u;
    
    index_res = 1;
    for i = 1:length(network)
        freq_feedback = zeros(size(network(i).A,1),1);
        freq_feedback(3:3:end) = -network(i).tg_con(:,4);
        area_index = freq_index(i):freq_index(i+1)-1;
        x_ = x(area_index);

        x_mec_index = freq_index(i)+1:freq_index(i+1)-2-network(i).res;
        u_mec_index = u_index(i):u_index(i+1)-1;

        dxdt(x_mec_index) = network(i).A*delta_x(x_mec_index) + network(i).B*delta_u(u_mec_index);
        %+ freq_feedback.*(delta_x(area_index(1)))  ;
        
        
        P_mech = network(i).C_mech*x_;
        
        % Frequency dynamics 
        if network(i).res ~= 0
            P_res = network(i).C_res*x_;
            dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(i,k))./network(i).inertia;


            dxdt(freq_index(i+1)-n_res-1:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+network(i).res,k);
            
        else
            dxdt(freq_index(i)) = ( (P_mech./x_(1))  - PL(i,k))./network(i).inertia;
        end

        
        % for j = 1:size(network(i).to_bus,1)
        % 
        %     ang_diff = 377*(x_(end) - x(angle_index(network(i).to_bus(j,1))));
        % 
        %     dxdt(freq_index(i)) = dxdt(freq_index(i)) + sin(ang_diff)/(network(i).to_bus(j,5))/network(i).inertia;
        % end

        % Frequency error dynamics
        %dxdt(angle_index(i)) = -delta_x(freq_index(i));

        index_res = index_res + network(i).res + 1;
      
    end

    C_mech*(x-x0);
    mask = abs(dxdt) < 1e-8;
    dxdt(mask) =  0;
end

