function dxdt = nonlinear_model(t,x,network,bus_ss,x0,u0,PL,Pres,Pt0,delta_u)
   

    delta_x = x - x0;


    dxdt = zeros(size(x0));

    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];
    
    
    
    u = u0 + delta_u;
    
    % u = max(u,0);

    index_res = 1;
    for i = 1:length(network)
        n_res = network(i).res;
        freq_feedback = zeros(size(network(i).A,1),1);
        freq_feedback(3:3:end) = -network(i).freq_feedback;
        area_index = freq_index(i):freq_index(i+1)-1;
        x_ = x(area_index);

        x_mec_index = network(i).mec_index;
        u_mec_index = network(i).u_index;

        dxdt(x_mec_index) = network(i).A*x(x_mec_index) + network(i).B*u(u_mec_index)+ freq_feedback.*(delta_x(area_index(1)))  ;
        
        
        P_mech = network(i).C_mech*x_;
        
        % Frequency dynamics 
        if n_res ~= 0
            P_res = network(i).C_res*x_;

            dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(i) - Pt0(i))./network(i).inertia;


            dxdt(freq_index(i+1)-n_res:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+n_res-1);
            index_res = index_res + n_res;
            
        else
            dxdt(freq_index(i)) = ( (P_mech./x_(1))  - PL(i) - Pt0(i) )./network(i).inertia;
        end


        for j = 1:size(network(i).to_bus,1)
                      
           ang_diff = (delta_x(angle_index(i)) - delta_x(angle_index(network(i).to_bus(j,1))));
            
           dxdt(freq_index(i)) = dxdt(freq_index(i)) - (sin(ang_diff)/(network(i).to_bus(j,5) * network(i).inertia));
        end

        % Frequency error dynamics
        dxdt(angle_index(i)) = 2*pi*60*delta_x(area_index(1));
    end

    dxdt;
end

