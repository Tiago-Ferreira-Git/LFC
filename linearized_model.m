function dxdt = linearized_model(t,x,network,bus_ss,x0,u0,delta_PL,Pres,delta_u,bus)
    
    
    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];

    % x(angle_index) = mod(x(angle_index) + pi, 2*pi) - pi;
    delta_x = x - x0;


    dxdt = zeros(size(x0));

    
    
    u = u0 + delta_u;
    
    u = max(u,0);
    

    index_res = 1;
    for i = 1:length(network)
        n_res = network(i).res;
        freq_feedback = zeros(size(network(i).A,1),1);
        freq_feedback(3:3:end) = -network(i).freq_feedback;
        area_index = freq_index(i):freq_index(i+1)-1;
        x_ = x(area_index);

        x_mec_index = network(i).mec_index;
        u_mec_index = network(i).u_index;

        dxdt(x_mec_index) = network(i).A*x(x_mec_index) + network(i).B*u(u_mec_index)+ freq_feedback.*(delta_x(area_index(1)));
       

        delta_P_mech = network(i).C_mech*delta_x(area_index);

        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);
        mask = network(i).to_bus(:,4) == 0;
        if any(mask)
            z_mod(mask) = z_mod(mask).*network(i).to_bus(mask,10);
        end

        theta_shift = network(i).to_bus(:,11);
        V_area = bus(network(i).to_bus(:,2),1);
        V_neighbour = bus(network(i).to_bus(:,3),1);

        phi = atan2(network(i).to_bus(:,5),(network(i).to_bus(:,4)));

        Tij = V_area.*V_neighbour.*cos(deg2rad(bus(network(i).to_bus(:,2),2)) - deg2rad(bus(network(i).to_bus(:,3),2)) - theta_shift + phi-pi/2)./z_mod;


        angle_bus = x(angle_index(i));
        angle_nei = x(angle_index(network(i).to_bus(:,1)));
               
        delta_ptie = Tij'*(angle_bus - angle_nei);


        % Frequency dynamics 
        if n_res ~= 0
            P_res = network(i).C_res*x_;

            dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(i) - Pt0(i))./network(i).inertia;


            dxdt(freq_index(i+1)-n_res:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+n_res-1);
            index_res = index_res + n_res;
            
        else
            dxdt(freq_index(i)) = ( -network(i).damping*delta_x(area_index(1)) + delta_P_mech  - delta_PL(i) + delta_ptie )./network(i).inertia;
        end

        % Frequency error dynamics
        dxdt(angle_index(i)) = -2*pi*60*delta_x(area_index(1));


      
    end
end

