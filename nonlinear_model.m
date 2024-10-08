function dxdt = nonlinear_model(t,x,network,bus_ss,x0,u0,PL,Pres,delta_u,bus,delta_load,K)
   

    


    dxdt = zeros(size(x0));

    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];
    
    
    % if t>1500
    %     t
    %     %delta_u = ones(size(delta_u)).*(sum(delta_load)/size(delta_u,1));
    % end
    
    
    %x(angle_index) = mod(x(angle_index) + pi, 2*pi) - pi;
    delta_x = x - x0;


    delta_u = -K*delta_x;

    u = u0 + delta_u;
    
    %u = max(u,0);
    

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
        



        theta_shift = network(i).to_bus(:,11);
        phi = atan2(network(i).to_bus(:,5),(network(i).to_bus(:,4)));
        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);

        V_area = bus(network(i).to_bus(:,2),1);
        V_neighbour = bus(network(i).to_bus(:,3),1);
        
        angle_bus = deg2rad(bus(network(i).to_bus(:,2),2)) - x(angle_index(i));
        angle_bus_0 = deg2rad(bus(network(i).to_bus(:,2),2));
        angle_nei = deg2rad(bus(network(i).to_bus(:,3),2)) - x(angle_index(network(i).to_bus(:,1)));
         angle_nei_0 = deg2rad(bus(network(i).to_bus(:,3),2));

        mask = bitand(network(i).to_bus(:,4) == 0, network(i).to_bus(:,10) ~= 0);
        if any(mask)
            z_mod(mask) = z_mod(mask).*network(i).to_bus(mask,10);
            
        end


        
        % ptie_teste = network(i).to_bus(:,end)'*(delta_x(angle_index(i)) - delta_x(angle_index(network(i).to_bus(:,1))));
        % 
        % ptie_teste_2 = network(i).to_bus(:,end)'*sin(delta_x(angle_index(i)) - delta_x(angle_index(network(i).to_bus(:,1))));


        Ptie_0 = sum((V_area.^2).*cos(phi)./z_mod + V_area.*V_neighbour.*sin(angle_bus_0 - angle_nei_0 - theta_shift + phi - pi/2)./z_mod);
        
        Ptie = sum((V_area.^2).*cos(phi)./z_mod + V_area.*V_neighbour.*sin(angle_bus - angle_nei - theta_shift + phi - pi/2)./z_mod);

        ptie = Ptie - Ptie_0;



                % Frequency dynamics 
        if n_res ~= 0
            P_res = network(i).C_res*x_;

            dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(i) - Pt0(i))./network(i).inertia;


            dxdt(freq_index(i+1)-n_res:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+n_res-1);
            index_res = index_res + n_res;
            
        else
            dxdt(freq_index(i)) = ( (P_mech./x_(1))  - PL(i) - Ptie )./network(i).inertia;
        end

        % Frequency error dynamics
        

        dxdt(angle_index(i)) = -2*pi*60*delta_x(area_index(1));
        % 
        % mask = abs(dxdt) < 1e-7;
        % dxdt(mask) = 0;

      
    end
end

