function dxdt = nonlinear_model(t,x,K,network,bus_ss,x0,u0,PL,Pres,Pt0,u_index,delta_u,debug,bus)
   

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
        n_res = network(i).res;
        freq_feedback = zeros(size(network(i).A,1),1);
        freq_feedback(3:3:end) = -network(i).tg_con(:,4);
        area_index = freq_index(i):freq_index(i+1)-1;
        x_ = x(area_index);

        x_mec_index = freq_index(i)+1:freq_index(i+1)-2-n_res;
        u_mec_index = u_index(i):u_index(i+1)-1;


        x_mec_index(1:3:end) = max(0,x_mec_index(1:3:end));

        dxdt(x_mec_index) = network(i).A*x(x_mec_index) + network(i).B*u(u_mec_index)+ freq_feedback.*(delta_x(area_index(1)))   ;
        
        
        
        P_mech = network(i).C_mech*x_;
        P_mech = max(P_mech,0);


        phi = atan2(network(i).to_bus(:,5),(network(i).to_bus(:,4)));
        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);

        V_area = bus(network(i).to_bus(:,2),1);
        V_neighbour = bus(network(i).to_bus(:,3),1);
        
        angle_bus = deg2rad(bus(network(i).to_bus(:,2),2)) - x(angle_index(i));
        angle_nei = deg2rad(bus(network(i).to_bus(:,3),2)) - x(angle_index(network(i).to_bus(:,1)));

        mask = network(i).to_bus(:,4) == 0;
        if any(mask)
            z_mod(mask) = z_mod(mask).*network(i).to_bus(mask,10);
            
        end
        
        Ptie = sum((V_area.^2).*cos(phi)./z_mod + V_area.*V_neighbour.*sin(angle_bus - angle_nei + phi - pi/2)./z_mod);
        
        % Frequency dynamics 
        if n_res ~= 0
            P_res = network(i).C_res*x_;
            dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(i,k) - Ptie)./network(i).inertia;


            dxdt(freq_index(i+1)-n_res:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+n_res-1,k);
            index_res = index_res + n_res;
            
        else
            dxdt(freq_index(i)) = ( (P_mech./x_(1))  - PL(i,k)  - Ptie)./network(i).inertia;
            %- Pt0(i)
        end

       
     
        
        %Ploss2(i)  = sum((V_bus.^2).*cos(phi)./z_mod + V_bus.*V_nei.*sin(angle_bus - angle_nei + phi - pi/2)./z_mod);
        
        % for j = 1:size(network(i).to_bus,1)
        % 
        %     if debug == 1
        %         %ang_diff = 2*pi*60*(delta_x(angle_index(i)) - delta_x(angle_index(network(i).to_bus(j,1))));
        %         Ptie = (V_area.^2).*cos(phi)./z_mod + V_area.*V_neighbour.*sin(angle_bus - 2*pi*60*x(angle_index(i)) - (angle_nei - 2*pi*60*x(angle_index(network(i).to_bus(j,1)))) + phi - pi/2)./z_mod;
        %     else
        %         %ang_diff = ((angle0(network(i).to_bus(j,2)) - x(angle_index(i))) - (angle0(network(i).to_bus(j,3)) - x(angle_index(network(i).to_bus(j,1)))));
        % 
        %         Ptie = sum((V_area.^2).*cos(phi)./z_mod + V_area.*V_neighbour.*sin(angle_bus - x(angle_index(i)) - (angle_nei - x(angle_index(network(i).to_bus(j,1)))) + phi - pi/2)./z_mod);
        % 
        %     end
        %     %Ptie = Ptie/(network(i).to_bus(j,5)); 
        % end

     
        % Frequency error dynamics
        if debug == 1
            dxdt(angle_index(i)) = -delta_x(freq_index(i));
        else
            dxdt(angle_index(i)) = -2*pi*60*delta_x(freq_index(i));
        end

        % mask = dxdt < 1e-12;
        % dxdt(mask) = 0;
      
    end
end

