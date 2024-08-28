function dxdt = nonlinear_model(t,x,x0,x_mec_index,area,areas,bus,ptie_dist,network,PL,P_res)
    
    dxdt = zeros(size(x0));
    
    
    delta_x = x - x0;
    
    
    k = find(t >= 0:3600:24*3600);
    k = k(1);
    
    
    
    freq_index = area.freq_index;
    angle_index = area.freq_index-1;
    angle_index = angle_index(2:end);
     %First step: compute u 
    
     %delta_u = -K*delta_x;
     u = zeros(size(area.B,2),1);
    
    x_mec_index(1:3:end) = max(0,x_mec_index(1:3:end));
            
    
    P_mech = area.C_mech*x;
    P_mech = max(P_mech,0);
    
    freq_feedback = zeros(size(area.A,1),length(areas));
    Ptie = zeros(length(areas),1);
    inertia = zeros(length(areas),1);
    for j = 1:length(areas)
        i = areas(j); 
        n_res = network(i).res;
    
        freq_feedback(freq_index(j)+3:3:freq_index(j+1)-2,j) = -network(i).tg_con(:,4) ;
        
        
    
        %Assume neighbours of neighbours are constant
        if j ~= 1
            mask = network(i).to_bus(:,1) == areas(1);
            network(i).to_bus(~mask,:) = [];
        end
        phi = atan2(network(i).to_bus(:,5),(network(i).to_bus(:,4)));
        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);
    
        V_area = bus(network(i).to_bus(:,2),1);
        V_neighbour = bus(network(i).to_bus(:,3),1);
        
        angle_bus = deg2rad(bus(network(i).to_bus(:,2),2)) - x(angle_index(j));
         
        [~,mask] = ismember(network(i).to_bus(:,1),areas');
        angle_nei = deg2rad(bus(network(i).to_bus(:,3),2)) - x(angle_index(mask));
    
        mask = network(i).to_bus(:,4) == 0;
        if any(mask)
            z_mod(mask) = z_mod(mask).*network(i).to_bus(mask,10);
            
        end
        Ptie(j) = sum((V_area.^2).*cos(phi)./z_mod + V_area.*V_neighbour.*sin(angle_bus - angle_nei + phi - pi/2)./z_mod) + ptie_dist(j);
        
        inertia(j) = network(i).inertia;
    end
    
    freq_index = freq_index(1:end-1);
     % Compute mechanical power dynamics
     dynamics = area.A*x + area.B*u + freq_feedback*delta_x(freq_index);
     dxdt(x_mec_index) = dynamics(x_mec_index); 
     %freq_feedback.*(delta_x(area_index(1)))   ;
    
    
    % Frequency dynamics 
    if n_res ~= 0
        P_res = network(i).C_res*x_;
        dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(k) - Ptie)./network(areas).inertia;
        
        error 'RENEWABLES NOT SUPPORTED YET'
    
        %dxdt(freq_index(i+1)-n_res:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+n_res-1,k);
        %index_res = index_res + n_res;
        
    else
        dxdt(freq_index) = ( (P_mech./x(freq_index))  - PL(:,k)  - Ptie)./inertia;
        
    end
end

