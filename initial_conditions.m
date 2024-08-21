function [x0,u0,Pt0,PL0,Ploss] = initial_conditions(ss_dim,n_machines,bus_ss,network,mpc)

    %Pt0 - Tie-line initial value per area
    %PL0 - Initial Consuption values per area


    x0 = zeros(ss_dim,1);
    u0 = zeros(n_machines,1);
    PL0 = zeros(length(network),1);
    Pt0 = zeros(length(network),1);
    Ploss = zeros(length(network),1);

    angle_index = cumsum(bus_ss);
    freq_index = [1 ; angle_index+1];
    
    x0(freq_index(1:end-1),1) = 1;
   
    j = 1;
    %%
    for i = 1:length(network)
        
        gen_index = freq_index(i)+1:freq_index(i+1)-2;
    
        x0_gen = zeros(network(i).machines*3+network(i).res,1);
        
    
        x0_gen(1:3:3*(network(i).machines)) = repelem(mpc.gen(network(i).mac_nr,2)./mpc.baseMVA,1,1);
        n = 3*network(i).machines;
        if network(i).res ~= 0
            x0_gen(n+1:1: n+ network(i).res) = mpc.gen(network(i).res_nr,2)./mpc.baseMVA;
        end
        x0(gen_index) = x0_gen;
        x0(angle_index(i)) =  deg2rad(mean(mpc.bus(network(i).mac_nr,9)));
        u0(j:j+network(i).machines-1) = mpc.gen(network(i).mac_nr,2)./mpc.baseMVA;
        j = j+network(i).machines;

        x0(angle_index(i)) =  0;

        PL0(i) = sum(mpc.bus(network(i).bus,3)./100);

        
        
        mask_outside_area = bitxor(ismember(mpc.branch(:,1),network(i).bus),ismember(mpc.branch(:,2),network(i).bus));
        mask = ismember(mpc.branch(:,1),network(i).bus,'rows');
        mask = bitand(mask_outside_area,mask);
        Pt0(i) = sum(mpc.branch(mask,14))./mpc.baseMVA;

        mask = ismember(mpc.branch(:,2),network(i).bus,'rows');
        mask = bitand(mask_outside_area,mask);
        Pt0(i) = Pt0(i) + sum(mpc.branch(mask,16))./mpc.baseMVA;

        
        %% Obtaining Ptie using equations (To use in nonlinear model)
        angle_bus = deg2rad(mpc.bus(network(i).to_bus(:,2),9));
        angle_nei = deg2rad(mpc.bus(network(i).to_bus(:,3),9));
        
        phi = atan2(network(i).to_bus(:,5),(network(i).to_bus(:,4)));

        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);
        V_bus = mpc.bus(network(i).to_bus(:,2),8);
        mask = network(i).to_bus(:,4) == 0;
        if any(mask)
            z_mod(mask) = z_mod(mask).*network(i).to_bus(mask,10);
            
        end
        
        V_nei = mpc.bus(network(i).to_bus(:,3),8);
        
       sum((V_bus.^2).*cos(phi)./z_mod + V_bus.*V_nei.*sin(angle_bus - angle_nei + phi - pi/2)./z_mod);

        %% Losses inside area

        mask_inside_area = bitand(ismember(mpc.branch(:,1),network(i).bus),ismember(mpc.branch(:,2),network(i).bus));
        Ploss(i) = sum(mpc.branch(mask_inside_area,14))./mpc.baseMVA + sum(mpc.branch(mask_inside_area,16))./mpc.baseMVA;
        PL0(i) = PL0(i) + Ploss(i);

    end

end

