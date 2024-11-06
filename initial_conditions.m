function [x0,u0,Pt0,PL0,Ploss] = initial_conditions(ss_dim,n_machines,bus_ss,network,mpc)
%initial_conditions   Given a MATPOWER case, this function obtains the initial conditions for state-space
%   representation .
%
%   Inputs:                          
%       
%       ss_dim - The dimension of the state-space.
%       n_machines - Number of machines in the network.
%       network - Vector of area objects with each area with buses assigned and the area parameters.
%       mpc - The MATPOWER case struct after running OPF.
%
%
%   Outputs:
%
%       x0 - The MATPOWER case struct after running OPF.
%       u0 - Number of areas in the network (deprecated).
%       Pt0 - A boolean flag that says if there are RESs in the system.
%       PL0 - Vector of area objects with each area with buses assigned and the area parameters (deprecated).
%       Ploss- Vector of area objects with each area with buses assigned and the area parameters (deprecated). 


    x0 = zeros(ss_dim,1);
    u0 = zeros(n_machines,1);
    PL0 = zeros(length(network),1);
    Pt0 = zeros(length(network),1);
    Ploss = zeros(length(network),1);

    angle_index = cumsum(bus_ss);
    freq_index = [1 ; angle_index+1];
    
    x0(freq_index(1:end-1),1) = 1;
   

    %x0(angle_index,1) = pi/30;
    %x0(freq_index(1:end-1),1) = 0;
    % x0(angle_index,1) = 0;
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


        %% Testing losses inside area

        mask_inside_area = bitand(ismember(mpc.branch(:,1),network(i).bus),ismember(mpc.branch(:,2),network(i).bus));
        Ploss(i) = sum(mpc.branch(mask_inside_area,14))./mpc.baseMVA + sum(mpc.branch(mask_inside_area,16))./mpc.baseMVA;
        PL0(i) = PL0(i) + Ploss(i);

    end

end

