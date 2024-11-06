function [A_global,B_global,C_global,D_global,W_global,C_mac,E,areas,network,bus_ss,ren_ss,k_ties] = get_global_ss(mpc,n_areas,flag_ren,network)
%get_global_ss   Given a MATPOWER case, this function handles grid partitioning and obtains a state-space
%   representation for the linearized dynamics as described in the thesis.
%
%   Inputs:                          
%       
%       mpc - The MATPOWER case struct after running OPF.
%       n_areas - Number of areas in the network.
%       flag_ren - A boolean flag that says if there are RESs in the system.
%       network - Vector of area objects with each area with buses assigned and the area parameters.
%
%
%   Outputs:
%
%       A_global - Continuous-time global transition matrix. 
%       B_global - Continuous-time global control matrix.
%       C_global - Continuous-time global output matrix.
%       D_global - Continuous-time global feedthrough matrix.
%       W_global - Continuous-time global disturbance matrix.
%       E - Sparsity matrix.
%       areas - Vector with the areas assigned to each bus.
%       network - Vector of area objects with each area with buses assigned and the area parameters.
%       bus_ss - Matrix that contains useful information about the state
%           space representation, such as the length of each area ss, number.
%       ren_ss - The index of RESs states.
%       k_ties - Tie-line coefficients.


    base_mva = 100;

    n_machine = sum(mpc.isgen);

    % Perform area partition unless it is passed as an argument
    if(nargin <= 3)
        % Area partition with only traditional generators
        if flag_ren
            areas = area_partitioning(mpc.branch,n_areas,mpc.gen(1:n_machine,1));
        else
            areas = area_partitioning(mpc.branch,n_areas,mpc.gen(:,1));
        end
        network = [];


        for i = 1:n_areas
            area_bus = find(areas==i);

            %Create area objects
            network = [network area(area_bus)];

            %Remove lines inside areas
            mask = bitand(ismember(mpc.branch(:,1),area_bus),ismember(mpc.branch(:,2),area_bus));
            mpc.branch(mask,:) = [];
            
        end
        
        % Convert machine parameters to the power base considered. For
        % reference see: Power System Control and Stability Paul M.Anderson
        % Chapter 2

        mpc.mac_con(:,16) = mpc.mac_con(:,16).*mpc.mac_con(:,3)/base_mva;
        
        % ADJUST 1/R TO STABILIZE THE SYSTEM
        mpc.tg_con(:,4) = mpc.tg_con(:,4)./5;


    
        %% Obtaining area inertia and damping
        for i=1:size(mpc.bus(:,1),1)
            
            mask = mpc.mac_con(:,2) == i;
            
            if any(mask)
                for j=1:n_areas
                    if ismember(i, network(j).bus)
                        network(j).inertia = network(j).inertia + mpc.mac_con(mask,16);
                        
                        network(j).machines = network(j).machines + 1;
                        network(j).tg_con = [network(j).tg_con ; mpc.tg_con(mpc.mac_con(mask,1),:)];
    
                        network(j).mac_base = [network(j).mac_base mpc.mac_con(mask,3)];

                        network(j).mac_bus = [network(j).mac_bus mpc.bus(i,1)];
                        network(j).mac_nr = [network(j).mac_nr ; find(mpc.mac_con(:,2) == i)];


                        network(j).freq_feedback = [network(j).freq_feedback; mpc.tg_con(mpc.mac_con(mask,1),4)];


                        break
                    end
                end
            end

        end


        for j=1:n_areas

            network(j).damping = sum(mpc.gen(network(j).mac_nr,2)./100); 
            network(j).inertia = network(j).inertia/network(j).machines;
            network(j).bias_factor =  network(j).damping + sum(network(j).tg_con(:,4));
        end

        %% Getting lines that connect to other areas
    
        lines = mpc.branch;
        
       
    
        for i = 1:size(network,2)
            network(i).to_bus = [];
            for j = 1:size(network,2)
                mask = bitand(ismember(lines(:,1),network(i).bus),ismember(lines(:,2),network(j).bus));
                if any(mask)
                    network(i).to_bus = [network(i).to_bus; repmat(j,sum(mask),1 ) lines(mask,:)];
                end
                
                mask = bitand(ismember(lines(:,1),network(j).bus),ismember(lines(:,2),network(i).bus));
                if any(mask)
                    line = lines(mask,:);
                    line(:,[1,2]) = line(:,[2,1]);
                    network(i).to_bus = [network(i).to_bus; repmat(j,sum(mask),1 ) line];
                end
            end
        
        end





    else
        % If the areas are assigned, just update the Damping coefficient
        areas = [];

        for j=1:n_areas
            network(j).damping = sum(mpc.gen(network(j).mac_nr,2)./100);           
        end
    end  



    
    
    %% Generate the ss for each area
    
    bus_ss = [];

    A_global = [];
    B_global = [];
    C_global = [];
    W_global = [];


    C_mac = [];

 
    
    n_areas = 1;
    if flag_ren
        ren_ss = zeros(1,size(mpc.gen(mpc.isolar_mask),1));
    else
        ren_ss = [];
    end



    %Each area is composed by the frequency state, the machine states, RESs
    %states and bus angle state

    for i=1:length(network)

        A_area = [];
        B_area = [];
        C_area = [];
        C_mac_area = [];
        
       

        n_ren = 0;

        % The state space representation differs according to the number of
        % machines

        for j = 1:size(network(i).tg_con,1)
            
            Ts = network(i).tg_con(j,6); Tc = network(i).tg_con(j,7); T3 = network(i).tg_con(j,8); T4 = network(i).tg_con(j,9); T5 = network(i).tg_con(j,10);

            a = Ts + Tc + T5;
            b = Tc*Ts+T5*(Tc+Ts);
            c = T5*Ts*Tc;
            d = T3 + T4;
            e = T3*T4;

            A_mech = [0 1 0;...
                 0 0 1;...
                 -1/c -a/c -b/c];
            B_mech = [0 0 1/c ]';

            C_mech = [1 d e];
            
            network(i).freq_feedback(j) = network(i).tg_con(j,4)/c;
           
            A_area = [A_area zeros(size(A_area,1),size(A_mech,2)) ; zeros(size(A_mech,1),size(A_area,2)) A_mech];
            B_area = [B_area zeros(size(B_area,1),size(B_mech,2)) ; zeros(size(B_mech,1),size(B_area,2)) B_mech];
            C_area = [C_area C_mech];
            C_mac_area = [C_mac_area zeros(size(C_mac_area,1),size(C_mech,2)) ; zeros(size(C_mech,1),size(C_mac_area,2)) C_mech];
    
        end
    

        % Machine Linear dynamics
        network(i).A = A_area;
        network(i).B = B_area;
        

        
        
        % Frequency state
        A_freq = -network(i).damping/network(i).inertia;
        B_freq = 1/network(i).inertia;
        
        % Governor frequncy feedback
        freq_feedback = zeros(size(B_area,1),1);
        freq_feedback(3:3:end) = -network(i).freq_feedback;
        
        A = [A_freq B_freq.*C_area ; freq_feedback A_area];

        B =  [ zeros(1,size(B_area,2)); B_area];
        
        
        
        network(i).mec_index = size(A_global,1)+2:size(A_global,1)+size(A,1);
        network(i).u_index = size(B_global,2)+1:size(B_global,2)+size(B_area,2);


        
        %Add renewables state(s)
        mask = ismember(network(i).bus,mpc.gen(mpc.isolar_mask));
        if any(mask) && flag_ren
            index = size(A_global,1) + size(A,1);
            
            generator_mask = zeros(size(mpc.gen,1),1);
            generator_mask(n_machine+1:end) = 1;

            network(i).res_bus = network(i).bus(mask);
            network(i).res_nr = find(bitand(ismember(mpc.gen(:,1),network(i).bus(mask)),generator_mask) == 1);
            

            n_ren = sum(mask);
            network(i).res = n_ren;

            network(i).res_index = index+1:1:index+n_ren;
            ren_ss(n_areas:n_areas+n_ren-1) = index+1:1:index+n_ren;
            
            for j = 1:n_ren
                A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
                A(1,end) = B_freq;
                Tpv = 1.8;
                Kpv = 1;
                A(end,end) = -1/Tpv;
    
                %u does not influence delta freq
                B =  [ B ; zeros(1,size(B,2)) ]; 

            end
        end

        %Add frequency integrator
        A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
        A(end,1) = 2*pi*60;
        B = [B; zeros(1,size(B,2))];

        
        W = zeros(size(B,1),n_ren+1);



        %Load disturbance
        W(1,1) = -B_freq;


        % Solar disturbance
        if n_ren ~= 0
            %location of ren states to area i
            index = ren_ss(n_areas):ren_ss(n_areas+n_ren-1);
            index  = index - size(A_global,1);
            W(index,2:end) = Kpv/Tpv*eye(n_ren);
            n_areas = n_areas + n_ren;
            network(i).C_res = zeros(1,size(A,1));
            network(i).C_res(end-n_ren:end-1) = 1;
            network(i).W_res = Kpv/Tpv*eye(n_ren);
        end
    

        
        % 4 outputs: Frequency, Sum of area mechanical power, Bus angle ,
        % tie-lines powerflow

        C = zeros(4,size(A,1));

        C(1,1) = 1;

        %Mechanical output
        C(2,:) = [0 C_area zeros(1,size(A,1)-size(C_area,2)-1)];

        %Angle
        C(3,end) = 1;

        
        network(i).C_mech = C(2,:);
        
        network(i).W = W;
        
        
        

        bus_ss = [bus_ss; i size(A,1) size(network(i).tg_con,1) n_ren];
        A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
        B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
        C_global = [C_global zeros(size(C_global,1), size(C,2)); zeros(size(C,1), size(C_global,2)) C];
        W_global = [W_global zeros(size(W_global,1), size(W,2)); zeros(size(W,1), size(W_global,2)) W];
        
        
        C_mac = [C_mac zeros(size(C_mac,1), size(C_mac_area,2)+2+n_ren); zeros(size(C_mac_area,1), size(C_mac,2)) zeros(size(C_mac_area,1),1) C_mac_area zeros(size(C_mac_area,1),1+n_ren)];
        
        
    end
    
    
    
    %% Add tie-lines to the state representation matrix A

    bus_sol = mpc.bus(:,9);
    bus_sol = deg2rad(bus_sol);



    %Sparsity matrix
    E = zeros(size(B_global,2),size(A_global,1));

    
    index_ss = cumsum(bus_ss(:,2))+1;
    index_ss = [1 ; index_ss];
    C_index = 4:4:size(C_global,1);
    tg_size = cumsum([ 1 ; bus_ss(:,3)]);
    for i = 1:size(network,2)


        neighbours = network(i).to_bus;

        E(tg_size(i):tg_size(i+1)-1,index_ss(i):index_ss(i+1)-1) = ones(bus_ss(i,3),bus_ss(i,2));
        T_i = 0;


        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);
        mask = bitand(network(i).to_bus(:,4) == 0, network(i).to_bus(:,10) ~= 0);
        if any(mask)
            z_mod(mask) = z_mod(mask).*network(i).to_bus(mask,10);
            
        end

        V_area = mpc.bus(network(i).to_bus(:,2),8);
        V_neighbour = mpc.bus(network(i).to_bus(:,3),8);

        Tij = V_area.*V_neighbour.*cos(bus_sol(network(i).to_bus(:,2)) - bus_sol(network(i).to_bus(:,3)))./z_mod;
        network(i).to_bus = [network(i).to_bus Tij];


        for j = 1:size(neighbours,1)

            T_ji =  Tij(j);
       
            
            A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) =  A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) +T_ji/network(i).inertia;
            

            E(tg_size(i):tg_size(i+1)-1,index_ss(neighbours(j,1)):index_ss(neighbours(j,1)+1)-1) = ones(bus_ss(i,3),bus_ss(neighbours(j,1),2));


           
            T_i =  T_i + T_ji;

            C_global(C_index(i),index_ss(neighbours(j,1)+1)-1) = C_global(C_index(i),index_ss(neighbours(j,1)+1)-1) -T_ji;
        end
       
        
        A_global(index_ss(i), index_ss(i+1)-1) = -T_i/network(i).inertia;
        k_ties(i) = T_i;
        C_global(C_index(i),index_ss(i+1)-1) = T_i;

    end

   

    D_global = zeros(size(C_global,1),size(B_global,2));

   
    


end

