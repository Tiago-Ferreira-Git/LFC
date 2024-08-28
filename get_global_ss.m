function [A_global,B_global,C_global,D_global,W_global,machine_ss,C_mac,u,E,areas,network,bus_ss,ren_ss,E_fs] = get_global_ss(mpc,n_areas,flag_ren,debug,network)
    
    base_mva = 100;
    ren_data = load('data\solar.mat');
    ren_data = ren_data.data;
    n_machine = size(mpc.gen,1)-length(ren_data.bus);
    if(nargin <= 4)
        if flag_ren
            areas = area_partitioning(mpc.branch,n_areas,mpc.gen(1:n_machine,1));
        else
            areas = area_partitioning(mpc.branch,n_areas,mpc.gen(:,1));
        end
        network = [];
        for i = 1:n_areas
            area_bus = find(areas==i);
            network = [network area(area_bus)];
            %Remove lines inside areas
            mask = bitand(ismember(mpc.branch(:,1),area_bus),ismember(mpc.branch(:,2),area_bus));
            mpc.branch(mask,:) = [];
            
        end
        
        
        mpc.mac_con(:,16) = mpc.mac_con(:,16).*mpc.mac_con(:,3)/base_mva;
        mpc.mac_con(:,17) = mpc.mac_con(:,17).*mpc.mac_con(:,3)/base_mva;
    



    
        %% obtaining area inertia and damping
        for i=1:size(mpc.bus(:,1),1)
            
            mask = mpc.mac_con(:,2) == i;
            
            if any(mask)
                for j=1:n_areas
                    if ismember(i, network(j).bus)
                        network(j).inertia = network(j).inertia + mpc.mac_con(mask,16);
                        
                        %network(j).damping = network(j).damping + mpc.mac_con(mask,17);
                        network(j).machines = network(j).machines + 1;
                        network(j).tg_con = [network(j).tg_con ; mpc.tg_con(mpc.mac_con(mask,1),:)];
    
                        network(j).mac_base = [network(j).mac_base mpc.mac_con(mask,3)];

                        network(j).mac_bus = [network(j).mac_bus mpc.bus(i,1)];
                        network(j).mac_nr = [network(j).mac_nr ; find(mpc.mac_con(:,2) == i)];


                        break
                    end
                end
            end

        end
        
        
        for j=1:n_areas
            %network(j).damping = network(j).damping/network(j).machines;
            network(j).damping = sum(mpc.gen(network(j).mac_nr,2)./100);
            network(j).inertia = network(j).inertia/network(j).machines;
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
        areas = [];
    end  



    
    
    %% Generate the ss for each area
    
    bus_ss = [];
    A_global = [];
    B_global = [];
    C_global = [];
    C_ren = [];
    W_global = [];
    %W_global_mech = [];

    machine_ss = [];

    C_mac = [];
    
    if flag_ren 
        n_res = size(ren_data.bus,1);
    else
        n_res = 0;
    end

   

    u = [];
    s = tf('s');
    n_areas = 1;
    if flag_ren
        ren_ss = zeros(1,size(ren_data.bus,1));
    else
        ren_ss = [];
    end


    n_machines = 0;
    for i=1:length(network)
        u_area = [];


        A_area = [];
        B_area = [];
        C_area = [];
        C_mac_area = [];
        

        n_ren = 0;
        for j = 1:size(network(i).tg_con,1)
            n_machines = n_machines + 1;
            %p_mech ss
            
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

            %%Applying machine loss faults
            if( Ts == 0 && Tc == 0 && T3 == 0 && T4 == 0 && T5 == 0 )
                A_mech = [0 0 0;...
                     0 0 0;...
                     0 0 0];
                B_mech = [0 0 0 ]';
    
                C_mech = [0 0 0];
                n = size(A_global,1)+1+(j-1)*3;
                machine_ss = [machine_ss ; (n:n+3-1)' ones(3,1).*n_machines];
            end
            
           
            A_area = [A_area zeros(size(A_area,1),size(A_mech,2)) ; zeros(size(A_mech,1),size(A_area,2)) A_mech];
            B_area = [B_area zeros(size(B_area,1),size(B_mech,2)) ; zeros(size(B_mech,1),size(B_area,2)) B_mech];
            C_area = [C_area C_mech];
            C_mac_area = [C_mac_area zeros(size(C_mac_area,1),size(C_mech,2)) ; zeros(size(C_mech,1),size(C_mac_area,2)) C_mech];
    
        end
    

        % Non-linear Simulation
        network(i).A_mech = A_area;
        network(i).B_mech = B_area;
        

        
        
    
        A_freq = -network(i).damping/network(i).inertia;
        B_freq = 1/network(i).inertia;
        
        freq_feedback = zeros(size(B_area,1),1);
        freq_feedback(3:3:end) = -network(i).tg_con(:,4);
        
        A = [A_freq B_freq.*C_area ; freq_feedback A_area];

        B =  [ zeros(1,size(B_area,2)); B_area];
        
        
        
        %Add renewables state
        mask = ismember(network(i).bus,ren_data.bus);
        if any(mask) && flag_ren
            index = size(A_global,1) + size(A,1);
            
            generator_mask = zeros(size(mpc.gen,1),1);
            generator_mask(n_machine+1:end) = 1;

            network(i).res_bus = network(i).bus(mask);
            network(i).res_nr = find(bitand(ismember(mpc.gen(:,1),network(i).bus(mask)),generator_mask) == 1);
            

            n_ren = sum(mask);
            network(i).res = n_ren;
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

        %Add error integrator 

        %unique_lines = unique(network(i).to_bus,'rows');
        %mask = bitor(ismember(network(i).bus,unique_lines(:,2),'rows'),ismember(network(i).bus,unique_lines(:,3),'rows'));
        


        A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
        if debug == 1
            A(end,1) = -1; 
        else
            A(end,1) = -2*pi*60;
        end
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
    

        


        C = zeros(3,size(A,1));

        C(1,1) = 1;
        C(2,:) = [0 C_area zeros(1,size(A,1)-size(C_area,2)-1)];
        %C(3,end) = 1;
        C(3,end) = -1;

        network(i).C_mech = C(2,:);

        network(i).W = W;
        
        network(i).C = zeros(2+network(i).machines+n_ren,size(A,1));

        network(i).A = A;
        network(i).B = B;
        network(i).C(1,1) = 1;
        network(i).C(2:end-1-n_ren,2:end-1-n_ren) = C_mac_area;
        network(i).C(end,end) = -1;

        

        bus_ss = [bus_ss; i size(A,1) size(network(i).tg_con,1) n_ren];
        A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
        B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
        C_global = [C_global zeros(size(C_global,1), size(C,2)); zeros(size(C,1), size(C_global,2)) C];
        W_global = [W_global zeros(size(W_global,1), size(W,2)); zeros(size(W,1), size(W_global,2)) W];
        
        %W_global_mech = [W_global_mech zeros(size(W_global_mech,1), size(W_mech,2)); zeros(size(W_mech,1), size(W_global_mech,2)) W_mech];

        C_mac = [C_mac zeros(size(C_mac,1), size(C_mac_area,2)+2+n_ren); zeros(size(C_mac_area,1), size(C_mac,2)) zeros(size(C_mac_area,1),1) C_mac_area zeros(size(C_mac_area,1),1+n_ren)];
        
        
    end
    
    
    
    %%

    bus_sol = mpc.bus(:,9);
    bus_sol = deg2rad(bus_sol);


    %Add tie-lines
    Adjacency_matrix = zeros(n_areas); 
    E = zeros(size(B_global,2),size(A_global,1));
    E_fs = zeros(size(B_global,2),size(A_global,1));

    %L is the E matrix for the integrators
    L = zeros(size(B_global,2),n_areas);
    index_ss = cumsum(bus_ss(:,2))+1;
    index_ss = [1 ; index_ss];
    C_index = 3:4:size(C_global,1);
    tg_size = cumsum([ 1 ; bus_ss(:,3)]);
    


    for i = 1:size(network,2)


        neighbours = network(i).to_bus;

        E(tg_size(i):tg_size(i+1)-1,index_ss(i):index_ss(i+1)-1) = ones(bus_ss(i,3),bus_ss(i,2));
        E_fs(tg_size(i):tg_size(i+1)-1,index_ss(i):index_ss(i+1)-1) = ones(bus_ss(i,3),bus_ss(i,2));
        T_i = 0;
        z_mod = vecnorm(network(i).to_bus(:,4:5),2,2);
        for j = 1:size(neighbours,1)
            
            % Global Matrix A 

            if debug == 1
                T_ji = 2*pi*60*cos(bus_sol(neighbours(j,3)) - bus_sol(neighbours(j,2)))/(z_mod(j)); 
            else
                T_ji = cos(bus_sol(neighbours(j,3)) - bus_sol(neighbours(j,2)))/(z_mod(j));
            end

            % T_ji =  T_ji/10000;

            %Default
            %A_global(index_ss(i+1)-1, index_ss(neighbours(j,1))) =  A_global(index_ss(i+1)-1, index_ss(neighbours(j,1)) ) -T_ji;
            
            %Changin ptie for error integral
            A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) =  A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) -T_ji/network(i).inertia;
            

            E(tg_size(i):tg_size(i+1)-1,index_ss(neighbours(j,1)):index_ss(neighbours(j,1)+1)-1) = ones(bus_ss(i,3),bus_ss(neighbours(j,1),2));


            Adjacency_matrix(i,neighbours(j,1)) = 1;
            L(1+tg_size(i):tg_size(i)+network(i).machines,neighbours(j,1)) = 1;

            T_i =  T_i + T_ji;
        end
        A_global(index_ss(i), index_ss(i+1)-1) = T_i/network(i).inertia;


        if size(Adjacency_matrix,1) ~= 1
            Adjacency_matrix(i,i) = size(unique(neighbours(:,1)),1);
        end
        L(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+network(i).machines,i) = 1;
        
        



        areas = unique(neighbours(:,1));
        n = size(network(i).A,1);
        network(i).freq_index = [1 1+ n];
        for j = 1:size(areas,1)
            mask = neighbours(:,1) == areas(j);
            % Local Matrix A 

            if debug == 1
                T_ji = sum(2*pi*60*cos(bus_sol(neighbours(mask,3)) - bus_sol(neighbours(mask,2)))./(z_mod(mask))); 
            else
                T_ji = sum(cos(bus_sol(neighbours(mask,3)) - bus_sol(neighbours(mask,2)))./(z_mod(mask)));
            end
            
            network(i).C = [network(i).C zeros(size(network(i).C,1),size(network(areas(j)).A,1)); zeros(size(network(areas(j)).C,1),size(network(i).A,2)) network(areas(j)).C];
            network(i).C_mech = [network(i).C_mech zeros(size(network(i).C_mech,1),size(network(areas(j)).A,1)); zeros(size(network(areas(j)).C_mech,1),size(network(i).A,2)) network(areas(j)).C_mech];
            

            network(i).A = [network(i).A zeros(size(network(i).A,1),size(network(areas(j)).A,1)) ; zeros(size(network(areas(j)).A,1),size(network(i).A,1)) network(areas(j)).A];
            network(i).A(1,end) = -T_ji/network(i).inertia;


            %Neighbour
            network(i).A(n+1,1) = -T_ji/network(areas(j)).inertia;
            network(i).A(n+1,end) = T_ji/network(areas(j)).inertia;
            T_i =  T_i + T_ji;


            network(i).B = [network(i).B zeros(size(network(i).B,1),size(network(areas(j)).B,2)) ; zeros(size(network(areas(j)).B,1),size(network(i).B,2)) network(areas(j)).B];
            network(i).freq_index = [network(i).freq_index network(i).freq_index(end)+ size(network(areas(j)).A,1)];
        end

        network(i).A(1,n) = T_i/network(i).inertia;
        [network(i).A,network(i).B] = discrete_dynamics(network(i).A,network(i).B,zeros(size(network(i).A,1),1),0.1);

        network(i).C = eye(size(network(i).A));
        
        G = eye(size(network(i).A));
    
    
        Q = eye(size(network(i).A));
        R = 10*eye(size(network(i).C,1));
    
        network(i).L = dlqe(network(i).A,G,network(i).C,Q,R);



    end
    
    


   

    D_global = zeros(size(C_global,1),size(B_global,2));

    %L = Adjacency_matrix;

    % g = graph(Adjacency_matrix);
    % figure
    % plot(g)
    % set(gcf,'renderer','Painters');
    % title='./fig/graph.png';
    % saveas(gca,title,'png');


   


    % s = ones(1,size(A_global,1));
    % s(1,cumsum(bus_ss(1:end-1,2))) = 1000;
    % S = diag(s);
    % A_global = S\A_global*S;
    % B_global = S\B_global; 
    % C_global = C_global*S; 


    E = logical(E);
    E_fs = logical(E);


end

