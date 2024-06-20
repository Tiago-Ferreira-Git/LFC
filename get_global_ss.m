function [A_global,B_global,C_global,D_global,W_global,u,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,n_areas,flag_ren,flag_integrator,network)
    
    base_mva = 100;
    ren_data = load('data\solar.mat');
    ren_data = ren_data.data;
    if(nargin <= 5)
        
        areas = area_partitioning(g.line,n_areas,g.mac.mac_con(:,2));
        
        network = [];
        for i = 1:n_areas
            area_bus = find(areas==i);
            network = [network area(area_bus)];
            %Remove lines inside areas
            mask = bitand(ismember(g.line(:,1),area_bus),ismember(g.line(:,2),area_bus));
            g.line(mask,:) = [];
            
        end
        
        
        g.mac.mac_con(:,16) = g.mac.mac_con(:,16).*g.mac.mac_con(:,3)/base_mva;
        g.mac.mac_con(:,17) = g.mac.mac_con(:,17).*g.mac.mac_con(:,3)/base_mva;
    
    
        %% obtaining area inertia and damping
        for i=1:size(g.bus.bus_int,1)
            
            mask = g.mac.mac_con(:,2) == i;
            if any(mask)
                for j=1:n_areas
                    if ismember(i, network(j).bus)
                        network(j).inertia = network(j).inertia + g.mac.mac_con(mask,16);
                        network(j).damping = network(j).damping + g.mac.mac_con(mask,17);
                        network(j).machines = network(j).machines + 1;
                        network(j).tg_con = [network(j).tg_con ; g.tg.tg_con(g.mac.mac_con(mask,1),:)];
    
                        network(j).mac_base = [network(j).mac_base g.mac.mac_con(mask,3)];

                        network(j).mac_bus = [network(j).mac_bus g.bus.bus_int(i)];
                        break
                    end
                end
            end

        end
        for j=1:n_areas
            network(j).damping = network(j).damping/network(j).machines;
            %network(j).tg_con = sum(network(j).tg_con,1)./network(j).machines;
            network(j).inertia = network(j).inertia/network(j).machines;
            % network(j).inertia  = network(j).inertia * ( 0.7 + rand*(1.5-0.7));
    
            % Retrieving lmod index
            [~,~,index_lmod] = intersect(network(j).bus,g.lmod.lmod_con(:,2),'stable');
            network(j).lmod = index_lmod;
        end
        %% Getting lines that connect to other areas
    
        lines = g.line;
        
       
    
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
    W_global = [];
    

    u = [];
    s = tf('s');
    n_areas = 1;
    if flag_ren
        ren_ss = zeros(1,size(ren_data.bus,1));
    else
        ren_ss = [];
    end
    for i=1:length(network)
        u_area = [];


        A_area = [];
        B_area = [];
        C_area = [];
        n_ren = 0;
        for j = 1:size(network(i).tg_con,1)
    
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
            end
            
           
            A_area = [A_area zeros(size(A_area,1),size(A_mech,2)) ; zeros(size(A_mech,1),size(A_area,2)) A_mech];
            B_area = [B_area zeros(size(B_area,1),size(B_mech,2)) ; zeros(size(B_mech,1),size(B_area,2)) B_mech];
            C_area = [C_area C_mech];
    
        end
    

        
    
        A_freq = -network(i).damping/network(i).inertia;
        B_freq = 1/network(i).inertia;
        
        freq_feedback = zeros(size(B_area,1),1);
        freq_feedback(3:3:end) = -network(i).tg_con(:,4);
        
        A = [A_freq B_freq.*C_area ; freq_feedback A_area];

        B =  [ zeros(1,size(B_area,2)); B_area];
        
        
        
        %Add renewables state
        if any(ismember(network(i).bus,ren_data.bus)) && flag_ren
            index = size(A_global,1) + size(A,1);
            
            n_ren = sum(ismember(network(i).bus,ren_data.bus)) ;
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
        A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
        A(end,1) = -1;
        B = [B; zeros(1,size(B_area,2))];

        
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
        end


        C = zeros(4,size(A,1));

        C(1,1) = 1;
        C(2,:) = [0 C_area zeros(1,size(A,1)-size(C_area,2)-1)];
        %C(3,end) = 1;
        C(4,end) = -1;

        

        bus_ss = [bus_ss; i size(A,1) size(network(i).tg_con,1) n_ren B_freq];
        A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
        B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
        C_global = [C_global zeros(size(C_global,1), size(C,2)); zeros(size(C,1), size(C_global,2)) C];
        W_global = [W_global zeros(size(W_global,1), size(W,2)); zeros(size(W,1), size(W_global,2)) W];
    end
    
    
    %%

    bus_sol = bus;
    bus_sol(:,3) = deg2rad(bus(:,3));


    %Add tie-lines
    Adjacency_matrix = zeros(n_areas); 
    E = zeros(size(B_global,2),size(A_global,1));

    %L is the E matrix for the integrators
    L = zeros(size(B_global,2),n_areas);
    index_ss = cumsum(bus_ss(:,2))+1;
    index_ss = [1 ; index_ss];
    C_index = 3:4:size(C_global,1);
    tg_size = cumsum([ 1 ; bus_ss(:,3)]);
    %tg_size = [1 ; tg_size];
    for i = 1:size(network,2)


        neighbours = network(i).to_bus;

        E(tg_size(i):tg_size(i+1)-1,index_ss(i):index_ss(i+1)-1) = ones(bus_ss(i,3),bus_ss(i,2));
        T_i = 0;
        for j = 1:size(neighbours,1)
            T_ji = 377*cos(bus_sol(g.bus.bus_int(neighbours(j,3)),3) - bus_sol(g.bus.bus_int(neighbours(j,2)),3))/(neighbours(j,5));

            % T_ji =  T_ji/10000;

            %Default
            %A_global(index_ss(i+1)-1, index_ss(neighbours(j,1))) =  A_global(index_ss(i+1)-1, index_ss(neighbours(j,1)) ) -T_ji;
            
            %Changin ptie for error integral
            A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) =  A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) -T_ji*bus_ss(i,5);
            

            E(tg_size(i):tg_size(i+1)-1,index_ss(neighbours(j,1)):index_ss(neighbours(j,1)+1)-1) = ones(bus_ss(i,3),bus_ss(neighbours(j,1),2));


            Adjacency_matrix(i,neighbours(j,1)) = 1;
            L(1+tg_size(i):tg_size(i)+network(i).machines,neighbours(j,1)) = 1;

            T_i =  T_i + T_ji;

            C_global(C_index(i),index_ss(neighbours(j,1)+1)-1) = C_global(C_index(i),index_ss(neighbours(j,1)+1)-1) -T_ji;
        end
        if size(Adjacency_matrix,1) ~= 1
            Adjacency_matrix(i,i) = size(unique(neighbours(:,1)),1);
        end
        L(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+network(i).machines,i) = 1;
        
        %Default
        %A_global(index_ss(i+1)-1, index_ss(i) ) = T_i;
        
        A_global(index_ss(i), index_ss(i+1)-1) = T_i*bus_ss(i,5);
        C_global(C_index(i),index_ss(i+1)-1) = T_i;

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


    


end

