function [A_global,B_global,C_global,D_global,W_global,u,E,L,areas,network,bus_ss,rows_NC] = get_global_ss(g,bus,nr_areas,flag_u,flag_ren)

    ren_data = load('data\solar.mat');
    ren_data = ren_data.data;

    if ~flag_u
        u = nan;
    end
    base_mva = 100;
    k = nr_areas; %three areas
    
    areas = area_partitioning(g.line,k,g.mac.mac_con(:,2));
    
    network = [];
    for i = 1:k
        area_bus = find(areas==i);
        network = [network area(area_bus)];
        
    end
    
    
    g.mac.mac_con(:,16) = g.mac.mac_con(:,16).*g.mac.mac_con(:,3)/base_mva;
    g.mac.mac_con(:,17) = g.mac.mac_con(:,17).*g.mac.mac_con(:,3)/base_mva;
    
    
    %% obtaining area inertia and damping
    for i=1:size(g.bus.bus_int,1)
        
        mask = g.mac.mac_con(:,2) == i;
        if any(mask)
            for j=1:k
                if ismember(i, network(j).bus)
                    network(j).inertia = network(j).inertia + g.mac.mac_con(mask,16);
                    network(j).damping = network(j).damping + g.mac.mac_con(mask,17);
                    network(j).machines = network(j).machines + 1;
                    network(j).tg_con = [network(j).tg_con ; g.tg.tg_con(g.mac.mac_con(mask,1),:)];
                    if flag_u   
                        network(j).tg_sig = [network(j).tg_sig  ;g.tg.tg_sig(g.mac.mac_con(mask,1),:)];
                    end
                    network(j).mac_base = [network(j).mac_base g.mac.mac_con(mask,3)];
                    break
                end
            end
        end
    
    end
    
    for j=1:k
        
        network(j).damping = network(j).damping/network(j).machines;
        %network(j).tg_con = sum(network(j).tg_con,1)./network(j).machines;
        network(j).inertia = network(j).inertia/network(j).machines;
        % Retrieving lmod index
        [~,~,index_lmod] = intersect(network(j).bus,g.lmod.lmod_con(:,2),'stable');
        network(j).lmod = index_lmod;
    end
    
    %% Getting lines that connect to other areas
    
    lines = g.line;
    
   

    %neighbour = 1:size(network,2);
    %parse lines such that it contain only the lines that connect each area
    for j = 1:k
        for i = 1:size(network(j).bus,1)
    
            mask = lines(:,1) == network(j).bus(i);
            if any(mask)
    
                mask_neighbour = ismember(lines(:,2), network(j).bus);
                lines(bitand(mask,mask_neighbour),:) = [];
    
            end
    
    
            mask = lines(:,2) == network(j).bus(i);
            if any(mask)
    
                mask_neighbour = ismember(lines(:,1), network(j).bus);
                lines(bitand(mask,mask_neighbour),:) = [];
    
            end
        end
    end
    
    
    %%
    for i = 1:size(network,2)
        to_areas = 1:k;
        mask = to_areas == i;
        to_areas = to_areas(~mask);
        % each line to an area
        for j = 1:size(lines,1)
            mask = network(i).bus == lines(j,1);
    
            if any(mask)
                for n=1:k-1
                    mask_neighbour =  network(to_areas(n)).bus == lines(j,2);
                    
                    if any(mask_neighbour)
                        network(i).to_bus = [network(i).to_bus; to_areas(n) lines(j,:)];
                    end
                        
            
                end
    
            end
    
            mask = network(i).bus == lines(j,2);
    
            if any(mask)
                for n=1:k-1
                    mask_neighbour =  network(to_areas(n)).bus == lines(j,1);
                    
                    if any(mask_neighbour)
                        to_line = lines(j,:);
                        to_line(:,[2,1]) = to_line(:,[1,2]);
                        network(i).to_bus = [network(i).to_bus; to_areas(n) to_line];
                        %network(i).to_bus(:,[2,3]) = network(i).to_bus(:,[3,2]);
                    end
                        
            
                end
    
            end
        end
    
    end
    
    
    %% Generate the ss for each area
    bus_ss = [];
    A_global = [];
    B_global = [];
    C_global = [];
    W_global = [];

    ren_ss = zeros(1,size(ren_data.bus,1));
    

    u = [];
    s = tf('s');
    k = 1;
    for i=1:length(network)
        u_area = [];
        % if flag_u
        %     if( 0 < size(network(i).lmod,1))
        %         u_area = [u_area ; -sum(g.lmod.lmod_sig(network(i).lmod,:),1)];
        %     else
        %         u_area = [u_area ; zeros(1,size(g.lmod.lmod_sig(1,:),2))];
        %     end
        % end

        A_area = [];
        B_area = [];
        C_area = [];

        for j = 1:size(network(i).tg_con,1)
    
            %p_mech ss
            
            Ts = network(i).tg_con(j,6); Tc = network(i).tg_con(j,7); T3 = network(i).tg_con(j,8); T4 = network(i).tg_con(j,9); T5 = network(i).tg_con(j,10);
            % T3 = 0;
            % T4 = 0;
            % T5=0;

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

            %sys = (1/(1 + s*Ts)) * ((1+s*T3)/(1+s*Tc)) * ((1+s*T4)/(1+s*T5));

            %[A_mech,B_mech,C_mech,~] = tf2ss(sys.Numerator{1},sys.Denominator{1});
            
            %Control signal input
            if flag_u
                u_area = [u_area ;  network(i).tg_sig(j,:)];
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
        
    

        
        
        
        %Add renewables state
        if any(ismember(network(i).bus,ren_data.bus)) && flag_ren
            A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
            A(1,end) = B_freq;
            Tpv = 1.8;
            Kpv = 1;
            A(end,end) = -1/Tpv;

            %u does not influence delta freq, p_tie
            B =  [ zeros(1,size(B_area,2)); B_area ; zeros(1,size(B_area,2)) ; zeros(1,size(B_area,2)) ];

            %Load disturbance and solar
            W = zeros(size(B,1),2);
            W(1) = B_freq;
            W(end-1) = Kpv/Tpv;
            ren_ss(k) = size(A_global,1)+size(A,1);
            k = k+1;
        else 
            %u does not influence delta freq, p_tie
            B =  [ zeros(1,size(B_area,2)); B_area ; zeros(1,size(B_area,2)) ];

            %Load disturbance and solar
            W = zeros(size(B,1),2);
            W(1) = B_freq;
        end

        
        %Add tie-line state
        A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
        A(1,end) = -B_freq;

        

        
        
    



        C = zeros(3,size(A,1));
        C(1,1) = 1;
        C(2,:) = [0 C_area zeros(1,size(A,1)-size(C_area,2)-1)];
        C(3,end) = 1;
        
        if flag_u
            u = [u;u_area];
        end
        

        bus_ss = [bus_ss; i size(A,1) size(network(i).tg_con,1)];
        A_global = [A_global zeros(size(A_global,1), size(A,2)) ; zeros(size(A,1), size(A_global,2)) A ];
        B_global = [B_global zeros(size(B_global,1), size(B,2)); zeros(size(B,1), size(B_global,2)) B];
        C_global = [C_global zeros(size(C_global,1), size(C,2)); zeros(size(C,1), size(C_global,2)) C];
        W_global = [W_global zeros(size(W_global,1), size(W,2)); zeros(size(W,1), size(W_global,2)) W];
    end
    
    
    %%

    bus_sol = bus;
    bus_sol(:,3) = deg2rad(bus(:,3));
    
    
    %Add tie-lines
    Adjacency_matrix = zeros(nr_areas); 
    E = zeros(size(B_global,2),size(A_global,1));

    %L is the E matrix for the integrators
    L = zeros(size(B_global,2),nr_areas);
    for i = 1:size(network,2)
        
    
        neighbours = network(i).to_bus;

        
        E(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i,3)),1+sum(bus_ss(1:i-1,2)):sum(bus_ss(1:i,2))) = ones(size(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i,3)),2),size(1+sum(bus_ss(1:i-1,2)):sum(bus_ss(1:i,2)),2));
        T_i = 0;
        for j = 1:size(neighbours,1)
            T_ji = 377*cos(bus_sol(g.bus.bus_int(neighbours(j,3)),3) - bus_sol(g.bus.bus_int(neighbours(j,2)),3))/(neighbours(j,5));
    
            
            
            A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:neighbours(j,1)-1,2))+1 ) =  A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:neighbours(j,1)-1,2))+1 ) -T_ji;
            %A_global(sum(bus_ss(1:i,2))+1 ,sum(bus_ss(1:i,2)) ) = T_ji;
    
            E(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i,3)),1+sum(bus_ss(1:neighbours(j,1)-1,2)):sum(bus_ss(1:neighbours(j,1),2))) = ones(size(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i,3)),2),size(1+sum(bus_ss(1:neighbours(j,1)-1,2)):sum(bus_ss(1:neighbours(j,1),2)),2));
            
            Adjacency_matrix(i,neighbours(j,1)) = 1;

            L(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+network(i).machines,neighbours(j,1)) = 1;
            
            T_i =  T_i + T_ji;
        end
        Adjacency_matrix(i,i) = size(unique(neighbours(:,1)),1);
        L(1+sum(bus_ss(1:i-1,3)):sum(bus_ss(1:i-1,3))+network(i).machines,i) = 1;
        A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:i-1,2))+1 ) = T_i;
    
    end
   
    D_global = zeros(size(C_global,1),size(B_global,2));

    %L = Adjacency_matrix;

    g = graph(Adjacency_matrix);
    plot(g)
    set(gcf,'renderer','Painters');
    title='./fig/graph.png';
    saveas(gca,title,'png');

    rows_NC = ren_ss;
end

