function [A_global,B_global,C_global,D_global,W_global,u,E,areas,network,bus_ss,ren_ss] = get_global_ss(g,bus,nr_areas,flag_ren,flag_integrator)

    ren_data = load('data\solar.mat');
    ren_data = ren_data.data;

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
        % network(j).inertia  = network(j).inertia * ( 0.7 + rand*(1.5-0.7));

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
    

    u = [];
    s = tf('s');
    k = 1;
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

            %sys = (1/(1 + s*Ts)) * ((1+s*T3)/(1+s*Tc)) * ((1+s*T4)/(1+s*T5));

            %[A_mech,B_mech,C_mech,~] = tf2ss(sys.Numerator{1},sys.Denominator{1});
            
            %Control signal input

            
           
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
            ren_ss(k:k+n_ren-1) = index+1:1:index+n_ren;
            
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

        if flag_integrator
            %Add error integrator 
            A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
            A(end,1) = -1;
            B = [B; zeros(1,size(B_area,2))];
        end

        %Add tie-line state
        % A = [A zeros(size(A,1),1); zeros(1,size(A,1)+1)];
        % A(1,end) = -B_freq;
        % B = [B; zeros(1,size(B_area,2))];


        
        W = zeros(size(B,1),n_ren+1);

        %Load disturbance
        W(1,1) = -B_freq;


        % Solar disturbance
        if n_ren ~= 0
            %location of ren states to area i
            index = ren_ss(k):ren_ss(k+n_ren-1);
            index  = index - size(A_global,1);
            W(index,2:end) = Kpv/Tpv*eye(n_ren);
            k = k + n_ren;
        end


        if flag_integrator
            C = zeros(4,size(A,1));

            C(1,1) = 1;
            C(2,:) = [0 C_area zeros(1,size(A,1)-size(C_area,2)-1)];
            %C(3,end) = 1;
            C(4,end) = -1;
        else
            C = zeros(3,size(A,1));
            C(1,1) = 1;
            C(2,:) = [0 C_area zeros(1,size(A,1)-size(C_area,2)-1)];
            C(3,end) = 1;
        end

        

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
    Adjacency_matrix = zeros(nr_areas); 
    E = zeros(size(B_global,2),size(A_global,1));

    %L is the E matrix for the integrators
    L = zeros(size(B_global,2),nr_areas);
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

    %Last p_tie is equal to p_tie = -sum(ptie,i) . Since the contribution to the
    %frequency is -p_tie*damping/inertia ->  sum(p_tie,i)*damping/inertia 
    % A_global(index_ss(end-1),cumsum(bus_ss(1:end-1,2))) = -A_global(index_ss(end-1),size(A_global,1));
    % A_global(:,end) = [];
    % A_global(end,:) = [];
    % bus_ss(end,2) = bus_ss(end,2)-1;
    % B_global(end,:) = [];
    % C_global(:,end) = [];
    % C_global(end,:) = zeros(1,size(A_global,2));
    % C_global(end,cumsum(bus_ss(1:end-1,2))) = -1;
    % W_global(end,:) = [];
    % E(:,end) = [];

    D_global = zeros(size(C_global,1),size(B_global,2));

    %L = Adjacency_matrix;

    g = graph(Adjacency_matrix);
    plot(g)
    set(gcf,'renderer','Painters');
    title='./fig/graph.png';
    saveas(gca,title,'png');


    % [~,S,~] = svd(A_global);
    
    %S = diag(max(abs(A_global),[],2));
    % A_global = S\A_global*S;
    % B_global = S\B_global; 
    % C_global = C_global*S; 




    % s = ones(1,size(A_global,1));
    % s(1,cumsum(bus_ss(1:end-1,2))) = 1000;
    % S = diag(s);
    % A_global = S\A_global*S;
    % B_global = S\B_global; 
    % C_global = C_global*S; 


    
    % cond(A_global)
    % 
    % max(svd(A_global)) 

   

% Wc = gram(Gn,'c');


end

