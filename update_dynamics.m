function [bus,area_gain_tie_lines,angles,q_gen,p_gen,q_load,p_load] = update_dynamics(areas,bus,line,n_areas,A_global,bus_ss,g,network,hour,bus_initial)

    ren_data = load('data\solar.mat');
    ren_data = ren_data.data;
    n_ren = size(ren_data.bus,1);

    %bus = bus_initial;
    %select the first bus for each area

    buses_to_increase_con = [];
    %zeros(n_areas,1);
    for i = 1:n_areas
        mask = find(areas == i);
        buses_to_increase_con = [buses_to_increase_con  ; mask(1:min(2,size(mask,1)))];

    end
    bus_prev = bus;
    

    x = 0:3600:3600*24;

    % bus(:,4) = bus(:,4)./2;
    % bus(:,5) = bus(:,5)./2;

    %increase loads
    %nominal_load_p_increase_profile = (gaussmf(x,[4*3600 12*3600])*0.4+0.1);
    %nominal_load_q_increase_profile = (gaussmf(x,[4*3600 12*3600])*0.6+0.1);

    %bus(buses_to_increase_con,6) = nominal_load_p_increase_profile(hour);
    %bus(buses_to_increase_con,7) = nominal_load_q_increase_profile(hour);
    
    %set gen
    bus(end-n_ren+1:end,4) = ren_data.data(hour,2:end)';
    bus(119:end,10) = 2;




    buses_with_generation = bitor(bus_initial(:,10)==2 , bus_initial(:,10)==1 );
    bus(buses_with_generation,11) = 3;
    bus(buses_with_generation,12) = -3;

    bus(119:end,11) = min(1,ren_data.data(hour,2:end)'./10)+0.001;
    bus(119:end,12) = -min(1,ren_data.data(hour,2:end)'./10)-0.001;


    %ren_data.data(hour,2:end)'
    
    tol = 1e-8;                       % tolerance for convergence
    itermax = 100;                    % maximum number of iterations
    acc = 1.0;                        % acceleration factor
    [bus,~,~] = loadflow(bus,line,tol,itermax,acc,'n',2);
    teste = bus(:,3) - bus_prev(:,3);
        
    bus(:,3) = deg2rad(bus(:,3));
    
    
    %Add tie-lines
    %A_global(6:6:end,:) = 0; 
    for i=1:n_areas
        A_global(sum(bus_ss(1:i,2)),:) = 0;
    end


    
    bus_sol = bus;
    bus_sol(:,3) = deg2rad(bus(:,3));


    %Add tie-lines

    index_ss = cumsum(bus_ss(:,2))+1;
    index_ss = [1 ; index_ss];

    area_gain_tie_lines = zeros(1,n_areas);

    for i = 1:size(network,2)


        neighbours = network(i).to_bus;
        T_i = 0;
        for j = 1:size(neighbours,1)
            T_ji = 377*cos(bus_sol(g.bus.bus_int(neighbours(j,3)),3) - bus_sol(g.bus.bus_int(neighbours(j,2)),3))/(neighbours(j,5));

            
            %Changin ptie for error integral
            A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) =  A_global(index_ss(i), index_ss(neighbours(j,1)+1)-1) -T_ji*bus_ss(i,5);
            

            T_i =  T_i + T_ji;
        end
        area_gain_tie_lines(i) = T_i;
        
        A_global(index_ss(i), index_ss(i+1)-1) = T_i*bus_ss(i,5);

    end


    

    bus(119:end,10) = 2;
    bus(:,3) = rad2deg(bus(:,3));
    bus = bus(:,1:13);
    angles = bus(:,3);
    p_gen = bus(buses_with_generation,4);
    q_gen = bus(buses_with_generation,5);

    p_load = bus(:,6);
    q_load = bus(:,7);
    
    size(bus(buses_to_increase_con,1),1)


end




