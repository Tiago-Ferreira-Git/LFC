function [A_global,bus,area_gain_tie_lines] = update_dynamics(areas,bus,bus_initial,line,n_areas,A_global,bus_ss,g,network,hour)

    ren_data = load('data\solar.mat');
    ren_data = ren_data.data;
    n_ren = size(ren_data.bus,1);

    %bus = bus_initial;
    %select the first bus for each area

    buses_to_increase_con = [];
    %zeros(n_areas,1);
    for i = 1:n_areas
        mask = find(areas == i);
        buses_to_increase_con = [buses_to_increase_con  ; mask(1:min(3,size(mask,1)))];

    end
    bus_prev = bus;
    

    x = 0:3600:3600*24;
    %gaussmf(x,[10*3600 12*3600])
    %nominal_load_increase_profile = gaussmf(x,[4*3600 12*3600])*2;
    nominal_load_p_increase_profile = gaussmf(x,[4*3600 12*3600])*0.6+0.1;
    nominal_load_q_increase_profile = gaussmf(x,[4*3600 12*3600])*0.8+0.1;

    bus(buses_to_increase_con,6) = nominal_load_p_increase_profile(hour);
    bus(buses_to_increase_con,7) = nominal_load_q_increase_profile(hour);
    
    buses_with_generation = bitor(bus(:,10)==2 , bus(:,10)==1 );
    bus(buses_with_generation,11) = 9999;
    bus(buses_with_generation,12) = -9999;

    
    bus(end-n_ren+1:end,4) = ren_data.data(hour,2:end)';
    bus(119:end,11) = min(1,ren_data.data(hour,2:end)'./10)+0.001;
    bus(119:end,12) = -min(1,ren_data.data(hour,2:end)'./10)-0.001;
    bus(119:end,10) = 2;

    %ren_data.data(hour,2:end)'
    size(bus(buses_to_increase_con,1),1)
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

    area_gain_tie_lines = zeros(1,n_areas);
    for i = 1:size(network,2)
        
        neighbours = network(i).to_bus;

         T_i = 0;
        for j = 1:size(neighbours,1)
            T_ji = 377*cos(bus(g.bus.bus_int(neighbours(j,3)),3) - bus(g.bus.bus_int(neighbours(j,2)),3))/(neighbours(j,5));

            A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:neighbours(j,1)-1,2))+1 ) =  A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:neighbours(j,1)-1,2))+1 ) -T_ji;
              
            T_i =  T_i + T_ji;
        end
        area_gain_tie_lines(i) = T_i;
        A_global(sum(bus_ss(1:i,2)), sum(bus_ss(1:i-1,2))+1 ) = T_i;
    
    end


    bus(:,3) = rad2deg(bus(:,3));
    bus = bus(:,1:13);
    
end




