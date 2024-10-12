function k_ties = compute_tielines(mpc,network)


    bus_sol = mpc.bus(:,9);
    bus_sol = deg2rad(bus_sol);


    k_ties = zeros(length(network),1);

    for i = 1:size(network,2)

        neighbours = network(i).to_bus;
        T_i = 0;
       
        for j = 1:size(neighbours,1)
            T_ji = cos(bus_sol(neighbours(j,3)) - bus_sol(neighbours(j,2)))/(neighbours(j,5));
            T_i =  T_i + T_ji;
        end

        k_ties(i) = T_i;
    end

end

