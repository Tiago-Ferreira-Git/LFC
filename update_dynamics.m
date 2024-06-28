function [mpc,A,B,C,D,W,machine_ss,ktie] = update_dynamics(mpc,network,flag_ren,hour,mpc_initial,idx,h)

    data = load('data\solar.mat');
    data = data.data;


    bus_index = 1:4:size(mpc.bus,1);


    x = 0:24;
    
    nominal_load_p_increase_profile = (gaussmf(x,[4 12])*1);
    nominal_load_q_increase_profile = (gaussmf(x,[4 12])*1);
    
    if flag_ren
        %Active power
        mpc.gen(idx,2) = data.data(hour,2:end)';
        mpc.gen(idx,2) = mpc.gen(idx,2)*mpc_initial.baseMVA;
    
        %Active power limits
        mpc.gen(idx,9) = mpc.gen(idx,2) + 0.000001 ;
        mpc.gen(idx,10) = mpc.gen(idx,2) - 0.000001;
    
    
        %Reactive Power Limitis
        mpc.gen(idx,4) = min(1,data.data(hour,2:end)'./10)+0.001;
        mpc.gen(idx,4) = mpc.gen(idx,4)*mpc_initial.baseMVA;
    
        mpc.gen(idx,5) = 0;
        mpc.gen(idx,5) = mpc.gen(idx,5)*mpc_initial.baseMVA;
    end

    mpc.bus(bus_index, 3) = mpc_initial.bus(bus_index, 3) *(1 + nominal_load_p_increase_profile(mod(hour,24)+1));
    mpc.bus(bus_index, 4) = mpc_initial.bus(bus_index, 4) *(1 + nominal_load_q_increase_profile(mod(hour,24)+1));
    

    [mpc,flag] = runopf(mpc,mpoption('verbose',0,'out.all',0));
    if ~flag
        error 'Power Flow did not converge.'
    end

    n_areas = size(network,2);
    [A,B,C,D,W,machine_ss,~,~,~,~,~,bus_ss,~] = get_global_ss(mpc,n_areas,flag_ren,network);

    [A,B,W] = discrete_dynamics(A,B,W,h);
    

    %ktie with C
    ktie = zeros(1,size(network,2));
    rows = 3:4:size(C,1);
    cols = cumsum(bus_ss(:,2));
    for i = 1:size(network,2)
        ktie(i) = C(rows(i),cols(i));
    end
end




