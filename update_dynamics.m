function [K] = update_dynamics(mpc,network,flag_ren,h,n_areas,simulation_seconds,k,Q,R,E,flag)

    % w_load__hour = load_('data\w_load_.mat');
    % w_load__hour = w_load__hour.w_load_(:,1:24);
    % 
    % w_res_hour = load_('data\w_res.mat');
    % w_res_hour = w_res_hour.w_res(:,1:24);

 

    load('data\load_profile.mat');

    if (size(mpc.bus,1) == 118)
        load('data\res_profile_118.mat');
    else
        load('data\res_profile.mat');
        max_res = 200;
        res.forecast = res.forecast*max_res;
        res.measured = res.measured*max_res;
        
    end




    n_res = sum(mpc.isolar_mask);

    % n_res = sum(bus_ss(:,4));
    %if the number of areas in the data available is lower than the areas 
    if size(load_.forecast,1) < length(mpc.bus) 
        load_.forecast = [repmat(load_.forecast, floor(length(mpc.bus)/size(load_.forecast,1)), 1); ...
                        load_.forecast(1:mod(length(mpc.bus), size(load_.forecast,1)),:)];
        load_.measured = [repmat(load_.measured, floor(length(mpc.bus)/size(load_.measured,1)), 1); ...
                load_.measured(1:mod(length(mpc.bus), size(load_.measured,1)),:)];
    else 
        load_.forecast = load_.forecast(1:length(mpc.bus),:);
        load_.measured = load_.measured(1:length(mpc.bus),:);
    end

    if size(res.forecast,1) < n_res 
        res.forecast = [repmat(res.forecast, floor(n_res/size(res.forecast,1)), 1); ...
                        res.forecast(1:mod(n_res, size(res.forecast,1)),:)];
        res.measured = [repmat(res.measured, floor(n_res/size(res.measured,1)), 1); ...
                res.measured(1:mod(n_res, size(res.measured,1)),:)];
    else 
        res.forecast = res.forecast(1:n_res,:);
        res.measured = res.measured(1:n_res,:);
    end



    simulation_hours = ceil(simulation_seconds/3600);


    n_machines = sum(mpc.isgen);

    

    %Test if all load and renewables converge
    mpc_initial = mpc;
    for i = 1:n_areas

        %Change loads
        %mpc_initial.bus(network(i).bus,3).*load_.forecast(network(i).bus,k);
        mpc.bus(network(i).bus,3) = mpc_initial.bus(network(i).bus,3).*load_.forecast(network(i).bus,k);
        load_measured(network(i).bus) = mpc_initial.bus(network(i).bus,3).*load_.measured(network(i).bus,k);     
        

        idx = network(i).res_nr;
        mpc.gen(idx,2) = res.forecast(idx-n_machines,k);

        %Active power limits
        mpc.gen(idx,9) = mpc.gen(idx,2) + 0.000001 ;
        mpc.gen(idx,10) = mpc.gen(idx,2) - 0.000001;
    
    
        %Reactive Power Limitis
        mpc.gen(idx,4) = min(1,res.forecast(idx-n_machines,k)./10)+0.001;
        mpc.gen(idx,4) = mpc.gen(idx,4);
    
        mpc.gen(idx,5) = 0;
        mpc.gen(idx,5) = mpc.gen(idx,5);
        
    end     
    [mpc,flag] = runopf(mpc,mpoption('verbose',0,'out.all',0));
    if ~flag
        error 'Load or Renewables profile made runopf not converge!'
    end
    debug = 2;

    [A_c,B_c,~,~,W_c,~,~,~,~,~,~,~,~,~,~] = get_global_ss(mpc,n_areas,flag_ren,debug,network);
    [A,B,~] = discrete_dynamics(A_c,B_c,W_c,h);

    if flag  
        K = dlqr(A,B,Q,R);
    else
        [K,~,~]  = LQROneStepLTI(A,B,Q,R,E,NaN);
    end

        
end
