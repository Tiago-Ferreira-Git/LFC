function [nominal,nominal_fault,w_m] = nominal_profiles(simulation_hours,n_gen,n_res,n_areas,mpc_initial,network,idx_initial,t_fault,mac_remove,bus_remove,flag_ren,h,dim_ss,bus_ss)

    
    nominal = zeros(simulation_hours,n_gen-n_res);
    nominal_fault = zeros(simulation_hours,n_gen-n_res);
    w_m = zeros(dim_ss,simulation_hours);

    
    mpc_ = mpc_initial;
    mpc_fault = mpc_initial;
    mask = true(1,n_gen-n_res);
    
    mac_idx = zeros(n_gen-n_res,1);
    idx = idx_initial;
    j = 1;
    for i = 1:n_areas
        mac_idx(j:j+network(i).machines-1,1) =  network(i).mac_nr;
        j = j + network(i).machines;
    end
    
    
    for i = 1:simulation_hours
        [mpc_,~,~,~,~,~,~,~] = update_dynamics(mpc_,network,flag_ren,i,mpc_initial,idx,h);
        nominal(i,:) = mpc_.gen(1:end-n_res,2);
    
    
        if i == t_fault
            mask = true(1,n_gen-n_res);
            mask(mac_remove) = 0;
            idx = idx_initial - size(mac_remove,2);
            [~,~,~,~,~,~,~,mpc_fault,~] = fault(mpc_fault,network,flag_ren,mac_remove,bus_remove,h);
        end
        [mpc_fault,~,~,~,~,~,~,~] = update_dynamics(mpc_fault,network,flag_ren,i,mpc_initial,idx,h);
        nominal_fault(i,:) = mpc_fault.gen(1:end-n_res,2);

        w_m(:,i) = initial_conditions(dim_ss,bus_ss,network,mpc_fault);
        %nominal_fault(i,~mask) = 0;

    end

    w_m = [zeros(dim_ss,1) diff(w_m,1,2)];
    
    figure
    hold on
    p3 = plot(1:simulation_hours,nominal./mpc_initial.baseMVA,'Color',[0.8500 0.3250 0.0980],'LineWidth',1.4);
    p1 = plot(1:simulation_hours,nominal_fault./mpc_initial.baseMVA,'Color',[0 0.4470 0.7410],'LineWidth',1.4);
    ylabel('$P_G$ (pu)','interpreter','latex');
    xlabel('$t \;[\mathrm{h}]$','Interpreter','latex');
    legend([p3(1) p1(1)],{'Nominal','Nominal Fault'})
    hold off
    set(gcf,'renderer','Painters');
    title='./fig/nominal.png';
    saveas(gca,title,'png');
    
    
    nominal_fault = nominal_fault(:,mac_idx);
    nominal = nominal(:,mac_idx);




end

