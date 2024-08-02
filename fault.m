function [A,B,C,W,machine_ss,E,network,mpc,bus_ss,x] = fault(mpc,network,flag_ren,mac_fault,tie_fault,h,x,ss_dim)


    

    if ~isempty(mac_fault)

        %Cannot remove the reference machine
        if(any(mpc.bus(mpc.gen(mac_fault,1),2) == 3))
            mask = ismemeber(mac_fault, find(mpc.bus(:,2) == 3));
            mac_fault(mask) = [];
            warning 'User tried to apply fault to the reference machine, ignored this command'
        end

        %mac_fault is the generator number
        for i = 1:size(mac_fault,2)
            %clear variables in the network structure
            x_index = 1;
            for j = 1:size(network,2)
                
                mask = network(j).mac_nr == mac_fault(i);
                if(any(mask))
                    network(j).tg_con(mask,:) = zeros(size(network(j).tg_con(1,:)));
                    mac_index = 3* (find(network(j).mac_nr == mac_fault(i)) - 1);
                    mask = x_index + mac_index +1 : x_index + mac_index +1 + 2;
                    x(mask) = 0;
                end
                x_index = x_index + ss_dim(j);
            end
        
        end
        
        %clear active and reactive power in the bus array
        %mpc.bus(mpc.gen(mac_fault,1),10) = 1;
        mpc.gen(mac_fault,2) = 0;
        mpc.gen(mac_fault,3) = 0;
        mpc.gen(mac_fault,4) = 1e-6;
        mpc.gen(mac_fault,5) = -1e-6;
        mpc.gen(mac_fault,9) = 0;
        mpc.gen(mac_fault,10) = 0;
        %mpc.gencost(mac_fault,:) = [];
        

    end

    if ~isempty(tie_fault)
        for i = 1:size(tie_fault,2)  
            for j = 1:size(network,2)
                %clear variables in the network structure
                mask = ismember(network(j).to_bus(:,2:3), tie_fault,'rows');
                if(any(mask))
                    network(j).to_bus(mask,:) = [];
                end
                mask = ismember(network(j).to_bus(:,2:3), tie_fault(:,[2,1]),'rows');
                if(any(mask))
                    network(j).to_bus(mask,:) = [];
                end
            end
        
        end

        mask = ismember(mpc.branch(:,1:2), tie_fault,'rows');
        mpc.branch(mask,:) = [];
        mask = ismember(mpc.branch(:,1:2), tie_fault(:,[2,1]),'rows');
        mpc.branch(mask,:) = [];
        
    end


    
    [mpc,flag] = runopf(mpc,mpoption('verbose',0,'out.all',0));
    %,mpoption('verbose',0,'out.all',0)
    if ~flag
        error 'Power Flow did not converge. Try running again or changing faults'
    end

    n_areas = size(network,2);
    [A,B,C,~,W,machine_ss,~,~,E,~,network,bus_ss,~] = get_global_ss(mpc,n_areas,flag_ren,2,network);
    
    [A,B,W] = discrete_dynamics(A,B,W,h);


      


end

