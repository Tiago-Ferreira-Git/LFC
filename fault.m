function [A,B,C,W,E,network,g,bus_ss] = fault(g,bus,network,flag_ren,mac_bus_fault,tie_fault,h)


    if ~isempty(mac_bus_fault)
        for i = 1:size(mac_bus_fault,2)
            %clear variables in the network structure
            for j = 1:size(network,2)
                mask = network(j).mac_bus == mac_bus_fault(i);
                if(any(mask))
                    network(j).tg_con(mask,:) = zeros(size(network(j).tg_con(1,:)));
                end
            end
        
        end
        
        %clear active and reactive power in the bus array
        bus(mac_bus_fault,4) = 0;
        bus(mac_bus_fault,5) = 0;
        bus(mac_bus_fault,10) = 3;

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
        
    end

    mask = ismember(g.line(:,1:2), tie_fault,'rows');
    g.line(mask,:) = [];
    mask = ismember(g.line(:,1:2), tie_fault(:,[2,1]),'rows');
    g.line(mask,:) = [];
    


    line = g.line;
    [bus,~,~] = loadflow(bus,line,1e-8,50,1.0,'n',2);
        

    n_areas = size(network,2);
    [A,B,C,D,W,~,E,~,network,bus_ss,~] = get_global_ss(g,bus,n_areas,flag_ren,1,network);
    

    [A,B,W] = discrete_dynamics(A,B,W,h);


      


end

