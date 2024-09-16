function dxdt = nonlinear_model(t,x,K,network,bus_ss,x0,u0,PL,Pres,Pt0,u_index,delta_u,debug)
   

    delta_x = x - x0;

    if nargin < 12
        
        delta_u = -K*delta_x;
    end


    dxdt = zeros(size(x0));

    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];
    
    
    
    u = u0 + delta_u;
    
    index_res = 1;
    for i = 1:length(network)
        n_res = network(i).res;
        freq_feedback = zeros(size(network(i).A,1),1);
        freq_feedback(3:3:end) = -network(i).freq_feedback;
        area_index = freq_index(i):freq_index(i+1)-1;
        x_ = x(area_index);

        x_mec_index = freq_index(i)+1:freq_index(i+1)-2-n_res;
        u_mec_index = u_index(i):u_index(i+1)-1;

        dxdt(x_mec_index) = network(i).A*x(x_mec_index) + network(i).B*u(u_mec_index)+ freq_feedback.*(delta_x(area_index(1)))  ;
        
        
        P_mech = network(i).C_mech*x_;
        
        % Frequency dynamics 
        if n_res ~= 0
            P_res = network(i).C_res*x_;

            error 'P_res(i,k)?' 
            dxdt(freq_index(i)) = (P_mech./x_(1) + P_res - PL(i) - Pt0(i))./network(i).inertia;


            dxdt(freq_index(i+1)-n_res:freq_index(i+1)-1) = network(i).W_res*Pres(index_res:index_res+n_res-1,k);
            index_res = index_res + n_res;
            
        else
            dxdt(freq_index(i)) = ( (P_mech./x_(1))  - PL(i) - Pt0(i) )./network(i).inertia;
        end


        for j = 1:size(network(i).to_bus,1)
            
            if debug == 1
                ang_diff = 2*pi*60*(delta_x(angle_index(i)) - delta_x(angle_index(network(i).to_bus(j,1))));
            else
                ang_diff = (delta_x(angle_index(i)) - delta_x(angle_index(network(i).to_bus(j,1))));
            end

            dxdt(freq_index(i)) = dxdt(freq_index(i)) + (sin(ang_diff)/(network(i).to_bus(j,5)))/network(i).inertia;
        end

        % Frequency error dynamics
        if debug == 1
            dxdt(angle_index(i)) = -delta_x(freq_index(i));
        else
            dxdt(angle_index(i)) = -2*pi*60*delta_x(freq_index(i));
        end

        
      
    end
end

