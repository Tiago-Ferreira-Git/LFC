function dxdt = nonlinear_model(t,x,K,network,bus_ss,x0,u0,w,A,B,W,delta_u)
   

    k = find(t >= 0:3600:24*3600);
    k = k(end);

    dxdt = zeros(size(x0));

    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];

    
    
    x = x - x0;

    if nargin < 12
        u = -K*x;
    else
        u = delta_u;
    end
    u = min(max(u,-0.1),0.1);
    dxdt = A*x + B*u + W*w(:,k);
    for i = 1:length(network)
        
        x_ = x(freq_index(i):freq_index(i+1)-1);
        %dxdt(freq_index(i):freq_index(i+1)-1,1) = network(i).A*x_ + network(i).B*u(u_index(i):u_index(i+1)-1) + network(i).W*w(w_index(i):w_index(i+1)-1,k);

        for j = 1:size(network(i).to_bus,1)

            ang_diff = 377*(x_(end) - x(angle_index(network(i).to_bus(j,1))));

            dxdt(freq_index(i)) = dxdt(freq_index(i)) + network(i).W(1,1)*sin(ang_diff)/(network(i).to_bus(j,5));
        end
    
    end



    dxdt =  dxdt;
end

