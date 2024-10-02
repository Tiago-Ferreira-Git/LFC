function dxdt = linearized_model_discrete(t,x,bus_ss,A,B,W,x0,u0,delta_PL,delta_u,K)
   

    delta_u = -K*(x-x0);
    
    angle_index = cumsum(bus_ss(:,2));
    freq_index = [1 ; angle_index+1];
    
 
    
    % x(angle_index) = mod(x(angle_index) + pi, 2*pi) - pi;
    delta_x = x - x0;
    

    dxdt = A*delta_x + B*delta_u + W*delta_PL;      
    
end

