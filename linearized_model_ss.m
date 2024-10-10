function dxdt = linearized_model_ss(t,x,bus_ss,A,B,W,x0,u0,delta_PL,delta_u,K)
   

    %delta_u = -K*(x);
     
    delta_x = x - x0;
    

    dxdt = A*delta_x + B*delta_u + W*delta_PL;      
    
end

