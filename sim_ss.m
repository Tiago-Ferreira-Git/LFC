function [y,u] = sim_ss(A,B,C,D,K,W,ref,d_input,t,x0)
    %sim_ss function to simulate the discrete ss
    %   
    %x_d = zeros(size(A,1),size(t,2));
    x = zeros(size(A,1),size(t,2));
    u = zeros(size(B,2),size(t,2));
    y = zeros(size(C,1),size(t,2));
    
    x(:,1) = x0;
    for k = 1:length(t)
        
        u(:,k) = ref(:,k)-K*x(:,k);
        x(:,k+1) = A*x(:,k) + B*u(:,k) + W*d_input(:,k);
        %x(:,k+1) = x(:,k) + h*x_d(:,k);
        y(:,k) = C*x(:,k);
    end
    y = y';
end

