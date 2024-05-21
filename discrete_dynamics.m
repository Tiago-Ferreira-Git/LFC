function [A,B,C,W,L,E] = discrete_dynamics(A_global,B_global,C_global,D_global,W_global,E_global,bus_ss,h,n_areas,L)
    
    
    G = expm([A_global W_global; zeros(size(W_global,2),size(A_global,2)) zeros(size(W_global,2),size(W_global,2))]*h);
    
    W_global_discrete = G(1:size(A_global,1),size(A_global,2)+1:end);
    
    
    sys = ss(A_global,B_global,C_global,D_global);
    sys = c2d(sys,h);
    
    A = zeros(size(A_global,1)+n_areas);
    
    B = zeros(size(B_global,1)+n_areas,size(B_global,2));
    %
    C = zeros(size(C_global,1),size(C_global,2)+n_areas);
    
        %
    W = zeros(size(W_global_discrete,1)+n_areas,size(W_global_discrete,2));
    
    E = zeros(size(E_global,1),size(E_global,2)+n_areas);
    
    
    A(1:size(A_global,1),1:size(A_global,1)) = sys.A;
    
    A(size(A_global,1)+1:end,size(A_global,1)+1:end) = eye(n_areas);


    %add error integrator

    for i=1:n_areas
        A(size(A_global,1)+i,sum(bus_ss(1:i-1,2))+1) = -h;
    end
    
    
    
    B(1:size(B_global,1),1:size(B_global,2)) = sys.B;
    
    C(1:size(C_global,1),1:size(C_global,2)) = sys.C;
    
    W(1:size(W_global,1),1:size(W_global,2))  = W_global_discrete;
    
    E(1:size(E_global,1),1:size(E_global,2)) = E_global;
    
    L(1:1+size(L,1):end) = 1;
    E(1:size(E_global,1),size(E_global,2)+1:end) = L;

end

