function [A,B,W,L] = discrete_dynamics(A_global,B_global,W_global,h,L)
    
    
    G = expm([A_global W_global; zeros(size(W_global,2),size(A_global,2)) zeros(size(W_global,2),size(W_global,2))]*h);
    
    W = G(1:size(A_global,1),size(A_global,2)+1:end);
    
    G = expm([A_global B_global; zeros(size(B_global,2),size(A_global,2)) zeros(size(B_global,2),size(B_global,2))]*h);
    
    B = G(1:size(A_global,1),size(A_global,2)+1:end);

    A = G(1:size(A_global,1),1:size(A_global,2));
    
    % sys = ss(A_global,B_global,C_global,D_global);
    % sys = c2d(sys,h);
    

end

