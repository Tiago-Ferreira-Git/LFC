function [A,B,W] = discrete_dynamics(A_global,B_global,W_global,h)
    
    
    G = expm([A_global W_global; zeros(size(W_global,2),size(A_global,2)) zeros(size(W_global,2),size(W_global,2))]*h);
    
    W = G(1:size(A_global,1),size(A_global,2)+1:end);
    
    G = expm([A_global B_global; zeros(size(B_global,2),size(A_global,2)) zeros(size(B_global,2),size(B_global,2))]*h);
    
    B = G(1:size(A_global,1),size(A_global,2)+1:end);

    A = G(1:size(A_global,1),1:size(A_global,2));
    

    % D_global = zeros(size(A,1),size(B_global,2));
    % 
    % rank(ctrb(A_global,B_global))
    % rank(ctrb(A,B))
    % 
    % sys = ss(A_global,B_global,eye(size(A,1)),D_global);
    % sys = c2d(sys,h);


end
