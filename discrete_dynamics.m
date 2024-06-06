function [A,B,W] = discrete_dynamics(A,B,W,h)
    
    
    G = expm([A W; zeros(size(W,2),size(A,2)) zeros(size(W,2),size(W,2))]*h);
    
    W = G(1:size(A,1),size(A,2)+1:end);
    
    G = expm([A B; zeros(size(B,2),size(A,2)) zeros(size(B,2),size(B,2))]*h);
    
    B = G(1:size(A,1),size(A,2)+1:end);

    A = G(1:size(A,1),1:size(A,2));
    

    % D_global = zeros(size(A,1),size(B_global,2));
    % 
    % rank(ctrb(A_global,B_global))
    % rank(ctrb(A,B))
    % 
    % sys = ss(A_global,B_global,eye(size(A,1)),D_global);
    % sys = c2d(sys,h);


end
