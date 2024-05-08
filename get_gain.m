function [K] = get_gain(A_global,B_global,E,R_,weight_freq)
   
    ctrbM = ctrb(A_global,B_global); % controlability matrix
    
    r = rank(ctrbM);
    Z = size(A_global,1);
    H = orth(ctrbM);
    V = null(H');
    W = [H V];
    A_hat = W\A_global*W;
    Bg_hat = W\B_global; 
    Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
    A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];
    
    
    %% Controller gain synthesis 
    % Compute LQR weight matrices as detailled in [1] for DTUC
    Q = [eye(r) zeros(r,Z-r)]*W'*diag(weight_freq)*...
        W*[eye(r) ;zeros(Z-r,r)];
    R = R_*eye(size(B_global,2));
    
    K = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q,R,E,3e3,1e-5,A_global,B_global,r,Z,W);
end

