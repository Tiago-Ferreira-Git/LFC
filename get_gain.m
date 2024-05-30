function [K] = get_gain(A,B,E,R_,weight_freq,W,n_NC)

    % r = size(A,1) - n_NC;
    % Z = size(A,1);


    ctrbM = ctrb(A,B); % controlability matrix
    r = rank(ctrbM);
    Z = size(A,1);
    H = orth(ctrbM);
    V = null(H');
    W = [H V];

    A_hat = W\A*W;
    Bg_hat = W\B; 
    Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
    A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];
    
    A11 = A_hat(1:r,1:r);
    A12 = A_hat(1:r,r+1:end);
    A21 = A_hat(r+1:end,1:r);
    A22 = A_hat(r+1:end,r+1:end);
    
    %% Controller gain synthesis 
    % Compute LQR weight matrices as detailled in [1] for DTUC
    Q = [eye(r) zeros(r,Z-r)]*W'*diag(weight_freq)*...
        W*[eye(r) ;zeros(Z-r,r)];
    R = R_*eye(size(B,2));
    
    K = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q,R,E,500e3,1e-5,A,B,r,Z,W);
end

