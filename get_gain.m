function [K] = get_gain(A,B,E,R_,weight_freq,W,r)
    Z = size(A,1);
    ctrbM = ctrb(A,B);
    if (nargin <= 5)
        warning 'Using uninformed transformation'
        ctrbM = ctrb(A,B); % controlability matrix
        r = rank(ctrbM)
        H = orth(ctrbM);
        V = null(H');
        W = [H V];
 
    end

    % max(A,[],'all')
    % min(A,[],'all')
    % max(ctrbM,[],'all')
    % min(ctrbM,[],'all')

    A_hat = W\A*W;
    Bg_hat = W\B; 
    Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
    A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];
    

    
    
    A11 = A_hat(1:r,1:r);
    A12 = A_hat(1:r,r+1:end);
    A21 = A_hat(r+1:end,1:r);

        
    % mask = abs(A21) < 1e-4;
    % A21(mask) = 0;

    A22 = A_hat(r+1:end,r+1:end);
    
    %% Controller gain synthesis 
    % Compute LQR weight matrices as detailled in [1] for DTUC
    Q = [eye(r) zeros(r,Z-r)]*W'*diag(weight_freq)*...
        W*[eye(r) ;zeros(Z-r,r)];
    R = R_*eye(size(B,2));

    K = LQROneStepLTI_augmented(A1_hat,Bg1_hat,Q,R,E,300e3,1e-5,A,B,r,Z,W);
end

