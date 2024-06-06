function [] = canonical_structure(A,B,n_areas)
    ctrbM = ctrb(A,B); % controlability matrix

    r = rank(ctrbM)
    %max(max(ctrbM,[],2),[],1)
    Z = size(A,1)
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
end

