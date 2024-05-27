function [] = canonical_structure(A,B,W,n_NC)
    
    r = size(A,1) - n_NC;
    Z = size(A,1);
    A_hat = W\A*W;
    Bg_hat = W\B; 
    Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
    A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];


    A_hat = W\A*W;
    Bg_hat = W\B; 
    Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
    A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];


    r
    rank(ctrb(A1_hat,Bg1_hat))
    
    A11 = A_hat(1:r,1:r);
    A12 = A_hat(1:r,r+1:end);
    A21 = A_hat(r+1:end,1:r);
    A22 = A_hat(r+1:end,r+1:end);



    [eigen_values,eigen_vectors] = eig(A);
    diag_A = eigen_vectors\A*eigen_vectors;
    diag_B = eigen_vectors\B;

end

