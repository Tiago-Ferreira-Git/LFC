function [A] = teste(A,B,E,R_,weight_freq,n_areas,bus_ss)
    
    original_size = size(A,1);
    A_ = zeros(original_size+n_areas);
    A_(1:original_size,1:original_size) = A;
    for i=1:n_areas
        A_(original_size+i,sum(bus_ss(1:i-1,2))+1) = -1;
    end

    A = A_;
    B = [B ; zeros(n_areas,size(B,2))];
    ctrbM = ctrb(A,B); % controlability matrix
    
    

    
    Z = size(A,1);
    H = orth(ctrbM);
    V = null(H');
    W = [H V];

    A_hat = W\A*W;
    Bg_hat = W\B; 
    Bg1_hat = [eye(r) zeros(r,Z-r)]*Bg_hat;
    A1_hat = [eye(r) zeros(r,Z-r)]*A_hat*[eye(r);zeros(Z-r,r)];



    A11 = A_hat(1:r,1:Z);
    A12 = A_hat(1:r,Z:end);
    A21 = A_hat(r:end,1:Z);
    A22 = A_hat(r:end,Z:end)
end

