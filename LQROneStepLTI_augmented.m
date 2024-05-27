%% LQROneStepLTI_augmented - Description
% This function computes the steady-state augmented one-step LQR regulator 
% gain for a window w. Method derived in [1].
% Input:    - A_hat, B_hat
%           - Q, R
%           - E: sparsity pattern
%           - itMax: maximum number of iterations until convergence 
%           - epslInf: minimum relative improvement
%           - A, B
%           - r,Z,W: as defined in [1]
% Output:   - K: nxo steady-state gain matrix
%           - P: nxn steady-state estimation error covariance matrix
% Important notes: 
%           - output gain corresponds to the control law: u(k)=-K(k)*x(k)
% WARNING: Returns Kinf = NaN and Pinf = NaN if convergence could not be reached 
function [K,P] = LQROneStepLTI_augmented(A_hat,B_hat,Q,R,E,itMax,epslInf,A,B,r,Z,W)
% Gain computation
n = size(E,2); % Get value of n from the size of A 
m = size(E,1); % Get value of n from the size of B  
P = Q; % terminal condition
Pprev = NaN;
it = itMax;
while it > 0 % LQ iterations
    K = zeros(m,n);
    S = R+B_hat'*(P)*B_hat;
    for i = 1:n
        L = zeros(n);
        L (i,i) = 1; % Generate matrix L_i
        M = zeros(m);
        for j = 1:m % Gererate matrix M_i
            if E(j,i) ~= 0
                M(j,j) = 1;
            end
        end
        % Compute the ith term of the summation 
        K = K + (eye(m)-M+M*S*M)\(M*(B_hat')*P*A_hat*([eye(r) zeros(r,Z-r)]/W)*L');
    end
    % Update P
    P_ =((W')\[eye(r);zeros(Z-r,r)]) *(P)* ([eye(r) zeros(r,Z-r)]/W);
    Q_ = ((W')\[eye(r);zeros(Z-r,r)]) *(Q)* ([eye(r) zeros(r,Z-r)]/W);
    P = Q_+K'*R*K+...
        (A-B*K)'*P_*(A-B*K);
    P = [eye(r) zeros(r,Z-r)]*W'*P*W*[eye(r);zeros(Z-r,r)];
    % Check convergence
    it = it-1;
    if abs(trace(P)-trace(Pprev))/trace(Pprev) < epslInf
        itMax - it
        break; 
    end
    Pprev = P;
    if it == 0
        fprintf("One-step did not converge.\n");
        P = NaN;
        K = NaN;
    end
end 
end
