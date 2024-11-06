function [A,B,W] = discrete_dynamics(A,B,W,h)
%discrete_dynamics   Effeciently discretize the continuous-time dynamics.
%
%   Inputs:                          
%       
%       A - Continuous-time transition matrix. 
%       B - Continuous-time control matrix.
%       W - Continuous-time disturbance matrix.
%       h - Sampling time in seconds.
%
%
%   Outputs:
%
%       A - Discrete-time transition matrix. 
%       B - Discrete-time control matrix.
%       W - Discrete-time disturbance matrix.
    

    
    G = expm([A W; zeros(size(W,2),size(A,2)) zeros(size(W,2),size(W,2))]*h);
    
    W = G(1:size(A,1),size(A,2)+1:end);


    G = expm([A B; zeros(size(B,2),size(A,2)) zeros(size(B,2))]*h);
    
    B = G(1:size(A,1),size(A,2)+1:end);

    A = G(1:size(A,1),1:size(A,2));
    

end
