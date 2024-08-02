for h = 10.^(-1:1:0)
    disp('h:')
    disp(h)
    for R_ = 10.^(-2:1:3)
        disp('R')
        disp(R_)
        clearvars -except path h R_; close all; clc;
        debug = 2;
        
        load('data/sim_118_30')
        % load('data/sim_14_3')
        % disp('Sampling time should be:')
        % pi/abs(min(real(eig(A_c))))
        % disp('Sampling time should be:')
        % pi/max(abs(eig(A_c)))
        
        
        %h = 0.2;
        
        [A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);
        
        
        
        % Controller gain synthesis 
        q = zeros(1,size(A,1));
        
        freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
        angle_index = cumsum(bus_ss(:,2));
        q(1,freq_index) = 40;    
        q(1,angle_index) = 1;
        if debug == 1
        
            q(1,freq_index) = 40;    
            q(1,angle_index) = 1;
        else
            
            q(1,freq_index) = 40;    
            q(1,angle_index) = 0.000005;
        
        end
        
       
        K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E);
        
        save(sprintf('data/K_%.3f_%.2f.mat',R_,h),"K");
    
    
    
    end
end