load('data/sim_118_30')

h = 0.1;


E_to = zeros(size(E));

index_ss = cumsum(bus_ss(:,2))+1;
index_ss = [1 ; index_ss];
tg_size = cumsum([ 1 ; bus_ss(:,3)]);

for i = 1:size(network,2)


    neighbours = network(i).to_bus;

    E_to(tg_size(i):tg_size(i+1)-1,index_ss(i):index_ss(i+1)-1) = ones(bus_ss(i,3),bus_ss(i,2));
    
    
    unique_areas = unique(network(i).to_bus(:,1),'rows');


    for j = unique_areas'
        neighbours = [neighbours ; network(j).to_bus];
    end

    unique_areas = unique(neighbours(:,1),'rows');

    for j = unique_areas'
        neighbours = [neighbours ; network(j).to_bus];
    end

    for j = 1:size(neighbours,1)
        E_to(tg_size(i):tg_size(i+1)-1,index_ss(neighbours(j,1)):index_ss(neighbours(j,1)+1)-1) = ones(bus_ss(i,3),bus_ss(neighbours(j,1),2));
    end


end


int_values = [0.01,0.1];
R_values = [10, 100000];

for int_value = int_values
    %10.^(-1:1:2)
    disp('int_value:')
    disp(int_value)
    for R_ = R_values
        %10.^(1:1:6)
        disp('R')
        disp(R_)
        clearvars -except path h R_ E_to int_value int_values R_values;
        debug = 2;
        
        load('data/sim_118_30')
        %eigen_values = eig(A_c);
        % for i = eigen_values'
        %     %i
        %    if rank([A_c-i.*eye(size(A_c))  B_c]) ~= size(A_c,1)
        %         i
        %    end
        %end

        % load('data/sim_14_3')
        % disp('Sampling time should be:')
        % pi/abs(min(real(eig(A_c))))
        % disp('Sampling time should be:')
        % pi/max(abs(eig(A_c)))
        
        
        %h = 0.2;
        
        [A,B,~] = discrete_dynamics(A_c,B_c,W_c,h);
        
        
        
        % Controller gain synthesis 
        q = zeros(1,size(A,1));
        
        freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
        angle_index = cumsum(bus_ss(:,2));

     
        
        q(1,freq_index) = 4;    
        q(1,angle_index) = int_value;
        
       
        K  = LQROneStepLTI(A,B,diag(q),R_*eye(size(B,2)),E_to);

        save(sprintf('data/third_order/K_%.3f_%.2f_%.6f.mat',R_,h,int_value),"K");
    
    
    
    end
end