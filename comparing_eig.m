close all;
clear all

load('data/sim_118_30')



% h = 0.02;
% C_ = eye(size(A,1));
% C_ = C;
% sys = ss(A_c,B_c,C_,zeros(size(C_,1),size(B,2)));
% 
% sys_d = ss(A,B,C_,zeros(size(C_,1),size(B,2)),h);
% 
% 
% figure('Position',4*[0 0 144 132])
% pzmap(sys)
% saveas(gca,'./fig/pzmap_continuous_open_loop.png','png');
% savefig('./fig/pzmap_continuous_open_loop.fig');
% 
% 
% figure('Position',4*[0 0 144 132])
% pzmap(sys_d)
% saveas(gca,'./fig/pzmap_discrete_open_loop.png','png');
% savefig('./fig/pzmap_discrete_open_loop.fig');
% 
% feedback_ss_discrete = ss(A-B*K,zeros(size(A_c,1),1),C_,zeros(size(C_,1),1),h);
% figure('Position',4*[0 0 144 132])
% pzmap(feedback_ss_discrete)
% saveas(gca,'./fig/pzmap_discrete_closed_loop.png','png');
% savefig('./fig/pzmap_discrete_closed_loop.fig');
% 
% feedback_ss_continuous = d2c(feedback_ss_discrete);
% figure('Position',4*[0 0 144 132])
% pzmap(feedback_ss_continuous)
% saveas(gca,'./fig/pzmap_continu_closed_loop.png','png');
% savefig('./fig/pzmap_continuous_closed_loop.fig');
% 
% 

% s = tf('s');
% G = C*inv(s.*eye(size(A)) - A )*B;

C_ = eye(size(A_c,1));
%C_ = C;


s = tf('s');



%

% ss = [1 ; cumsum(bus_ss(:,2))+1];
% u_ss = [1 ; cumsum(bus_ss(:,3))+1];
% y_ss = 1:3:(n_areas+1)*3;
% for k = 1:n_areas
%     A_area = A_c(ss(k):ss(k+1)-1,ss(k):ss(k+1)-1);
%     B_area = B_c(ss(k):ss(k+1)-1,u_ss(k):u_ss(k+1)-1);
%     C_area = C(y_ss(k):y_ss(k)-1 , ss(k):ss(k+1)-1);
% 
% 
%     % C_area = 1; 
%     % for o = 1:network(k).machines
%     %     C_area = [C_area zeros(size(C_area,1),3) ; zeros(1,size(C_area,2)) [1 1.25 0]];
%     % end
%     % C_area = [C_area zeros(size(C_area,1),1) ; zeros(1,size(C_area,2)) -1];
%     transfer_function = C_area/(s.*eye(size(A_area))-A_area)*B_area;
%     k;
%     for i = 1:size(C_area,1)
%         for j = 1:size(transfer_function.Numerator,2)
%             if any(real(roots(transfer_function.Numerator{i,j})) > 0.001)
%                 roots(transfer_function.Numerator{i,j})
% 
%                 f = figure;
%                 f.Position(3:4) = f.Position(3:4)*1.728;
%                 f.Position(2) = f.Position(2)-200;
%                 hold on
%                 plot(real(roots(transfer_function.Numerator{i,j})),imag(roots(transfer_function.Numerator{i,j})),"x")
%                 axis equal
%                 grid on
%                 xlabel("Re(z)")
%                 ylabel("Im(z)")
%                 if i == 1
%                     title(sprintf('Area %d , Frequency Machine Input %d',k,j))
%                 elseif i == size(C_area,1)                    
%                     title(sprintf('Area %d , Error Machine Input %d',k,j))
%                 else
%                     title(sprintf('Area %d , Machine Output %d Machine Input %d',k,i-1,j))
% 
%                 end
% 
% 
% 
%             end    
% 
% 
%         end
% 
%         i;
%     end
% end
% Controller gain synthesis 
q = zeros(1,size(A_c,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
angle_index = cumsum(bus_ss(:,2));

q(1,freq_index) = 40;    
q(1,angle_index) = 1;


for h = 10.^(-1:1:-1)
    [A,B,~] = discrete_dynamics(A_c,B_c,W_c,h);
    for R_ = 10.^(-2:1:3)
        disp('R:')
        disp(R_)
        load(sprintf('data/K_%.3f_%.2f.mat',R_,h));
        % K_neighbour = zeros(size(K));
        % K_neighbour(~logical(E_fs)) = K(~logical(E_fs));
        % K_local = zeros(size(K));
        % K_local(logical(E_fs)) = K(logical(E_fs));


        %sys = ss(A_c,B_c,C_,zeros(size(C_,1),size(B,2)));

        %sys_d = ss(A,B,C_,zeros(size(C_,1),size(B,2)),h);
        %feedback_ss_discrete = ss(A-B*K,zeros(size(A_c,1),1),C_,zeros(size(C_,1),1),h);
        %feedback_ss_continuous = d2c(feedback_ss_discrete);

        %figure('Position',4*[0 0 144 132])
        %pzmap(feedback_ss_continuous)
        %saveas(gca,'./fig/pzmap_continu_closed_loop.png','png');
        %savefig('./fig/pzmap_continuous_closed_loop.fig');

        eigen_discrete = eig(A-B*K);
        eigen_cont = log(eigen_discrete)/h;
        K_dec = K;

        figure
        hold on
        plot(real(eigen_cont),imag(eigen_cont),"x")
        axis equal
        grid on
        xlabel("Re(z)")
        ylabel("Im(z)")
        %title(sprintf('Decentralized , R = %.3f',R_))


        K_cen = dlqr(A,B,diag(q),0.1*eye(size(B,2)));
        eigen_discrete = eig(A-B*K_cen);
        eigen_cont = log(eigen_discrete)/h;


        plot(real(eigen_cont),imag(eigen_cont),"x")
        %title(sprintf('Centralized , R = %.3f',R_))
        legend({'Decentralized','Centralized'},'Location','best')
        title=sprintf('./fig/email/poles_h_%.2f_R_%.2f.fig',h,R_);
        savefig(title);
        set(gcf,'renderer','Painters');
        title=sprintf('./fig/email/poles_h_%.2f_R_%.2f.png',h,R_);
        saveas(gca,title,'png');


        
        eigen_discrete = eig(A);
        eigen_cont = eig(A_c);
        figure
        hold on
        plot(real(eigen_cont),imag(eigen_cont),"x")
        axis equal
        grid on
        xlabel("Re(z)")
        ylabel("Im(z)")
        title=sprintf('./fig/email/poles_openloop.png',h,R_);
        saveas(gca,title,'png');




    end
end

% C_area = [];
% for k = 1:n_areas
% 
%     C_area = [C_area zeros(size(C_area,1),1) ;  zeros(1,size(C_area,2)) 1]; 
%     for i = 1:network(k).machines
%         C_area = [C_area zeros(size(C_area,1),3) ; zeros(1,size(C_area,2)) [1 1.25 0]];
%     end
%     C_area = [C_area zeros(size(C_area,1),1) ; zeros(1,size(C_area,2)) -1];
% 
% 
% end
% 
% 
% 
% s = tf('s');
% a = inv(s.*eye(size(A_c)) - A_c);
% 
% %%
% disp('Full state feedback')
% transfer_function_fullstatefeedback = eye(size(A_c))*a*B_c;
% 
% disp('Measurements')
% transfer_function_meas = C_area*a*B_c;
% 
% disp('Measurements of output power per area')
% transfer_function_meas_area = C*a*B_c;


%%
clearvars -except bus_ss

load('D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Implementation\RES in PST\analysis\area partition\data\zeros.mat')

transfer_function = transfer_function_fullstatefeedback;

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
angle_index = cumsum(bus_ss(:,2));

o = 1;
for i = 1:size(transfer_function.Numerator,1)

    mask = i < freq_index(o+1) || i >= freq_index(o);
    if sum(i+1 < freq_index(o+1) || i+1 >= freq_index(o)) == 0
        o = o +1 ;
    end
    mask = find(mask == 1 );
    for j = 1:size(transfer_function.Numerator,2)
        if any(real(roots(transfer_function.Numerator{i,j})) > 0.001)
            roots(transfer_function.Numerator{i,j});

            f = figure;
            f.Position(3:4) = f.Position(3:4)*1.728;
            f.Position(2) = f.Position(2)-200;
            hold on
            plot(real(roots(transfer_function.Numerator{i,j})),imag(roots(transfer_function.Numerator{i,j})),"o")
            plot(real(roots(transfer_function.Denominator{i,j})),imag(roots(transfer_function.Denominator{i,j})),"x")
            axis equal
            grid on
            xlabel("Re(z)")
            ylabel("Im(z)")
            if ismember(i,freq_index)
                k = find(i == freq_index);
                title(sprintf('Area %d  Zeros $\\frac{\\omega_{%d}}{u_{%d}}$',k,k,j),'Interpreter','latex','FontSize',35)
            elseif ismember(i,angle_index)                  
                title(sprintf('Error frequency Input %d',j))
            else
                title(sprintf('Machine Output %d Machine Input %d',i,j))

            end



        end    


    end
end





