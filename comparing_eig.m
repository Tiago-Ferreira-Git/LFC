close all;
clear all

load('data/sim_118_30')


h = 0.1;

[A,B,W] = discrete_dynamics(A_c,B_c,W_c,h);

[K_method,E_fs] = slow_ss(mpc,debug,network,h);


C_ = zeros(n_areas*2,size(A_c,1));
C_(1:2:end) = C(1:3:end,:);
C_(2:2:end) = -C(3:3:end,:);


eigen_discrete = eig(A - B*K_method*C_);
eigen_cont_method = log(eigen_discrete)/h;


%%Centralized
figure
plot(real(eigen_cont_method), imag(eigen_cont_method),'x')
xlabel("Real axis")
ylabel("Imaginary axis")
savefig('./fig/poles_my_method.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/poless_my_method.png','png');


q = zeros(1,size(A,1));

freq_index = [1 ;cumsum(bus_ss(1:end-1,2))+1];
angle_index = cumsum(bus_ss(:,2));
q(1,freq_index) = 4;
q(1,angle_index) =  1;

R_ = 10000;



[K_cent,S] = dlqr(A,B,diag(q),R_*eye(size(B,2)));
eigen_discrete = eig(A - B*K_cent);
eigen_cont_centralized = log(eigen_discrete)/h;



figure
plot(real(eigen_cont_centralized), imag(eigen_cont_centralized),'x')
xlabel("Real axis")
ylabel("Imaginary axis")
savefig('./fig/poless_centralized.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/poless_centralized.png','png');



%%Decentralized
load('D:\OneDrive - Universidade de Lisboa\Aulas\Tese\Implementation\RES in PST\analysis\area partition\data\testing\K_0.01_0.10_100000.mat')
eigen_discrete = eig(A - B*K);
eigen_cont_decentralized = log(eigen_discrete)/h;



figure
plot(real(eigen_cont_decentralized), imag(eigen_cont_decentralized),'x')
xlabel("Real axis")
ylabel("Imaginary axis")
savefig('./fig/poless_decentralized.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/poless_decentralized.png','png');

%%
figure
hold on
plot(real(eigen_cont_method),imag(eigen_cont_method),'x');
plot(real(eigen_cont_centralized),imag(eigen_cont_centralized),'x');
xlim([-0.03 0])
xlabel("Real axis")
ylabel("Imaginary axis")
legend({'New Method','Centralized'},'Interpreter','latex','Location','best')
savefig('./fig/closed_loop_zoomed_poles.fig');
set(gcf,'renderer','Painters');
saveas(gca,'./fig/closed_loop_zoomed_poles.png','png');






