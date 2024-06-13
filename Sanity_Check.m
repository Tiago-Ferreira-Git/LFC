
load Area.mat



%Add tie-line
A(end,1) = 1e3;
A(1,end) = -1/network(1).inertia;
A12 = zeros(size(A));
A12(end,1) = -1e3;
A21 = A12;


%Global matrix (Output f1)
A_ =  [A A12; A21 A];
B_= [B zeros(size(B));zeros(size(B)) B ];
C_ = zeros(1,size(A_,1));
C_(1,1) = 1;


%Get transfer function
s = tf('s');
global_ = C_*inv(s*eye(size(A_)) - A_)*B_;


P1 = [1 0 0 0] *inv(s*eye(size(A(1:4,1:4))) - A(1:4,1:4))*B(1:4);
P2 = P1;

P1_tie = -A(end,1)/(s*(s*network(1).inertia + network(1).damping));
P2_tie = P1_tie;

G1 = -(P1*(1-P2_tie))/(1-P1_tie-P2_tie);

disp('Poles Calulated by Matlab of G1(s)')
roots(G1.Denominator{1})

disp('Poles Calulated by Matlab of P1(s)')
roots(P1.Denominator{1})
disp('Supposedly added Poles from tie-lines Calulated by Matlab')
roots([network(1).inertia network(1).damping 2e3])


disp('Poles Calulated by Matlab from state-space')
roots(global_.Denominator{1,1})

disp('Poles of ss (should be equal to G1(s))')
roots(global_.Denominator{1})


G2 = -(P1_tie*(P2))/(1-P1_tie-P2_tie);
roots(G2.Denominator{1})
roots(global_.Denominator{1,2})



-A(end,1)*A(1,end)