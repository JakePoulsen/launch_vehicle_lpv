clear all; close all; clc

% Taken from Table 10.1 on page 170 in Simplicio.

Q_d = 52.48*10^3;
S_ref = 7.14;
C_N_a = 2.007;
N_a = Q_d*S_ref*C_N_a;
J_N = 46*10^5;
CG = 8.919;
CP = 18.36;
l_CP = abs(CG-CP);
mu_a = (N_a/J_N)*l_CP;
g = 9.81;
m = 93.5*10^3;
theta_0 = 0.8125; % Derived from table
F_g = m*g;
V = 571.4; % For ascent alpha_0 = 0, so this should be accurate. V*cos(alpha_0)

% Values in B-matrix

F_TVC = 2386*10^3; % beta_0 is 0 during ascent. F_TVC*cos(beta_0)
PVP = 0.96; % point in origin of yz-plane
l_CG = abs(CG-PVP);
l_off = 1; % Considered to be reasonably


l_off_z_1 = 1;
l_off_z_2 = 1;

l_PVP_z_1 = sqrt(l_off_z_1^2+l_CG^2);
l_PVP_z_2 = sqrt(l_off_z_2^2+l_CG^2);

mu_PVP_1 = (F_TVC*l_PVP_z_1)/J_N;
mu_PVP_2 = (F_TVC*l_PVP_z_2)/J_N;

beta_off_z_1 = atan(l_off_z_1/l_CG);
beta_off_z_2 = atan(l_off_z_2/l_CG);

% Matrices
% 
A = [   0 1 0 0;
        mu_a -(l_CP*mu_a)/(V) 0 (mu_a/V);
        0 0 0 1;
        -(N_a+F_g*sin(theta_0))/m (l_CP*N_a)/(m*V)+V 0 -(N_a)/(m*V)];
% 
% B = [0  0;
%     -mu_PVP_1*cos(beta_off_z_1)  -mu_a/V;
%     0  0;
%     -(F_TVC*cos(beta_off_z_1))/m  (N_a)/(m*V)];
% 
% Bc = [0;-mu_PVP*cos(beta_off_z_1);0;-(F_TVC*cos(beta_off_z_1))/m];
%
max_pitch=2
max_pitch_rate=2
 C = diag([(max_pitch*pi/180)^-1 (max_pitch_rate*pi/180)^-1 (5000)^-1 (50)^-1]);
% 
 D = 0;


%% LPV Quadratic Stability
clc
close all



A0 = [0 1 0 0;(Q_d*S_ref*C_N_a*l_CP)/J_N -(l_CP*mu_a)/(V) 0 (Q_d*S_ref*C_N_a*l_CP)/(J_N*V);0 0 0 1;0 V 0 0];
A1 = [0 0 0 0;0 0 0 0;0 0 0 0; -(N_a+m*g*sin(theta_0))/m (l_CP*N_a)/(m*V) 0 -(N_a)/(m*V)];

m = 137*10^3-88*10^3; % Total mass - fuel mass in P80

Ap1 = [0 0 0 0;0 0 0 0;0 0 0 0; -(N_a+m*g*sin(theta_0))/m (l_CP*N_a)/(m*V)+V 0 -(N_a)/(m*V)];
Bp =@(m) [0 0;-mu_PVP_1*cos(beta_off_z_1) -mu_PVP_2*cos(beta_off_z_2);0 0;-(F_TVC*cos(beta_off_z_1))/m -(F_TVC*cos(beta_off_z_2))/m];
Bp1 = Bp(m)

Ep1 = [0;mu_a/V;0;(N_a)/(m*V)*0.1];
%Ep1 = zeros(4)
Fp1 = zeros(size(A,1),size(Ep1,2));
m = 137*10^3; % Total mass

Ap2 = [0 0 0 0;0 0 0 0;0 0 0 0; -(N_a+m*g*sin(theta_0))/m (l_CP*N_a)/(m*V)+V 0 -(N_a)/(m*V)];
Bp2 = Bp(m)

Ep2 = [0;mu_a/V;0;(N_a)/(m*V)*0.1];
%Ep2 = zeros(4)
Fp2 = zeros(4,1);

Af1 = A0+Ap1;
Af2 = A0+Ap2;

n = size(A,1) % Number of states
m = size(Bp(0),2) % Number of inputs w/o disturbance

X = sdpvar(n,n,'symmetric');
Y1 = sdpvar(m,n,'full');
Y2 = sdpvar(m,n,'full');
gamma = sdpvar(1);

lmi11 =@(A,B,Y) (A*X+B*Y)+(A*X+B*Y)';
lmi12 = @(E) E;
lmi13 = (C*X)';
lmi22 = @(E) -gamma*eye(size(E,2));
lmi23 = @(F) F';
lmi33 = @(A) -gamma*eye(size(A,1));

lmi_qc = @(A,B,Y,E,F) [lmi11(A,B,Y) lmi12(E) lmi13; lmi12(E)' lmi22(E) lmi23(F); lmi13' lmi23(F)' lmi33(A)] <= 0
lmi = [lmi_qc(Af1,Bp1,Y1,Ep1,Fp1), lmi_qc(Af2,Bp2,Y2,Ep2,Fp2)]
lmi = [lmi, X >= eye(n)*0.05];
lmi = [lmi,Y1.^2 <= 100];
lmi = [lmi,Y2.^2 <= 100];


opt = sdpsettings('solver','mosek','verbose',1);
d = optimize(lmi,gamma,opt);
%P_feas = double(P)
%P_eig = eig(P_feas)
X = double(X);
Y1 = double(Y1);
Y2 = double(Y2);

K1 = Y1*inv(X)
K2 = Y2*inv(X)

Ac1 = (Af1+Bp1*K1);
Ac2 = (Af2+Bp2*K2);
C=eye(4)
sys_p1 = ss(Ac1,Ep1,C,D);
sys_p2 = ss(Ac2,Ep2,C,D);

eig_p1 = eig(Ac1)
eig_p2 = eig(Ac2)

% Plotting
tfinal = 110;
sgtitle('Quadratic stabilization - State feedback','fontweight','bold','fontsize',16);
subplot(1,2,1)
step(sys_p1,tfinal)
title('Step function: Empty fuel tank')
subplot(1,2,2)
step(sys_p2,tfinal)
title('Step function: Full fuel tank')