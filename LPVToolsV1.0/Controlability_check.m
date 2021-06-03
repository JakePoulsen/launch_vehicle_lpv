%% Math for controlability of multi engine rocket
% Author: Frederik Rasmussen
clear;close;clc
% Variables
CG = 8.919     % Nominal CG for the rocket
r_off = 1  % Engine of set from center

A=[0,1,0,0;
    1.35566407480734,-0.0278007408165645,0,0.00142626415024444;
    0,0,0,1;
    -8.61415950171208,950.536067054665,0,-0.00185035166555339]

thrust=2449900 /2

% Values in B-matrix

F_TVC = 2449900 /2; % beta_0 is 0 during ascent. F_TVC*cos(beta_0)
PVP = 0; % point in origin of yz-plane
l_CG = abs(CG-PVP);
l_off = 1; % Considered to be reasonably
mass=92100;
J=2.329004577500000e+06

l_off_z_1 = l_off;
l_off_z_2 = l_off;

l_PVP_z_1 = sqrt(l_off_z_1^2+l_CG^2);
l_PVP_z_2 = sqrt(l_off_z_2^2+l_CG^2);

mu_PVP_1 = (F_TVC*l_PVP_z_1)/J;
mu_PVP_2 = (F_TVC*l_PVP_z_2)/J;

beta_off_z_1 = atan(l_off_z_1/l_CG);
beta_off_z_2 = atan(l_off_z_2/l_CG);

%rotaiton without translation

F_1= @(beta) sin(beta)*thrust
F_2= @(beta) sin(beta)*thrust

M_1= @(beta) F_1(beta)*l_PVP_z_1
M_2= @(beta) F_2(beta)*l_PVP_z_2

b=0:0.01:deg2rad(6.5)
m1=M_1(b)
m2=M_2(-b)
m_total=m1+m2
hold on
bp=rad2deg(b)
plot(bp,m_total)
plot(bp,m1)
plot(bp,m2)
title("Torque per radian")
figure
hold on
plot(bp,m_total./m1)
plot(bp,m_total./-m2)
title("Torque ratio")
figure
plot(bp,m2./m1)
title("Torque ratio")

ratio=l_PVP_z_1/l_PVP_z_2
%%




Bp =@(m)    [0                              0;
            -mu_PVP_1*cos(beta_off_z_1)     -mu_PVP_2*cos(beta_off_z_2);
            0                               0;
            -(F_TVC*cos(beta_off_z_1))/m    -(F_TVC*cos(beta_off_z_2))/m];
B = Bp(mass)

CTR=ctrb(A,B(:,1))
rank(CTR)
C=eye(length(A))

[Abar,Bbar,Cbar,T,k]=ctrbf(A,B(:,1),C)
k=place(A,B,[-4 -3 -2 -1]*1)

eig(A-B*k)