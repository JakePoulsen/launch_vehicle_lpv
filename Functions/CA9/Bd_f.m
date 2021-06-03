function Bd = Bd_f(V,M,mass,J_N,rho,CG,CP,S_ref)
C_N = @(alpha,M) Lift_Coeff( alpha, M)*cos(alpha) + Drag_Coeff(alpha,M)*sin(alpha);
C_N_a = @(alpha,M, d_alpha) (C_N(alpha+d_alpha,M) - C_N(alpha,M))./d_alpha;

l_CP = CP-CG;
d_alpha=deg2rad(0.001);
alpha=0;
C_N_a = C_N_a(alpha,M, d_alpha);
Q = 0.5*rho*V^2;
N_a=Q*S_ref*C_N_a;
mu_a =(N_a/J_N)*l_CP;

Bd=[0;-mu_a/V;0;(N_a)/(mass*V)];