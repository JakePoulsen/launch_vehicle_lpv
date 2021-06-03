function A=A_fun(V,M,pitch,m,J_N,rho,CG,CP,g,S_ref)
    C_N = @(alpha,M) Lift_Coeff( alpha, M)*cos(alpha) + Drag_Coeff(alpha,M)*sin(alpha);
    C_N_a = @(alpha,M, d_alpha) (C_N(alpha+d_alpha,M) - C_N(alpha,M))./d_alpha;
    Q = @(rho,v)0.5*rho*v^2;
    N_a= @(Q_d,S_ref,C_N_a) Q_d*S_ref*C_N_a;
    mu_a = @(N_a,J_N,l_CP) (N_a/J_N)*l_CP;
    F_g=@(g,m)g*m;

    l_CP = @(CG, CP) (CP-CG);

    d_alpha=deg2rad(0.001);
    alpha=0;
    C_N_a = C_N_a(alpha,M,d_alpha);
    N_a=N_a(Q(rho,V),S_ref,C_N_a);
    l_CP=l_CP(CG, CP);
    mu_a=mu_a(N_a,J_N,l_CP);
    F_g=F_g(g,m);
    
    A=     [0                           1                   0   0;
            mu_a                        -(l_CP*mu_a)/(V)    0   (mu_a/V);
            0                           0                   0   1;
            -(N_a+F_g*sin(pitch))/m   (l_CP*N_a)/(m*V)+V  0   -(N_a)/(m*V)];
end