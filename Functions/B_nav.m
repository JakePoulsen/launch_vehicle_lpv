function B_u=B_rigid_pitch(x_CG,x_PVP,x_n,T_TVC,J_y,J_ny,m,m_n,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP)

I_N = J_ny+m_n*(x_PVP - x_n)^2;
ap = -T_TVC/m;
k2 = (m_n*(x_PVP - x_n))/m;
k1 = -(T_TVC*(x_CG - x_PVP))/(J_y);
k3 = (1/J_y)*(m_n*(x_PVP - x_n)*(x_CG - x_PVP)-I_N);


B_u = [
0,0;
k1,k3;
0,0;
ap,k2;
0,0;
0,0;
-Psi_z1_PVP*T_TVC,I_N*Psi_theta1_PVP - Psi_z1_PVP*m_n*(x_PVP - x_n);
-Psi_z2_PVP*T_TVC,I_N*Psi_theta2_PVP - Psi_z2_PVP*m_n*(x_PVP - x_n);
];
end