function B_u=B_rigid_pitch(x_CG,x_PVP,x_n,T_TVC,J_y,J_ny,m,m_n,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP)
B_u = [
0,0;
-(T_TVC*(x_CG - x_PVP))/J_y,-(J_ny + m_n*(x_PVP - x_n)^2 - m_n*(abs(x_PVP - x_n)^2)^(1/2)*(x_CG - x_PVP))/J_y;
0,0;
-T_TVC/m,-(m_n*(abs(x_PVP - x_n)^2)^(1/2))/m;
0,0;
0,0;
-Psi_z1_PVP*T_TVC,J_ny*Psi_theta1_PVP - Psi_z1_PVP*m_n*(abs(x_PVP - x_n)^2)^(1/2);
-Psi_z2_PVP*T_TVC,J_ny*Psi_theta2_PVP - Psi_z2_PVP*m_n*(abs(x_PVP - x_n)^2)^(1/2);
];
end