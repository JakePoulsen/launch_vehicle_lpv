function A_pitch = A_pitch(S_ref,rho,cn_diff,v_x,x_CG,x_CP,J_y,m,T_TVC,x_PVP,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP,zeta,omega_n1,omega_n2)
A_pitch = [
0,1,0,0,0,0,0,0;
-(S_ref*rho*cn_diff*v_x^2*(x_CG-x_CP))/(2*J_y),(S_ref*rho*cn_diff*v_x*(x_CG-x_CP)^2)/(2*J_y),0,-(S_ref*rho*cn_diff*v_x*(x_CG-x_CP))/(2*J_y),(T_TVC*(Psi_z1_PVP + Psi_theta1_PVP*(x_CG - x_PVP)))/J_y,(T_TVC*(Psi_z2_PVP + Psi_theta2_PVP*(x_CG - x_PVP)))/J_y,0,0;
0,0,0,1,0,0,0,0;
(S_ref*rho*cn_diff*v_x^2)/(2*m),-(S_ref*rho*cn_diff*v_x*(x_CG-x_CP))/(2*m),0,(S_ref*rho*cn_diff*v_x)/(2*m),(Psi_theta1_PVP*T_TVC)/m,(Psi_theta2_PVP*T_TVC)/m,0,0;
0,0,0,0,0,0,1,0;
0,0,0,0,0,0,0,1;
0,0,0,0,-omega_n1^2,0,-2*omega_n1*zeta,0;
0,0,0,0,0,-omega_n2^2,0,-2*omega_n2*zeta;
];
end