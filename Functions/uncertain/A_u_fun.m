function A_u=A_rigid_pitch(x_CG,x_CP,x_PVP,J_y,S_ref,m,mach,rho,v_x,g,T_TVC,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP,zeta,omega_n1,omega_n2,idx,cd_coeff,cl_coeff,aoa)
A_u = [
0,1,0,0,0,0,0,0;
-((S_ref*polyval(cd_coeff{idx}, aoa)*rho*v_x^2*x_CG)/2 - (S_ref*polyval(cd_coeff{idx}, aoa)*rho*v_x^2*x_CP)/2)/J_y,((S_ref*rho*v_x*x_CG*(x_CG - x_CP)*polyval(cd_coeff{idx}, aoa))/2 - (S_ref*rho*v_x*x_CP*(x_CG - x_CP)*polyval(cd_coeff{idx}, aoa))/2)/J_y,0,-((S_ref*rho*v_x*x_CG*polyval(cd_coeff{idx}, aoa))/2 - (S_ref*rho*v_x*x_CP*polyval(cd_coeff{idx}, aoa))/2 + (S_ref*rho*v_x^2*x_CG*polyval(polyder(cl_coeff{idx}), aoa))/2 - (S_ref*rho*v_x^2*x_CP*polyval(polyder(cl_coeff{idx}), aoa))/2)/J_y,T_TVC*(Psi_z1_PVP + Psi_theta1_PVP*(x_CG - x_PVP)),T_TVC*(Psi_z2_PVP + Psi_theta2_PVP*(x_CG - x_PVP)),0,0;
0,0,0,1,0,0,0,0;
-(S_ref*polyval(cd_coeff{idx}, aoa)*rho*v_x^2)/(2*m),-(v_x - (S_ref*rho*v_x*(x_CG - x_CP)*polyval(cd_coeff{idx}, aoa))/2)/m,0,-((S_ref*rho*v_x^2*polyval(polyder(cl_coeff{idx}), aoa))/2 + (S_ref*rho*v_x*polyval(cd_coeff{idx}, aoa))/2)/m,(Psi_theta1_PVP*T_TVC)/m,(Psi_theta2_PVP*T_TVC)/m,0,0;
0,0,0,0,0,0,1,0;
0,0,0,0,0,0,0,1;
0,0,0,0,-omega_n1^2,0,-2*omega_n1*zeta,0;
0,0,0,0,0,-omega_n2^2,0,-2*omega_n2*zeta;
];
end