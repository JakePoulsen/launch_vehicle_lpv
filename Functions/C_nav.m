function C_u=C_rigid_pitch(x_CG,x_INS,rho,v_x,Psi_z1_INS,Psi_z2_INS,Psi_theta1_INS,Psi_theta2_INS)

Q = (0.5*rho*v_x^2);

C_u = [
(0.5*rho*v_x^2),0,0,(0.5*rho*v_x),0,0,0,0;
1,0,0,0,Psi_theta1_INS,-Psi_theta2_INS,0,0;
0,1,0,0,0,0,Psi_theta1_INS,-Psi_theta2_INS;
x_CG - x_INS,0,1,0,Psi_z1_INS,Psi_z2_INS,0,0;
0,x_CG - x_INS,0,1,0,0,Psi_z1_INS,Psi_z2_INS;
];

end