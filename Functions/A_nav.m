function A_pitch = A_pitch(rho,S_ref,cn_fun,cl_fun,cd_fun,v_x,x_CG,x_CP,J_y,m,T_TVC,x_PVP,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP,zeta,omega_n1,omega_n2)

Q = 0.5*rho*v_x^2;
C_N = @(aoa,mach) cl_fun*cos(aoa) + cd_fun*sin(aoa);
C_N_a = @(alpha,M, d_alpha) (C_N(alpha+d_alpha,M) - C_N(alpha,M))./d_alpha;

N = (0.5*rho*v_x^2)*S_ref*(cn_fun);
acc = (T_TVC-cd_fun)/m;
a1 = -N/(m*v_x);
a2 = -a1*(x_CG-x_CP);
a3 = -acc+a1*v_x;
a4 = (x_CG-x_CP)*N/(J_y*v_x);
a5 = -a4*(x_CG-x_CP);
a6 = a4*v_x;

A_pitch = [
0,1,0,0,0,0,0,0;
a6,a5,0,a4,(T_TVC*(Psi_z1_PVP + Psi_theta1_PVP*(x_CG - x_PVP)))/J_y,(T_TVC*(Psi_z2_PVP + Psi_theta2_PVP*(x_CG - x_PVP)))/J_y,0,0;
0,0,0,1,0,0,0,0;
a3,a2,0,a1,(Psi_theta1_PVP*T_TVC)/m,(Psi_theta2_PVP*T_TVC)/m,0,0;
0,0,0,0,0,0,1,0;
0,0,0,0,0,0,0,1;
0,0,0,0,-omega_n1^2,0,-2*omega_n1*zeta,0;
0,0,0,0,0,-omega_n2^2,0,-2*omega_n2*zeta;
];
end