function A_matrix=A_pitch(x_CG,x_CP,x_PVP,J_y,S_ref,m,mach,rho,v_x,g,T_TVC,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP,zeta,omega_n1,omega_n2)

clc

A_r_fun_struct = load('Mats/save_A_pitch.mat');
    
A_matrix = double(subs(A_r_fun_struct.A_pitch));
   
end
