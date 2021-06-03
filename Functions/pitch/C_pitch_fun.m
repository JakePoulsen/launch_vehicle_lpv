function C_matrix=C_pitch(x_CG,x_INS,rho,v_x,Psi_z1_INS,Psi_z2_INS,Psi_theta1_INS,Psi_theta2_INS)

clc

C_fun_struct = load('Mats/save_C_pitch.mat');
    
C_matrix = double(subs(C_fun_struct.C_pitch));
    
end
