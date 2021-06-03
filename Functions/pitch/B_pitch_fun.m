function B_matrix=B_pitch(x_CG,x_PVP,x_n,T_TVC,J_y,J_ny,m,Psi_z1_PVP,Psi_z2_PVP,Psi_theta1_PVP,Psi_theta2_PVP)

clc

B_fun_struct = load('Mats/save_B_pitch.mat');
    
B_matrix = double(subs(B_fun_struct.B_pitch));
    
end
