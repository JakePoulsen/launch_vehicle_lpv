function Bd_matrix=Bd_pitch(x_CG,x_PVP,S_ref,J_y,m,mach,rho,v_x)

clc

Bd_fun_struct = load('Mats/save_Bd_pitch.mat');
    
Bd_matrix = double(subs(Bd_fun_struct.Bd_pitch));
    
end
