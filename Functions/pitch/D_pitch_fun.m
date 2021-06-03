function D_matrix=D_pitch()

clc

D_fun_struct = load('Mats/save_D_pitch.mat');
    
D_matrix = double(subs(D_fun_struct.D_pitch));
    
end
