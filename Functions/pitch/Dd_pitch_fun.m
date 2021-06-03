function Dd_matrix=Dd_pitch(rho,v_x)

clc

Dd_fun_struct = load('Mats/save_Dd_pitch.mat');
    
Dd_matrix = double(subs(Dd_fun_struct.Dd_pitch));
    
end
