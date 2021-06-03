function B = B_f(mass,J,thrust,CG,pvp)
l_CG = CG-pvp;


F_TVC=thrust;


% l_off_z_1 = offset_1;





mu_PVP_1 = (F_TVC*l_CG)/J;


% beta_off_z_1 = atan(l_off_z_1/l_CG);



B = [0;-mu_PVP_1*cos(0);0;-(F_TVC)/mass];
 
