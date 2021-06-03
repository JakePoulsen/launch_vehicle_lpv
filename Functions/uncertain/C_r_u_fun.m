function C_r_u=C_rigid_pitch(x_CG,x_INS,rho,v_x)
C_r_u = [
(rho*v_x^2)/2,0,0,(rho*v_x)/2;
1,0,0,0;
0,1,0,0;
x_CG - x_INS,0,1,0;
0,x_CG - x_INS,0,1;
];
end