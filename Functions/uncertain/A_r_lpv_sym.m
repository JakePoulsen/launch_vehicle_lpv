function A_sym=A_rigid_pitch(x_CG,x_CP,J_y,S_ref,m,mach,rho,v_x,g,cd_coeff,cl_coeff,aoa)
A_sym = [
0,1,0,0;
-((S_ref*cd_coeff*rho*v_x^2*x_CG)/2 - (S_ref*cd_coeff*rho*v_x^2*x_CP)/2)/J_y,((S_ref*rho*v_x*x_CG*(x_CG - x_CP)*cd_coeff)/2 - (S_ref*rho*v_x*x_CP*(x_CG - x_CP)*cd_coeff)/2)/J_y,0,-((S_ref*rho*v_x*x_CG*cd_coeff)/2 - (S_ref*rho*v_x*x_CP*cd_coeff)/2 + (S_ref*rho*v_x^2*x_CG*diff(cl_coeff, aoa))/2 - (S_ref*rho*v_x^2*x_CP*diff(cl_coeff, aoa))/2)/J_y;
0,0,0,1;
-(S_ref*cd_coeff*rho*v_x^2)/(2*m),-(v_x - (S_ref*rho*v_x*(x_CG - x_CP)*cd_coeff)/2)/m,0,-((S_ref*rho*v_x^2*diff(cl_coeff, aoa))/2 + (S_ref*rho*v_x*cd_coeff)/2)/m;
];
end