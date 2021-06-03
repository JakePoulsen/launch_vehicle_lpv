function A_r_u=A_rigid_pitch(x_CG,x_CP,J_y,S_ref,m,mach,rho,v_x,g,idx,cd_coeff,cl_coeff,aoa)
A_r_u = [
0,1,0,0;
-((S_ref*polyval(cd_coeff{idx}, aoa)*rho*v_x^2*x_CG)/2 - (S_ref*polyval(cd_coeff{idx}, aoa)*rho*v_x^2*x_CP)/2)/J_y,((S_ref*rho*v_x*x_CG*(x_CG - x_CP)*polyval(cd_coeff{idx}, aoa))/2 - (S_ref*rho*v_x*x_CP*(x_CG - x_CP)*polyval(cd_coeff{idx}, aoa))/2)/J_y,0,-((S_ref*rho*v_x*x_CG*polyval(cd_coeff{idx}, aoa))/2 - (S_ref*rho*v_x*x_CP*polyval(cd_coeff{idx}, aoa))/2 + (S_ref*rho*v_x^2*x_CG*polyval(diff(cl_coeff{idx}), aoa))/2 - (S_ref*rho*v_x^2*x_CP*polyval(diff(cl_coeff{idx}), aoa))/2)/J_y;
0,0,0,1;
-(S_ref*polyval(cd_coeff{idx}, aoa)*rho*v_x^2)/(2*m),-(v_x - (S_ref*rho*v_x*(x_CG - x_CP)*polyval(cd_coeff{idx}, aoa))/2)/m,0,-((S_ref*rho*v_x^2*polyval(diff(cl_coeff{idx}), aoa))/2 + (S_ref*rho*v_x*polyval(cd_coeff{idx}, aoa))/2)/m;
];
end