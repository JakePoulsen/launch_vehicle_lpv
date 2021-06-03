function Bd_r_u=Bd_rigid_pitch(x_CG,x_PVP,S_ref,J_y,m,mach,rho,v_x,x_CP,cd_coeff,idx,aoa)
Bd_r_u = [
,0;
,-((S_ref*rho*v_x*x_CG*polyval(cd_coeff{idx}, aoa))/2 - (S_ref*rho*v_x*x_CP*polyval(cd_coeff{idx}, aoa))/2)/J_y;
,0;
,-(S_ref*rho*v_x*polyval(cd_coeff{idx}, aoa))/(2*m);
,0;
,0;
,0;
,0;
];
end