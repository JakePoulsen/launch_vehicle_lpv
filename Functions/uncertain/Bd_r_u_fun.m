function Bd_r=Bd_rigid_pitch(x_CG,x_PVP,S_ref,J_y,m,rho,v_x,idx,cd_coeff,aoa,x_CP)
Bd_r = [
,0;
,-((S_ref*rho*v_x*x_CG*polyval(cd_coeff{idx}, aoa))/2 - (S_ref*rho*v_x*x_CP*polyval(cd_coeff{idx}, aoa))/2)/J_y;
,0;
,-(S_ref*rho*v_x*polyval(cd_coeff{idx}, aoa))/(2*m);
];
end