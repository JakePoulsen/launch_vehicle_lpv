function B_r_u=B_rigid_pitch(x_CG,x_PVP,T_TVC,J_y,m)
B_r_u = [
,0;
,-(T_TVC*(x_CG - x_PVP))/J_y;
,0;
,-T_TVC/m;
];
end