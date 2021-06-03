function Bd_r_u=Bd_rigid_pitch(rho,S_ref,cn_fun,m,v_x,x_CG,x_CP,J_y)

Q = 0.5*rho*v_x^2;
N = Q*S_ref*cn_fun;

a1 = -N/(m*v_x);
a4 = (x_CG-x_CP)*N/(J_y*v_x);

Bd_r_u = [
,0;
,-a4;
,0;
,-a1;
,0;
,0;
,0;
,0;
];
end