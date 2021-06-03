function Dd_r_u=Dd_rigid_pitch(rho,v_x)

Q = 0.5*rho*v_x^2;

Dd_r_u = [
,-(Q/v_x);
,0;
,0;
,0;
,0;
];
end