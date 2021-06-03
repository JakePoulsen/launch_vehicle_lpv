function Bn_r=Bn_rigid_pitch(J_ny,m_n,x_PVP,x_n,x_CG,J_y,m)
Bn_r = [
,0;
,-(J_ny + m_n*(x_PVP - x_n)^2 - m_n*(abs(x_PVP - x_n)^2)^(1/2)*(x_CG - x_PVP))/J_y;
,0;
,-(m_n*(abs(x_PVP - x_n)^2)^(1/2))/m;
];
end