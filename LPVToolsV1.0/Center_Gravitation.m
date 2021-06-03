function [CG J] = Center_Gravitation(m_fuel, J_rocket, cg_rocket, m_rocket, m0_fuel, l_1st_stage, x_1st_stage, r_1st_stage)
%Cg_r is the center of gravitation of the rocket.
%m_r is the mass of the rocket without fuel
%Cg_f is the center of gravitation of the fuel
%m_f is the mass of the fuel
%J is the inertia
%M is the total mass of the rocket and the fuel

m_procent_fuel=m_fuel/m0_fuel; % procent of max fuel

l_fuel= l_1st_stage*m_procent_fuel; %length of fuel

cg_1st_stage_fuel=(x_1st_stage)+(l_fuel)/(2)%make variable
CG_x=(m_rocket*cg_rocket+m_fuel*cg_1st_stage_fuel)/(m_rocket+m_fuel);
d_fuel=abs(cg_1st_stage_fuel-CG_x);
Jx_fuel=0.5*m_fuel*r_1st_stage^2;

Jyz_fuel=1/4*m_fuel*r_1st_stage^2+1/12*m_fuel...
    *l_fuel^2+m_fuel*d_fuel^2;



J_fuel = diag([Jx_fuel,Jyz_fuel,Jyz_fuel]);
J=J_fuel+J_rocket;
CG=[CG_x 0 0]';


    return;
end

