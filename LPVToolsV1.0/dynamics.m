function [xdot, mdot, rdot, vdot, qdot, wdot] = dynamics(x,T_eci,aeng,Anoz,P,F_eci,M_body,g_eci,J,Jdot,Mjet)

% INPUTS
% x      - state vector
% T_eci  - thrust force in ECI frame
% aeng   - mass depletion constant (Isp*g0)^-1
% Anoz   - engine nozzle area
% P      - atmospheric pressure
% F_eci  - sum of non-propulsive forces in ECI frame
% M_body - sum of moments in body frame
% g_eci  - gravity acceleration in ECI frame
% J      - inertia matrix
% Jdot   - inertia matrix derivative
% Mjet   - engine jet damping moment

% State vector
m=x(1);
r=x(2:4);
v=x(5:7);
q=x(8:11);
w=x(12:14);

% Propulsion
T=norm(T_eci);
mdot=-aeng*(T+P*Anoz);

% Translational ECI
F=T_eci+F_eci+m*g_eci;
vdot=F/m;
rdot=v;

% Quaternions
qdot=1/2*[q(4) -q(3)  q(2);
          q(3)  q(4) -q(1);
         -q(2)  q(1)  q(4);
         -q(1) -q(2) -q(3)]*w;

% Rotational BODY
M=M_body-Mjet-Jdot*w;
wdot=J\(-sksim(w)*J*w+M);

% State derivative -> to be integrated
xdot=[mdot;rdot;vdot;qdot;wdot];

end


function M=sksim(w)

M=[   0 -w(3) w(2);
    w(3)   0 -w(1);
   -w(2) w(1)   0];

end