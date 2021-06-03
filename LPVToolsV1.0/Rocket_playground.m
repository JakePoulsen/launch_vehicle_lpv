% VEGA Rocket
m_stage_gross = [95796, 25751,10948];% 1st, 2nd,3d
% First stage(Solid Fuel)
m_prop = 88365;       % [kg] Propellant mass
Isp    = 280 ;        % [s]  Specific impulse
d      = 3;           % [m]  Diameter
g0     = 9.81;        % [m/s^2] Constant at its sea-level value
m0  = 137000;         % [kg] Initial mass
A   = pi*d^2/4;       % [m^2]Frontal area
Cd  = 0.5 ;             % Drag coefficient,assumed to have the constant value
rh0 = 1.225;          % [kg/m^3]
H0 = 7500;            % [m] Density scale height
Re = 6378e3;          % [m] Earth's radius
hgr_turn = 200;       % [m] Rocket starts the gravity turn when h = hgr_turn
tburn = 106.8;        % [s] Fuell burn time, first stage
md = (m_prop)/tburn;  % [kg/s]Propellant mass flow rate
T   = md*(Isp*g0);    % [N] Thrust (mean)
mf = m0 - m_prop;     % [kg] Final mass of the rocket(first stage is empty)
t0 = 0;               % Rocket launch time
tf = t0 + tburn;      % The time when propellant is completely burned
%and the thrust goes to zero
t_range     = [t0,tf];  % Integration interval
%% Dynamic equations
%refrence frames
clf
clc
Omega_earth=(2*pi)/(86164.1)%Full earth rotation=86164.1 seconds
C_ie_t=@(t) [cos(Omega_earth*t)     sin(Omega_earth*t)  0;...
             -sin(Omega_earth*t)    cos(Omega_earth*t)  0;...
             0                      0                   0];
earth_r=6364.609*10^3
L_lla=[5.305 geod2geoc(52.834,0) 0]     
p_el=lla2ecef(L_lla)
earth_r=norm(p_el)
C_el_ned= dcmecef2ned(L_lla(1), L_lla(2))'
L_az=deg2rad(-0.0165)
C_el=C_el_ned*angle2dcm( pi, L_az,0,'XZY' )
lp_r0=[0 0 0]
lp_q0=zeros(4)


base=[0 0 1; 0 1 0; 1 0 0]
base=eye(3)
base_l=C_el*base
col=eye(3);


c1=angle2dcm( pi*0,0,0)
c2=angle2dcm( 0,pi*0.5,0)
c=c1*c2
b_q0=c

base_b=base_l*b_q0
plotRefBase([0 0 0],base,col)
plotRefBase(p_el/norm(p_el),base_l,col)
plotRefBase(p_el/norm(p_el),base_b,col)
%plotRefBase(p_el/norm(p_el),base_l1,col)
[X Y Z]=sphere;
surf(X,Y,Z,'FaceAlpha',0.3, 'EdgeColor','none');
lims=[-1 1]*2;
xlabel('x')
ylabel('y')
zlabel('z')
xlim(lims)
ylim(lims)
zlim(lims)
%%
Thrust=2200*10^3;       % [kN] first state Vega
g0     = 9.81;          % [m/s^2] Constant at its sea-level value
r_noz=1.5;              % [m] nozzle radius
Isp=280;                % [s] specific impulse
% INPUTS

aeng=(Isp*g0)^-1                % aeng   - mass depletion constant (Isp*g0)^-1
Anoz=pi*1.5^2                   % Anoz   - engine nozzle area
P=1                             % P      - atmospheric pressure
F_eci=0                           % F_eci  - sum of non-propulsive forces in ECI frame
M_body= [1,1,1]   % M_body - sum of moments in body frame
    % g_eci  - gravity acceleration in ECI frame
    % J      - inertia matrix
    % Jdot   - inertia matrix derivative
    % Mjet   - engine jet damping moment


%% Plot reference base
function plot = plotRefBase(origin,base,col)
for i = 1:3
quiver3(origin(1),origin(2),origin(3),base(1,i),base(2,i),base(3,i),'color',col(i,:))
hold on
end
end

