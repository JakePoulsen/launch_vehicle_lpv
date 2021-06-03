clear all;
%% Frame transformations
Omega_earth=(2*pi)/(86164.1)%Full earth rotation=86164.1 seconds
earth_r=6364.609*10^3
Omega_I=Omega_earth*[0 0 1]'
C_ie_t=@(t) [cos(Omega_earth*t)     sin(Omega_earth*t)  0;...
    -sin(Omega_earth*t)    cos(Omega_earth*t)  0;...
    0                      0                   1];
t0=0;
C_IE0=C_ie_t(t0);
L_lla=[5.305 geod2geoc(52.834,0) 0]
p_el=lla2ecef(L_lla)
earth_r=norm(p_el)
C_LE_ned= dcmecef2ned(L_lla(1), L_lla(2))'
L_az=deg2rad(-0.0165);
C_LE=C_LE_ned*angle2dcm( pi, L_az,0,'XZY' )
lp_r0=p_el
lp_q0=zeros(4)
b_q0=angle2dcm( 0, pi*0.5, 0) % initial body in launchpad frame

sim_dcm=angle2dcm( 0, 0 ,pi*0.5)
sim_quat=dcm2quat(sim_dcm)
sim_cam_offset=[0 0 50]
%% Parameters
% mass
m_payload=1500;
m_fairing=540;
m_payload_adapter=78;
m_upper_stage=688;%dry/NTO/UDMH
m_3rd_stage=12000;%total
m_2nd_stage=26300;%total
m_1st_stage=7408;%dry
m_1st_stage_fuel=88400;%fuel

m_rocket=m_fairing+m_payload_adapter+m_upper_stage+m_3rd_stage+...
    m_2nd_stage+m_1st_stage

m0=m_rocket+m_1st_stage_fuel

% position
r0=C_IE0'*lp_r0'
dcm0=C_IE0'*C_LE*b_q0 % initial direction cosine matrix
%dcm0=eye(3)
q0=dcm2quat(dcm0)'
v0=cross(Omega_I,r0)
w0=flip(dcm0'*Omega_I*-1)

x0=[m0;r0;v0;q0;w0]

Thrust=2200*10^3;       % [kN] first state Vega
g0 = 9.81;          % [m/s^2] Constant at its sea-level value
r_noz=1.5*(7.5/12);              % [m] nozzle radius
Isp=280;                % [s] specific impulse
Anoz=pi*r_noz^2;                   % Anoz   - engine nozzle area
aeng=(Isp*g0)^-1;                % aeng   - mass depletion constant (Isp*g0)^-1
P=1;                             % P      - atmospheric pressure
%% Measurements

% Main engine thrusters
thruster_radius=1;
P_pvp1=[0 thruster_radius 0];
P_pvp2=[0 0 thruster_radius];
P_pvp3=[0 -thruster_radius 0];
P_pvp4=[0 0 -thruster_radius];

% Interstage Distance
interstage_12=2;
interstage_23=1.1;
interstage_3u=0.6;

% Payload
r_payload=1.5/2;
l_payload=5;

% Fairing
r_fairing=2.6/2;
l_fairing=7.88;

% Payload
r_payload_adapter=2.18/2;
l_payload_adapter=1.461;

% Upper stage 
r_upper_stage = 2.18/2;
l_upper_stage = 2.04;

% Third stage
r_3rd_stage=1.9/2;
l_3rd_stage=4.12;

% Second stage
r_2nd_stage=1.9/2;
l_2nd_stage=8.39;

% First stage
r_1st_stage=3/2;
l_1st_stage=11.2;

l_total=l_fairing+l_upper_stage+l_3rd_stage+l_2nd_stage+l_1st_stage-interstage_12-interstage_23-interstage_3u

% Begining x-coordinate
x_thruster=0
x_1st_stage=0
x_2nd_stage=l_1st_stage-interstage_12
x_3rd_stage=x_2nd_stage+l_2nd_stage-interstage_23
x_upper_stage=x_3rd_stage+l_3rd_stage-interstage_3u
x_payload_adapter=x_upper_stage+l_upper_stage
x_fairing=x_upper_stage+l_upper_stage
x_payload=x_payload_adapter+l_payload_adapter

% Center of gravtity
cg_1st_stage=x_1st_stage+l_1st_stage/2
cg_2nd_stage=x_2nd_stage+l_2nd_stage/2
cg_3rd_stage=x_3rd_stage+l_3rd_stage/2
cg_upper_stage=x_upper_stage+l_upper_stage/2
cg_payload_adapter=x_payload_adapter+l_payload_adapter/2
cg_fairing=x_fairing+l_fairing/2
cg_payload=x_payload+l_payload/2
cg_1st_stage_fuel=x_1st_stage+l_1st_stage/2%make variable
cg_rocket=(cg_1st_stage*m_1st_stage...
    +cg_2nd_stage*m_2nd_stage+cg_3rd_stage*m_3rd_stage...
    +cg_upper_stage*m_upper_stage+cg_payload_adapter*m_payload_adapter...
    +cg_fairing*m_fairing+cg_payload*m_payload)...
    /(m_rocket)
cg_total=(cg_rocket*(m_rocket)+cg_1st_stage*m_1st_stage_fuel)/m0

% Distance from center of gravity
d_1st_stage=abs(cg_1st_stage-cg_total)
d_2nd_stage=abs(cg_2nd_stage-cg_total)
d_3rd_stage=abs(cg_3rd_stage-cg_total)
d_upper_stage=abs(cg_upper_stage-cg_total)
d_payload_adapter=abs(cg_payload_adapter-cg_total)
d_fairing=abs(cg_fairing-cg_total)
d_payload=abs(cg_payload-cg_total)
d_1st_stage_fuel=abs(cg_1st_stage_fuel-cg_total)


%% Inertia

% Intertia from center of mass
Icm_1st_stage=1/4*m_1st_stage*r_1st_stage^2+1/12*m_1st_stage*l_1st_stage^2
Icm_2nd_stage=1/4*m_2nd_stage*r_2nd_stage^2+1/12*m_2nd_stage*l_2nd_stage^2
Icm_3rd_stage=1/4*m_3rd_stage*r_3rd_stage^2+1/12*m_3rd_stage*l_3rd_stage^2
Icm_upper_stage=1/4*m_upper_stage*r_upper_stage^2+1/12*m_upper_stage*l_upper_stage^2
Icm_payload_adapter=1/4*m_payload_adapter*r_payload_adapter^2+1/12*m_payload_adapter*l_payload_adapter^2
Icm_faring=1/4*m_fairing*r_fairing^2+1/12*m_fairing*l_fairing^2
Icm_payload=1/4*m_payload*r_payload^2+1/12*m_payload*l_payload^2
Icm_1st_stage_fuel=1/4*m_1st_stage_fuel*r_1st_stage^2 ...
    +1/12*m_1st_stage_fuel*l_1st_stage^2%make dynamic

% Inertia added form parallel axis theorem
Ipa_1st_stage=m_1st_stage*d_1st_stage^2
Ipa_2nd_stage=m_2nd_stage*d_2nd_stage^2
Ipa_3rd_stage=m_3rd_stage*d_3rd_stage^2
Ipa_upper_stage=m_upper_stage*d_upper_stage^2
Ipa_payload_adapter=m_payload_adapter*d_payload_adapter^2
Ipa_fairing=m_fairing*d_fairing^2
Ipa_payload=m_payload*d_payload^2
Ipa_1st_stage_fuel=m_1st_stage_fuel*d_1st_stage_fuel^2

Jx_rocket=0.5*m_payload*r_payload^2+...
    m_fairing*r_fairing^2+...
    0.5*m_payload_adapter*r_payload_adapter^2+...
    0.5*m_upper_stage*r_upper_stage^2+...
    0.5*m_3rd_stage*r_3rd_stage^2+...
    0.5*m_2nd_stage*r_2nd_stage^2+...
    m_1st_stage*r_1st_stage^2;

Jyz_rocket=Icm_1st_stage+Icm_2nd_stage+Icm_3rd_stage+Icm_upper_stage...
    +Icm_payload_adapter+Icm_faring+Icm_payload...
    +Ipa_1st_stage+Ipa_2nd_stage+Ipa_3rd_stage+Ipa_upper_stage...
    +Ipa_payload_adapter+Ipa_fairing+Ipa_payload

J_rocket=diag([Jx_rocket,Jyz_rocket,Jyz_rocket])

%total
Jx_total=Jx_rocket+0.5*m_1st_stage_fuel*r_1st_stage^2


Jyz_total=Jyz_rocket+Icm_1st_stage_fuel+Ipa_1st_stage_fuel

J_total=diag([Jx_total,Jyz_total,Jyz_total])

%%
%load trajectory
Trajectory=load('guidanceLV.mat');
traj_t=Trajectory.out_time;
traj_pitch=Trajectory.out_pitch;
traj_thrust=Trajectory.out_thrust;
%% Controller
% pitch
K1_pitch_1=     [36.6899   10.1404    0    1.0732*0.1;...
    -144.4206  -49.2549   0   -0.7775*0.1];
K1_pitch_2=    [  220.5717    40.9126   -0.0020   -0.1549;
    -120.3241   -10.8772    0.0020    0.1554]*-0.02;

K1_pitch_3= [   -1.0605   -0.032    0.00002    0.0011;
                -1.0605   -0.032    0.00002    0.0011]*5

K1_pitch_m=  [-1  -0.2    0    0;...
    -1  -0.2    0  0]*5;
K1_pitch =K1_pitch_3*1

%yaw
K1_yaw =K1_pitch_m*1

%roll
K1_roll=[-1   -0.2]