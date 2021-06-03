clear all
clc
%Y = atmosphere(h, vel, CL)
%Y = [TM; rhoE; Mach; Kn; asound; d; PR; MU; RE; KT];
%Kinetic Temperature
%Density (kg/m^3)
% Mach Number
% Knudsen Number
% Speed of Sound (m/s)
% free-molecule flow % continuum flow % transition flow
% Static Pressure (N/m^2)
% Dynamic Viscosity Coeff. (N.s/m^2)
% Reynold’s Number
% KT = 2.64638e-3*TM^1.5/(TM + 245.4*10^(-12/TM));
% KT is the coefficient of thermal conductivity of the perfect gas
height = 200000;
x = 1:200000;
d = 0;
p = 0;

for i = 1:200000
    Y = atmosphere(i, 1, 1);
    p(i) = Y(7);
end

f2 = figure;
semilogy(x,p)

