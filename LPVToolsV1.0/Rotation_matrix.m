
phi=pi
theta=pi/2
psi=0
syms theta psi phi
Rz=[cos(phi) sin(phi) 0; -sin(phi) cos(phi) 0; 0 0 1]
Rx=[1 0 0; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)]
Ry=[cos(psi) 0 sin(psi); 0 1 0; -sin(psi) 0 cos(psi)]
R=Ry*Rz
R*[1 0 0]'
%%
C=angle2dcm(theta, psi,0,'XZY' )*[1 0 0]'
%%
base=eye(3)
base_l=R*base
col=eye(3);

plotRefBase([0 0 0],base,col)
plotRefBase([0 0 0],base_l,col)
%plotRefBase(p_el/norm(p_el),base_l,col)
%% Plot reference base
function plot = plotRefBase(origin,base,col)
for i = 1:3
quiver3(origin(1),origin(2),origin(3),base(1,i),base(2,i),base(3,i),'color',col(i,:))
hold on
end
end