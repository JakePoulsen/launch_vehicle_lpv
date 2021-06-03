function [pitch,yaw,roll] = PYR(base)
x_vec=base(1:3,1)
z_vec=base(1:3,3)

a=[1 0 0]'
b=[x_vec(1) 0 x_vec(3)]'
pitch=atan2(norm(cross(a,b)), dot(a,b))*sign(-x_vec(3))

Ry=[cos(pitch) 0 sin(pitch); 0 1 0; -sin(pitch) 0 cos(pitch)]
a=Ry*[1 0 0]'
b=x_vec
yaw=atan2(norm(cross(a,b)), dot(a,b))*sign(x_vec(2))

Rz=[cos(yaw) sin(yaw) 0; -sin(yaw) cos(yaw) 0; 0 0 1]
a=Ry*[0 0 1]'
b=Rz*z_vec
b=z_vec
roll=atan2(norm(cross(a,b)), dot(a,b))*sign(-b(2))