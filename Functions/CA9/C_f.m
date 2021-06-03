function C = C_f(V,rho,states)

Q = 0.5*rho*V^2;

 C = [Q 0 0 Q/V;eye(states)];