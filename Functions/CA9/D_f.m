function D = D_f(V,rho,states,thrusters)

Q = 0.5*rho*V^2;

D = [-Q/V zeros(1,thrusters);zeros(states,thrusters+1)];