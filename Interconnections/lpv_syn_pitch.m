function [Krb Gamma] = lpv_pitch(LVic,v_x)
b1 = basis(1,0); % ganger den på P0. (1*P0)
b2 = basis(v_x,1); % min parameter jeg ganger på: v_x*P1, samlet set bliver det p0 + v_x*P1
b3 = basis(v_x^2,2*v_x);
Xb = [b1; b2; b3];
% Xb = [b1];
Yb = Xb;
% opt = lpvsynOptions('BackOffFactor',1.02);

[Krb,Gamma,Info] = lpvsyn(LVic,4,1,Xb,Yb);

