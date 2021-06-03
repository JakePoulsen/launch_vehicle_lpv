function Faero =  aerodynamic_force(Q, CD, CL)
S_ref=7.14
Faero = -(Q) * S_ref * [CD; 0; CL]
    return;
end