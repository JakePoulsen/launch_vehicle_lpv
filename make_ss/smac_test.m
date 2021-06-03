syms frequency damping
A=[0 1;-frequency^2 -2*damping*frequency];
B=[0;frequency^2];
C=[1 0;0 1];
D=[0;0];

abcd=sym2gss([A B;C D])