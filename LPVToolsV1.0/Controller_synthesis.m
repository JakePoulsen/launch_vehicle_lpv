%% Controller synthesis
%
clear;close;clc;
%% Params
% Functions
C_N = @(alpha,M) Lift_Coeff( alpha, M)*cos(alpha) + Drag_Coeff(alpha,M)*sin(alpha)
C_N_a = @(alpha,M, d_alpha) (C_N(alpha+d_alpha,M) - C_N(alpha,M))./d_alpha
Q= @(rho,v)0.5*rho*v^2
N_a= @(Q_d,S_ref,C_N_a) Q_d*S_ref*C_N_a
mu_a = @(N_a,J_N,l_CP) (N_a/J_N)*l_CP
F_g=@(g,m)g*m
l_PVP = @(CG, pvp) abs(pvp-CG)
l_CP = @(CG, CP) abs(CP-CG)


% Ranges
Q_r= [0 59197] % Dynamic pressure
V_r= [1 1900]   % Velocity
rho_r=[1.2 5.58e-4] % Air density
J_r=[4.658e+6 3.155+6] % Inertia
CG_r=[8.508 14.05]  % Center of gravity
CP_r=[30 28]   % Center of pressure
g_r = [9.78 9.61] % Gravity 
m_r = [1.354e+05 4.88e+04] % Mass
pitch_r =[pi/2 0] % Pitch
Thrust_r = [2.8118e+6 2.088e+6] % thrust
M_r = [0 5.7]

% Static
S_ref = 7.14;

% Nominal values
Q_n=max(Q_r)% Dynamic pressure
V_n=mean(V_r)   % Velocity
rho_n=mean(rho_r) % Air density
J_n=mean(J_r )% Inertia
CG_n=min(CG_r  )% Center of gravity
CP_n=min(CP_r  )% Center of pressure
g_n=mean(g_r  )% Gravity 
m_n=mean(m_r )% Mass
pitch_n=mean(pitch_r )% Pitch
Thrust_n=mean(Thrust_r  )% thrust
M_n=mean(M_r)

A_nom=A_f(V_n,M_n,pitch_n,m_n,J_n,Q_n,CG_n,CP_n,g_n,S_ref,C_N_a,N_a,l_CP,mu_a,F_g); 
variables=9;
%% Varying effect
index = zeros(2^variables,1);
eigens= zeros(2^variables,4);
singulars =zeros(2^variables,4);
A_array = cell(2^variables,1);
values = zeros(2^variables,variables+1);
valuenames=["V","M","pitch","m","J","Q","CG","CP","g"]
i=0;
for Q = Q_r
for V = V_r
%for rho = rho_r
for J = J_r
for CG = CG_r
for CP = CP_r
for g = g_r
for m = m_r
for pitch = pitch_r
%for thrust = thrust_r
for M = M_r
    
   i=i+1;
   index(i)=i;
   values(i,:)=[i,V,M,pitch,m,J,Q,CG,CP,g];
   A=A_f(V,M,pitch,m,J,Q,CG,CP,g,S_ref,C_N_a,N_a,l_CP,mu_a,F_g); 
   
   A_array{i}=A;
   eigens(i,:)=eig(A);   
   [~, s, ~]=svd(A);   
   singulars(i,:)=diag(s);                                    
                                    
                                    
                                    
end
end
end
end
end
end
end
end
end
%end
%end


subplot(3,1,1)
plot(index,abs(eigens(:,:)),'*')
legend

subplot(3,1,2)
plot(index,abs(singulars(:,:)),'*')
legend

subplot(3,1,3)
norm_values=values(:,2:end)./[max(V_r),max(M_r),max(pitch_r),max(m_r),max(J_r),max(Q_r),max(CG_r),max(CP_r),max(g_r)];
plot(index,norm_values,'*')
legend("V","M","pitch","m","J","Q","CG","CP","g")

% [V_high V_low]= splitHighLow(values,2,max(V_r))
% [M_high M_low]= splitHighLow(values,3,max(M_r))
% [pitch_high pitch_low]= splitHighLow(values,4,max(pitch_r))
% [m_high m_low]= splitHighLow(values,5,max(m_r))
% [J_high J_low]= splitHighLow(values,6,max(J_r))
% [Q_high Q_low]= splitHighLow(values,7,max(Q_r))
% [CG_high CG_low]= splitHighLow(values,8,max(CG_r))
% [CP_high CP_low]= splitHighLow(values,9,max(CP_r))
% [g_high g_low]= splitHighLow(values,10,max(g_r))

HL={}
for i = 1:variables

    [high low]= splitHighLow(values,1+i);
    hi=high(:,1);
    li=low(:,1);
    hs=singulars(hi,:)
    ls=singulars(li,:)
    singular_diff_norm=abs((hs-ls)./ls)
    subplot(variables,1,i)
    plot(singular_diff_norm,'.')
    ylabel(valuenames(i))
    HL(i,:)={high,low};
    %HL(i,2)=low;
    
end
sgtitle('Singular value scaling')


function [high low] =splitHighLow(M,col,high_val)
high_val=max(M(:,col))
high_index=M==high_val
high=M(high_index(:,col),:)
low=M(~high_index(:,col),:)
%high==low
end

function [Max, Max_i]=match_vector(search_vector,search_matrix)
matcher = repmat([0 search_vector],size(search_matrix,1),1);
match = search_matrix==matcher
[Max, Max_i]=max(sum(match,2))
end

function A=A_f(V,M,pitch,m,J_N,Q_d,CG,CP,g,S_ref,C_N_a,N_a,l_CP,mu_a,F_g)
    d_alpha=deg2rad(0.001);
    alpha=0;
    C_N_a = C_N_a(alpha,M,d_alpha) ;
    N_a=N_a(Q_d,S_ref,C_N_a);
    l_CP=l_CP(CG, CP);
    mu_a=mu_a(N_a,J_N,l_CP);
    F_g=F_g(g,m);
    
    A=     [0                           1                   0   0;
            mu_a                        -(l_CP*mu_a)/(V)    0   (mu_a/V);
            0                           0                   0   1;
            -(N_a+F_g*sin(pitch))/m   (l_CP*N_a)/(m*V)+V  0   -(N_a)/(m*V)];
end