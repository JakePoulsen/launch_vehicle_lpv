clear;close;clc;
a_max=10
m_max=4
a=0:a_max
m=0:0.5:m_max
[ag mg]=meshgrid(a,m)
drag= zeros(length(a),length(m))
for ai = 1:size(a,2)
    
    for mi = 1:size(m,2)
        
        drag(ai,mi)=Drag_Coeff(a(ai),m(mi))
        lift(ai,mi)=Lift_Coeff(a(ai),m(mi))
        cp(ai,mi)=Xcenterpressure(a(ai),m(mi))
    end
end
hold on

subplot(3,1,1)

ax=gca
plot(drag)
ylabel('C_D')
c=colormap(jet(length(m)))
ax.ColorOrder = c;
xlabel("\alpha")

subplot(3,1,2)
ax=gca
plot(lift)
ylabel('C_L')
c=colormap(jet(length(m)))
ax.ColorOrder = c;
xlabel("\alpha")

subplot(3,1,3)
ax=gca
plot(cp)
ylabel('CP')
c=colormap(jet(length(m)))
ax.ColorOrder = c;
xlabel("\alpha")

sgtitle('Lift, Drag and Center of Pressure')