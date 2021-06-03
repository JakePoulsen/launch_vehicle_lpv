function CD = Drag_Coeff(Alpha_eff ,Mach)



CD = (Alpha_eff * 0.0035 + 0.21) * (1 + Mach/5.714285);
    return;
end
