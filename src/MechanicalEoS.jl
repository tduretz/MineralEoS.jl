function mechanical_pressure(V, materials, phase)
    # Birch-Murnaghan EOS
    V0 = materials.V0[phase]
    K0 = materials.K[phase]/1e9
    Kp = materials.∂K∂P[phase]
    Kpp = materials.∂2K∂P2[phase]
    f = ((V0/(V))^(2/3) -1)/2
    P = 3*K0*f*(1+2*f)^(5/2) * (1 + 3/2*(Kp -4)*f + 3/2*(K0*Kpp + (Kp - 4)*(Kp-3) + 35/9)*f^2 )
    # P = 3*K0*f*(1+2*f)^(5/2) * (1 + 3/2*(Kp -4)*f )
    return P
end

function mechanical_pressure(V, materials)
    # Birch-Murnaghan EOS
    V0 = materials.V0
    K0 = materials.K/1e9
    Kp = materials.∂K∂P
    Kpp = materials.∂2K∂P2
    f = ((V0/(V))^(2/3) - 1.0)/2
    P = 3*K0*f*(1+2*f)^(5/2) * (1 + 3/2*(Kp -4)*f + 3/2*(K0*Kpp + (Kp - 4)*(Kp-3) + 35/9)*f^2 )
    # P = 3*K0*f*(1+2*f)^(5/2) * (1 + 3/2*(Kp -4)*f )
    return P
end