@inline function mechanical_pressure(model::BM4, V, materials)
    # Birch-Murnaghan EOS
    V0  = materials.V0
    K0  = materials.K/1e9
    Kp  = materials.∂K∂P
    Kpp = materials.∂2K∂P2
    η   = (V0/(V))^(2/3)
    f   = 1/2 * (η - 1.0)
    P   = 3*K0*f*(1+2*f)^(5/2) * (1 + 3/2*(Kp -4)*f + 3/2*(K0*Kpp + (Kp - 4)*(Kp-3) + 35/9)*f^2 )
    return P
end

@inline function mechanical_pressure(model::BM3, V, materials)
    # Birch-Murnaghan EOS
    V0 = materials.V0
    K0 = materials.K/1e9
    Kp = materials.∂K∂P
    η  = (V0/(V))^(2/3)
    f  = 1/2 * (η - 1.0) 
    P  = 3*K0*f*(1+2*f)^(5/2) * (1 + 3/2*(Kp -4)*f )
    return P
end

@inline function mechanical_pressure(model::BM2, V, materials)
    # Birch-Murnaghan EOS
    V0 = materials.V0
    K0 = materials.K/1e9
    η  = (V0/(V))^(2/3)
    f  = 1/2 * (η - 1.0) 
    P  = 3*K0*f*(1+2*f)^(5/2) * (1 + 2f )^(5/2)
    return P
end