function U_Einstein(T, θE, R)
    return R*θE/(exp(θE/T) - 1)
end

function thermal_pressure(V, T, materials, phase)
    γ0    = materials.γ0[phase]
    θE0   = materials.θE[phase]
    T0    = materials.T0[phase]
    V0    = materials.V0[phase]
    q     = materials.q[phase]
    Natom = materials.Natom[phase]
    R     = materials.R
    γ     = γ0 * ((V)/V0)^q
    θE    = θE0*((V)/V0)^(-γ)
    sca   = 1e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa) 
    P     = sca * 3*Natom*γ/V*(U_Einstein(T, θE, R) - U_Einstein(T0, θE, R))
    return P  
end

function thermal_pressure(V, T, materials)
    γ0    = materials.γ0
    θE0   = materials.θE
    T0    = materials.T0
    V0    = materials.V0
    q     = materials.q
    Natom = materials.Natom
    R     = materials.R
    γ     = γ0 * ((V)/V0)^q
    θE    = θE0*((V)/V0)^(-γ)
    sca   = 1e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa) 
    P     = sca * 3*Natom*γ/V*(U_Einstein(T, θE, R) - U_Einstein(T0, θE, R))
    return P  
end