const Ζ4 = pi^4/90

@inline function U_Einstein(T, θE, R)
    return R*θE/(exp(θE/T) - 1)
end

@inline function thermal_pressure(::Val{:Einstein}, V, T, materials)
    γ0    = materials.γ0
    θE0   = materials.θE
    T0    = materials.T0
    V0    = materials.V0
    q     = materials.q
    Natom = materials.Natom
    R     = materials.R
    γ     = γ0 * ((V)/V0)^q
    if materials.qcompromised
        θE    = θE0
    else
        θE    = θE0*((V)/V0)^(-γ)
    end
    sca   = 1e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa) 
    P     = sca * 3*Natom*γ/V * (U_Einstein(T, θE, R) - U_Einstein(T0, θE, R))
    return P  
end

@inline function I3_series(y; tol=1e-12, kmax=10_000)
    S = 6 * Ζ4
    for k in 1:kmax
        ek = exp(-k*y)
        term = ek * (y^3/k + 3y^2/k^2 + 6y/k^3 + 6/k^4)
        S -= term
        if abs(term) < tol
            return S
        end
    end
    @warn "Series did not converge before kmax"
    return S
end

@inline function U_Debye(T, θD, R)
    return R*T*Debye(θD/T)
end

@inline function Debye(x)
    return 3/x^3*I3_series(x) 
end

@inline function thermal_pressure(::Val{:Debye}, V, T, materials)
    γ0    = materials.γ0
    θD0   = materials.θD
    T0    = materials.T0
    V0    = materials.V0
    q     = materials.q
    Natom = materials.Natom
    R     = materials.R
    γ     = γ0 * (V/V0)^q
    θD    = θD0 * exp(- (γ0 / q) * ((V / V0)^q - 1.0))
    sca   = 1e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa) 
    P     = sca * 3*Natom*γ/V * (U_Debye(T,  θD, R) - U_Debye(T0, θD, R)) 
    return P  
end