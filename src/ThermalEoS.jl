using MuladdMacro

const Ζ4 = pi^4 / 90
@inline expm1(x) = exp(x) - 1

@inline function U_Einstein_no_opt(T, θE, R)
    return R * θE / (exp(θE / T) - 1)
end

@inline function U_Einstein(T, θE, R)
    x = θE / T
    return R * θE / expm1(x)
end

@inline function thermal_pressure_no_opt(model::Einstein, V, T, materials)
    γ0 = materials.γ0
    θE0 = materials.θE
    T0 = materials.T0
    V0 = materials.V0
    q = materials.q
    Natom = materials.Natom
    R = materials.R
    γ = γ0 * (V / V0)^q
    if materials.qcompromised
        θE = θE0
    else
        θE = θE0 * (V / V0)^(-γ)
    end
    sca = 1.0e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa)
    P = sca * 3 * Natom * γ / V * (U_Einstein_no_opt(T, θE, R) - U_Einstein_no_opt(T0, θE, R))
    return P
end

@inline @muladd function thermal_pressure(model::Einstein, V, T, materials)
    γ0 = materials.γ0
    θE0 = materials.θE
    T0 = materials.T0
    V0 = materials.V0
    q = materials.q
    Natom = materials.Natom
    R = materials.R
    V̄ = V / V0
    γ = γ0 * fastpow(V̄, q)
    if materials.qcompromised
        θE = θE0
    else
        θE = θE0 * fastpow(V̄, -γ)
    end
    sca = 1.0e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa)
    P = sca * 3 * Natom * γ / V * (U_Einstein(T, θE, R) - U_Einstein(T0, θE, R))
    return P
end

########################################################################################

@inline function I3_series_no_opt(y; tol = 1.0e-12, kmax = 10_000)
    S = 6 * Ζ4
    for k in 1:kmax
        ek = exp(-k * y)
        term = ek * (y^3 / k + 3y^2 / k^2 + 6y / k^3 + 6 / k^4)
        S -= term
        if abs(term) < tol
            return S
        end
    end
    @warn "Series did not converge before kmax"
    return S
end

@inline @muladd function I3_series(y; tol = 1.0e-12, kmax = 10_000)
    S = 6 * Ζ4
    for k in 1:kmax
        ek = exp(-k * y)

        invk = 1 / k
        invk2 = invk * invk
        invk3 = invk2 * invk
        invk4 = invk2 * invk2
        y2 = y * y
        y3 = y2 * y
        term = ek * (y3 * invk + 3 * y2 * invk2 + 6 * y * invk3 + 6 * invk4)

        S -= term
        if abs(term) < tol
            return S
        end
    end
    @warn "Series did not converge before kmax"
    return S
end

@inline function U_Debye_no_opt(T, θD, R)
    return R * T * Debye_no_opt(θD / T)
end

@inline function U_Debye(T, θD, R)
    return R * T * Debye(θD / T)
end

@inline function Debye_no_opt(x)
    return 3 / (x^3) * I3_series_no_opt(x)
end

@inline function Debye(x)
    return 3 / (x * x * x) * I3_series(x)
end

@inline function thermal_pressure_no_opt(model::Debye, V, T, materials)
    γ0 = materials.γ0
    θD0 = materials.θD
    T0 = materials.T0
    V0 = materials.V0
    q = materials.q
    Natom = materials.Natom
    R = materials.R
    γ = γ0 * (V / V0)^q
    θD = θD0 * exp(- (γ0 / q) * ((V / V0)^q - 1.0))
    sca = 1.0e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa)
    P = sca * 3 * Natom * γ / V * (U_Debye_no_opt(T, θD, R) - U_Debye_no_opt(T0, θD, R))
    return P
end

@inline @muladd function thermal_pressure(model::Debye, V, T, materials)
    γ0 = materials.γ0
    θD0 = materials.θD
    T0 = materials.T0
    V0 = materials.V0
    q = materials.q
    Natom = materials.Natom
    R = materials.R
    V̄ = fastpow(V / V0, q)
    γ = γ0 * V̄
    θD = θD0 * exp(- (γ0 / q) * (V̄ - 1.0))
    sca = 1.0e-3 # 1e6/1e9 : (m3 -> cm3) / (GPa -> Pa)
    P = sca * 3 * Natom * γ / V * (U_Debye(T, θD, R) - U_Debye(T0, θD, R))
    return P
end
