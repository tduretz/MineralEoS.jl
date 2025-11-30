using Enzyme

@inline function compute_density_derivative(EoS, P, dP, T, params, opts)
    return Enzyme.autodiff(Enzyme.Forward, (EoS, P, T, params) -> density_volume(EoS, P, T, params; options = opts), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
end

@inline function compute_density_derivative(EoS, P, dP, T, params)
    return Enzyme.autodiff(Enzyme.Forward, (EoS, P, T, params) -> density_volume(EoS, P, T, params), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
end

@inline function compute_eff_bulk_modulus(EoS, P, dP, T, params, opts)
    s = Enzyme.autodiff(Enzyme.ForwardWithPrimal, (EoS, P, T, params) -> density_volume(EoS, P, T, params; options = opts), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
    ρ = s[2][1]
    ∂ρ∂P = s[1][1]
    return ρ / ∂ρ∂P
end

@inline function compute_eff_bulk_modulus(EoS, P, dP, T, params)
    s = Enzyme.autodiff(Enzyme.ForwardWithPrimal, (EoS, P, T, params) -> density_volume(EoS, P, T, params), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
    ρ = s[2][1]
    ∂ρ∂P = s[1][1]
    return ρ / ∂ρ∂P
end
