@inline function compute_density_derivative(EoS, P, T, params, opts)
    return value_and_derivative(P -> density_volume(EoS, P, T, params; options = opts), P)
end

@inline function compute_density_derivative(EoS, P, T, params)
    return value_and_derivative(P -> density_volume(EoS, P, T, params), P)
end

@inline function compute_eff_bulk_modulus(EoS, P, T, params, opts)
    ρ, ∂ρ∂P = value_and_derivative(P -> density_volume(EoS, P, T, params; options = opts), P)
    return ρ / ∂ρ∂P
end

@inline function compute_eff_bulk_modulus(EoS, P, T, params)
    ρ, ∂ρ∂P = value_and_derivative(P -> density_volume(EoS, P, T, params), P)
    return ρ / ∂ρ∂P
end
