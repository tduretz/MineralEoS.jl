@inline compute_density_derivative(EoS, P, T, params) = compute_density_derivative(EoS, SA[P, T], params)
@inline compute_density_derivative(EoS, P, T, params, opts) = compute_density_derivative(EoS, SA[P, T], params, opts)

@inline function compute_density_derivative(EoS, PT::SVector{2}, params)
    J = DifferentiationInterface.jacobian(PT -> density_volume(EoS, PT, params), AutoForwardDiff(), PT)
    dρdP, dρdT = J[1, 1], J[1, 2]
    return dρdP, dρdT
end

@inline function compute_density_derivative(EoS, PT::SVector{2}, params, opts)
    J = DifferentiationInterface.jacobian(PT -> density_volume(EoS, PT, params; options = opts), AutoForwardDiff(), PT)
    dρdP, dρdT = J[1, 1], J[1, 2]
    return dρdP, dρdT
end

@inline compute_eff_bulk_modulus(EoS, P, T, params) = compute_eff_bulk_modulus(EoS, SA[P, T], params)
@inline compute_eff_bulk_modulus(EoS, P, T, params, opts) = compute_eff_bulk_modulus(EoS, SA[P, T], params, opts)

@inline function compute_eff_bulk_modulus(EoS, PT::SVector{2}, params)
    f, J = value_and_jacobian(PT -> density_volume(EoS, PT, params), AutoForwardDiff(), PT)
    ρ = f[1]
    dρdP = J[1, 1] 
    K = ρ / dρdP
    return K
end

@inline function compute_eff_bulk_modulus(EoS, PT::SVector{2}, params, opts)
    f, J = value_and_jacobian(PT -> density_volume(EoS, PT, params; options = opts), AutoForwardDiff(), PT)
    ρ = f[1]
    dρdP = J[1, 1] 
    K = ρ / dρdP
    return K
end



