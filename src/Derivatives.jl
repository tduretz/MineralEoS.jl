@inline function compute_density_derivative(EoS, P, T, params, opts)
    return DifferentiationInterface.derivative(P -> density_volume(EoS, P, T, params; options = opts), AutoForwardDiff(), P)[1]
end

@inline function compute_density_derivative(EoS, P, T, params)
    return DifferentiationInterface.derivative(P -> density_volume(EoS, P, T, params), AutoForwardDiff(), P)[1]
end

@inline function compute_eff_bulk_modulus(EoS, P, T, params, opts)
    f, fprime = value_and_derivative(P -> density_volume(EoS, P, T, params; options = opts), AutoForwardDiff(), P)
    return f[1] / fprime[1]
end

@inline function compute_eff_bulk_modulus(EoS, P, T, params)
    f, fprime = value_and_derivative(P -> density_volume(EoS, P, T, params), AutoForwardDiff(), P)
    return f[1] / fprime[1]
end
