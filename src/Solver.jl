@inline function residual(V, P, T, materials, options)
    # @show  thermal_pressure(Val(options.thermal_model), V, T, materials)
    return P - mechanical_pressure(options.mechanical_model, V, materials) - thermal_pressure(options.thermal_model, V, T, materials)
end

@inline function density_volume(
        EoS::ComplexEoS, P, T, materials;
        options = (thermal_model = Einstein(), mechanical_model = BM3()),
        niter = 20, tol = 1.0e-12
    )

    density_volume(
        EoS, SA[P, T], materials;
        options = (thermal_model = options.thermal_model, mechanical_model = options.mechanical_model),
        niter = niter, tol = tol
    )

end

@inline function density_volume(
        ::ComplexEoS, PT, materials;
        options = (thermal_model = Einstein(), mechanical_model = BM3()),
        niter = 20, tol = 1.0e-12
    )

    P, T = PT[1], PT[2]
    P /= 1.0e9
    ρ0 = materials.ρ0
    V0 = materials.V0
    iter = 0

    # Initial guess
    V = V0
    V̄ = 1.0
    r0 = 1.0
    iter = 0
    err = 1.0

    while (iter < niter && err > tol)
        iter += 1
        # Evaluate the function and the Jacobian: r, ∂r∂V
        r, J = value_and_derivative(x -> residual(x, P, T, materials, options), AutoForwardDiff(), V)
        if iter == 1
            r0 = r
        end
        err = abs(r / r0)
        # Newton update
        V -= r / J
    end
    ρ = ρ0 * V0 / V
    return SA[ρ, V]
end

density_volume(EoS::SimpleEoS, P, T, materials) = density_volume(EoS, SA[P, T], materials)

@inline function density_volume(::SimpleEoS, PT::SVector{2}, materials)
    (; ρ0, α, K, V0, T0) = materials
    P, T = PT[1], PT[2]
    ρ = @fastmath ρ0 * exp(P / K - α * (T - T0))
    V = V0 * ρ0 / ρ
    return SA[ρ, V]
end
