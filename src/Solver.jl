using Enzyme

@inline function residual(V, P, T, materials, options)
    # @show  thermal_pressure(Val(options.thermal_model), V, T, materials)
    return P - mechanical_pressure(Val(options.mechanical_model), V, materials) - thermal_pressure(Val(options.thermal_model), V, T, materials)
end

@inline function density_volume(::Val{:complex}, P, T, materials;   
    options = (thermal_model=:Einstein, mechanical_model=:BM3),
    niter = 20, tol = 1e-12)

    P    /= 1e9
    ρ0    = materials.ρ0
    V0    = materials.V0
    iter  = 0

    # Initial guess
    V     = V0
    r0    = 1.0
    iter  = 0
    err   = 1.0
    
    while (iter<niter && err>tol)
        iter += 1
        # Evaluate the function and the Jacobian: r, ∂r∂V
        J = Enzyme.jacobian(Enzyme.ForwardWithPrimal, residual, V, Const(P), Const(T), Const(materials), Const(options) )  
        r = J.val[1]
        if iter==1 r0 = r end
        err         = abs(r/r0)
        # @show r, J.derivs[1]
        # Newton update
        V -= r/J.derivs[1]
    end
    ρ = ρ0*V0/V
    return ρ, V
end

@inline function density_volume(::Val{:simple}, P, T, materials)
    ρ0 = materials.ρ0
    α  = materials.α
    K  = materials.K
    V0 = materials.V0
    T0 = materials.T0
    ρ  = ρ0* exp(P/K  - α*(T-T0))
    V  = V0*ρ0/ρ
    return ρ, V
end