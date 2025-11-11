using Enzyme

function residual(V, P, T, materials, phase)
    # @show P, mechanical_pressure(V, materials, phase), thermal_pressure(V, T, materials, phase)
    return P - mechanical_pressure(V, materials, phase) - thermal_pressure(V, T, materials, phase)
end

function density_volume_BM3(P, T, materials, phase; niter = 20, tol = 1e-12)

    P    /= 1e9
    ρ0    = materials.ρ0[phase]
    V0    = materials.V0[phase]

    niter = 20
    tol   = 1e-12
    iter  = 0

    # Initial guess
    V     = V0
    r0    = 1.0
    iter  = 0
    err   = 1.0
    
    while (iter<niter && err>tol)
        iter += 1
        # Evaluate the function and the Jacobian: r, ∂r∂V
        J = Enzyme.jacobian(Enzyme.ForwardWithPrimal, residual, V, Const(P), Const(T), Const(materials), Const(phase) )  
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


function residual(V, P, T, materials)
    return P - mechanical_pressure(V, materials) - thermal_pressure(V, T, materials)
end

function density_volume_BM3(P, T, materials; niter = 20, tol = 1e-12)

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
        J = Enzyme.jacobian(Enzyme.ForwardWithPrimal, residual, V, Const(P), Const(T), Const(materials) )  
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

function density_volume_exp(P, T, materials, phase)
    ρ0 = materials.ρ0[phase]
    α  = materials.α[phase]
    K  = materials.K[phase]
    T0 = materials.T0[phase]
    V0 = materials.V0[phase]
    ρ  = ρ0* exp(P/K  - α*(T-T0))
    V  = V0*ρ0/ρ
    return ρ, V
end  

function density_volume_exp(P, T, materials)
    ρ0 = materials.ρ0
    α  = materials.α
    K  = materials.K
    V0 = materials.V0
    T0 = materials.T0
    ρ  = ρ0* exp(P/K  - α*(T-T0))
    V  = V0*ρ0/ρ
    return ρ, V
end

function density_volume(P, T, materials, phase; EoS=:BM3)
    if EoS === :BM3
        return density_volume_BM3(P, T, materials, phase)
    elseif EoS === :exp
        return density_volume_exp(P, T, materials, phase)
    else
        error("The requested EoS is not implemented (EoS = :BM3 or EoS = :exp)")
    end
end

function density_volume(P, T, materials; EoS=:BM3)
    if EoS === :BM3
        return density_volume_BM3(P, T, materials)
    elseif EoS === :exp
        return density_volume_exp(P, T, materials)
    else
        error("The requested EoS is not implemented (EoS = :BM3 or EoS = :exp)")
    end
end