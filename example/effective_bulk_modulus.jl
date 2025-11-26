using MineralEoS, Enzyme

let 
    # Define characteristic units
    scales = (L=1e0, t=1e-2, σ=1e9, T=1)
    # scales = (L=1e0, t=1e-0, σ=1e0, T=1)

    ρc     = scales.σ * scales.L * scales.t^2.0 / scales.L^3

    # Get material parameters from database
    params = assign_EoS_parameters(:OlivineFo90; sc=scales)

    # Define P, T
    P = 0e9   / scales.σ
    T = 850.0 / scales.T

    # 1. Let's compute the effective bulk modulus value for the simple EoS 

    # Set EoS type (:complex or :simple)
    EoS  = Val(:simple)

    # Compute density
    ρ, V = density_volume(EoS, P, T, params)

    # Compare K with ∂ρ∂P
    dP = 1.0 

    @time s    = Enzyme.autodiff( Enzyme.Forward, (EoS, P, T, params) -> density_volume(EoS, P, T, params), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
    ∂ρ∂P = s[1][1]
    Keff = ρ / ∂ρ∂P

    @info "Compare K and ρ*∂P∂ρ: simple EoS" 
    @show params.K * scales.σ, Keff * scales.σ

    # 2. Let's compute the effective bulk modulus value for the complex EoS 

    # Set EoS type (:complex or :simple)
    EoS  = Val(:complex)

    # Fine tuning
    opts = (thermal_model=:Debye, mechanical_model=:BM2)

    # Compute density
    ρ, V = density_volume(EoS, P, T, params; options=opts)

    # Compare K with ∂ρ∂P
    dP = 1.0 

    @time  s    = Enzyme.autodiff( Enzyme.Forward, (EoS, P, T, params) -> density_volume(EoS, P, T, params; options=opts), Const(EoS), Duplicated(P, dP), Const(T), Const(params))
    ∂ρ∂P = s[1][1]
    Keff = ρ / ∂ρ∂P

    @info "Compare K and ρ*∂P∂ρ: complex EoS" 
    @show params.K * scales.σ, Keff * scales.σ
end