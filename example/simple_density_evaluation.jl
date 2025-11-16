let
    # Define characteristic units
    # scales = (L=1e0, t=1e-2, σ=1e9, T=1)
    scales = (L=1e0, t=1e-0, σ=1e1, T=1)

    ρc     = scales.σ * scales.L * scales.t^2.0 / scales.L^3

    # Get material parameters from database
    params = assign_EoS_parameters(:OlivineFo90; sc=scales)

    # Define P, T
    P = 0e9  / scales.σ
    T = 850.0 / scales.T
    # T = 298.0 / scales.T

    # Call routine

    @info "Default: Einstein model"
    ρ, V = density_volume(Val(:complex), P, T, params)

    # Scale back
    @show ρ*ρc, V*scales.L^3

    @info "Option: Debye model"
    opts = (thermal_model=:Debye, mechanical_model=:BM2)
    ρ, V = density_volume(Val(:complex), P, T, params; options=opts)

    # Scale back
    @show  ρ*ρc, V*scales.L^3
    @show 3143.55
end