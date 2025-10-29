let
    # Define characteristic units
    scales = (L=1e0, t=1e-2, σ=1e9, T=1)
    ρc     = scales.σ * scales.L * scales.t^2.0 / scales.L^3

    # Get material parameters from database
    params = assign_EoS_parameters(:OlivineFo90; sc=scales)

    # Define P, T
    P = 0.0e9  / scales.σ
    T = 1100.0 / scales.T

    # Call routine
    ρ, V = density_volume(P, T, params)

    # Scale back
    ρ*ρc
end