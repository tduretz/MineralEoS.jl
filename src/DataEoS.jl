# For the data please refer to EOSFIT7c written by
# RJ ANGEL, J GONZALEZ-PLATAS, M ALVARO, ML MAZZUCCHELLI

# Base.@kwdef mutable struct DataEoS
#     # Gas constant
#     R = 8.314
#     #-----------------------------------------#
#     # Birch-Murnaghan
#     ρ0      ::Union{Float64, Missing} = missing
#     V0      ::Union{Float64, Missing} = missing
#     K       ::Union{Float64, Missing} = missing
#     ∂K∂P    ::Union{Float64, Missing} = missing
#     ∂2K∂P2  ::Union{Float64, Missing} = missing
#     #-----------------------------------------#
#     # Einstein's model for Pthermal
#     θE      ::Union{Float64, Missing} = missing
#     γ0      ::Union{Float64, Missing} = missing
#     T0      ::Union{Float64, Missing} = missing
#     q       ::Union{Float64, Missing} = missing
#     Natom   ::Union{Float64, Missing} = missing
# end

function reference_expansivity(p)
    NA = 6.02214076e23 # Avogadro
    kB = 1.380649e-23
    u0 = p.θE / p.T0
    return 3 * p.γ0 * p.Natom * kB * u0^2 * (exp(u0) / (exp(u0) - 1)^2) / (p.K * p.V0 * 1.0e-6 / NA)
end

function assign_EoS_parameters(mineral; sc = (σ = 1.0, L = 1.0, t = 1.0, T = 1.0))

    m = sc.σ * sc.L * sc.t^2.0
    J = m * sc.L^2.0 / sc.t^2.0

    if mineral === :Diamond
        # EoS parameters
        params = (
            # Birch-Murnaghan
            ρ0 = 3515.0,
            V0 = 3.42,
            K = 444.0e9,
            ∂K∂P = 4.0,
            ∂2K∂P2 = -0.0088,
            #---------------#
            θE = 1500.0,
            θD = 2047.7707006369426,
            γ0 = 0.9726,
            T0 = 298.15,
            q = 1.0,
            Natom = 1.0,
            α = 2.672e-6,
            qcompromised = true,
        )
        # Compute reference expansivity
        α0 = reference_expansivity(params)

    elseif mineral === :OlivineFo90
        # EoS parameters
        params = (
            # Birch-Murnaghan
            ρ0 = 3200.0,
            V0 = 43.89,
            K = 126.302e9,
            ∂K∂P = 4.54207,
            ∂2K∂P2 = -0.0374,
            # Einstein's model for Pthermal
            θE = 471.0,
            θD = 643.0,
            γ0 = 1.044,
            T0 = 298.0,
            q = 1.88,
            Natom = 7.0,
            α = 1.0e-6, # dummy
            qcompromised = false,
        )
        # Compute reference expansivity
        α0 = reference_expansivity(params)

    else
        error("Mineral phase not yet implemented!")
    end

    # Apply dimensional scaling
    params_scaled = (
        # Gas constant
        R = 8.31415 / (J / sc.T),
        # Bolzmann constant
        kB = 1.380649e-23 / (J / sc.T),
        # Birch-Murnaghan
        ρ0 = params.ρ0 / (m / sc.L^3),
        V0 = params.V0 / sc.L^3,
        K = params.K / sc.σ,
        ∂K∂P = params.∂K∂P,
        ∂2K∂P2 = params.∂2K∂P2 / (1 / sc.σ),
        # Einstein's model for Pthermal
        θE = params.θE / sc.T,
        θD = params.θD / sc.T,
        γ0 = params.γ0,
        T0 = params.T0 / sc.T,
        q = params.q,
        Natom = params.Natom,
        α = α0 / (1 / sc.T),
        qcompromised = params.qcompromised,
    )
    return params_scaled
end
