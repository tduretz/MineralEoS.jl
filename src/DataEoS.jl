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

function assign_EoS_parameters(mineral; sc=(σ=1.0, L=1.0, t=1.0, T=1.0))

    m   = sc.σ * sc.L * sc.t^2.0
    J   = m * sc.L^2.0 / sc.t^2.0

    if mineral===:Diamond 
        params = (
            R      = 8.31415  / (J/sc.T),
            #---------------#
            ρ0     = 3515.0   / (m/sc.L^3),
            V0     = 3.42     / sc.L^3,
            K      = 444.0e9  / sc.σ ,
            ∂K∂P   = 4.0,
            ∂2K∂P2 = -0.0088  / (1/sc.σ),
            #---------------#
            θE     = 1500.0   / sc.T,
            γ0     = 0.9726,
            T0     = 298.15   /  sc.T,
            q      = 1.0,
            Natom  = 1.0,
        )
    elseif mineral===:OlivineFo90
        params = (
            R      = 8.31415  / (J/sc.T),

            # Birch-Murnaghan
            ρ0     = 3250     / (m/sc.L^3),
            V0     = 43.8900  / sc.L^3,
            K      = 126.30e9 / sc.σ ,
            ∂K∂P   = 4.54,
            ∂2K∂P2 = -0.0374  / (1/sc.σ),

            # Einstein's model for Pthermal
            θE     = 471.0    /  sc.T,
            γ0     = 1.044,
            T0     = 298.15   /  sc.T,
            q      = 1.88,
            Natom  = 7.0,
        )
    else
        error("Mineral phase not yet implemented!")
    end
    return params
end