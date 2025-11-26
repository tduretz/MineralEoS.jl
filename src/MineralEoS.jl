module MineralEoS

abstract type AbstractEoS end
struct SimpleEoS  <: AbstractEoS end
struct ComplexEoS <: AbstractEoS end
export SimpleEoS, ComplexEoS

abstract type AbstractMechanical end
struct BM4 <: AbstractMechanical end
struct BM3 <: AbstractMechanical end
struct BM2 <: AbstractMechanical end
export BM2, BM3, BM4

abstract type AbstractThermal end
struct Einstein <: AbstractThermal end
struct Debye    <: AbstractThermal end
export Einstein, Debye

# fast exponential from Albert
@inline fastpow(x::Number, n::Integer) = x^n

@inline function fastpow(x::Number, n::AbstractFloat)
    isinteger(n) && return x^Int(n)
    x > 0 && return exp(log(x) * n)
    return x^n
end
export fastpow

include("DataEoS.jl")
export assign_EoS_parameters

include("MechanicalEoS.jl")
export mechanical_pressure, mechanical_pressure_no_opt

include("ThermalEoS.jl")
export thermal_pressure, thermal_pressure_no_opt, I3_series, I3_series_no_opt

include("Solver.jl")
export density_volume

include("Derivatives.jl")
export compute_density_derivative, compute_eff_bulk_modulus

end # module MineralEoS
