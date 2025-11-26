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

include("DataEoS.jl")
export assign_EoS_parameters

include("MechanicalEoS.jl")
export mechanical_pressure

include("ThermalEoS.jl")
export thermal_pressure, I3_series

include("Solver.jl")
export density_volume

include("Derivatives.jl")
export compute_density_derivative

end # module MineralEoS
