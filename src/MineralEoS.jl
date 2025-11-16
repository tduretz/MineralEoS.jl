module MineralEoS

include("DataEoS.jl")
export assign_EoS_parameters

include("MechanicalEoS.jl")
export mechanical_pressure

include("ThermalEoS.jl")
export thermal_pressure, I3_series

include("Solver.jl")
export density_volume

end # module MineralEoS
