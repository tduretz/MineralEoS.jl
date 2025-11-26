using MineralEoS, Chairmarks

let 
    params = assign_EoS_parameters(:OlivineFo90)
    V      = 1e-6
    T      = 1000.

    # Mechanical pressure 

    @info "BM2"
    a = @b mechanical_pressure_no_opt($(BM2(), V, params)...) 
    @show a
    a = @b mechanical_pressure($(BM2(), V, params)...) 
    @show a

    check = mechanical_pressure((BM2(), V, params)...) ≈ mechanical_pressure_no_opt((BM2(), V, params)...) 
    @show check

    @info "BM3"
    a = @b mechanical_pressure_no_opt($(BM3(), V, params)...) 
    @show a
    a = @b mechanical_pressure($(BM3(), V, params)...) 
    @show a

    check = mechanical_pressure((BM3(), V, params)...) ≈ mechanical_pressure_no_opt((BM3(), V, params)...) 
    @show check

    @info "BM4"
    a = @b mechanical_pressure_no_opt($(BM4(), V, params)...) 
    @show a
    a = @b mechanical_pressure($(BM4(), V, params)...) 
    @show a

    check = mechanical_pressure((BM4(), V, params)...) ≈ mechanical_pressure_no_opt((BM4(), V, params)...) 
    @show check

    # Thermal pressure 
    
    @info "Einstein"
    a = @b thermal_pressure_no_opt($(Einstein(), V, T, params)...)
    @show a
    a = @b thermal_pressure($(Einstein(), V, T, params)...)
    @show a

    check = thermal_pressure_no_opt((Einstein(), V, T, params)...) ≈ thermal_pressure((Einstein(), V, T, params)...)
    @show check

    @info "Debye"
    a = @b thermal_pressure_no_opt($(Debye(), V, T, params)...)
    @show a
    a = @b thermal_pressure($(Debye(), V, T, params)...)
    @show a

    check = thermal_pressure_no_opt((Debye(), V, T, params)...) ≈ thermal_pressure((Debye(), V, T, params)...)
    @show check
    if check === false
        @show thermal_pressure_no_opt((Debye(), V, T, params)...)
        @show thermal_pressure((Debye(), V, T, params)...)
    end

end