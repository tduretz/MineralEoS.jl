using MineralEoS, Chairmarks

let 
    params = assign_EoS_parameters(:OlivineFo90)
    V      = 1e-6

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
    # ...

end