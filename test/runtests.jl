using Test
using MineralEoS

@testset "Diamond - SimpleEoS()           " begin
    EoS = SimpleEoS()
    params = assign_EoS_parameters(:Diamond)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 3.420016933179745
    @test ρ ≈ 3514.9825965403193
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 3.3889675068800655
    @test ρ ≈ 3547.186562159455
end
@testset "Diamond - ComplexEoS()          " begin
    EoS = ComplexEoS()
    params = assign_EoS_parameters(:Diamond)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 3.4200170958592153
    @test ρ ≈ 3514.982429343638
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 3.4091933464867825
    @test ρ ≈ 3526.1420454161403
end

@testset "OlivineFo90 - SimpleEoS()       " begin
    EoS = SimpleEoS()
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 43.89235389473845
    @test ρ ≈ 3199.828387805742
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 43.10351630031558
    @test ρ ≈ 3258.3884577178155
end

@testset "OlivineFo90 - ComplexEoS()      " begin
    EoS = ComplexEoS()
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 43.89235739941078
    @test ρ ≈ 3199.828132309098
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(EoS, P, T, params)
    @test V ≈ 43.258393151067864
    @test ρ ≈ 3246.7225379714073
end
@testset "∂ρ∂P, dρdT  - SimpleEoS()       " begin
    EoS    = SimpleEoS()
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    dρdP, dρdT = compute_density_derivative(EoS, P, T, params)
    @test abs( (dρdP - params.ρ0/params.K) ) / dρdP < 1e-4
    @test abs( (dρdT - params.ρ0*params.α) ) / dρdT < 1e-4
end
@testset "∂ρ∂P, dρdT - ComplexEoS()       " begin
    EoS    = ComplexEoS()
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    dρdP, dρdT = compute_density_derivative(EoS, P, T, params)
    @test dρdP ≈ 2.5342682726163707e-8
    @test dρdT ≈ -0.08606190522067628
end

@testset "Eff. bulk modulus - SimpleEoS() " begin
    EoS = SimpleEoS()
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    Keff = compute_eff_bulk_modulus(EoS, P, T, params)
    @test Keff ≈ params.K
end
@testset "Eff. bulk modulus - ComplexEoS()" begin
    EoS = ComplexEoS()
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    Keff = compute_eff_bulk_modulus(EoS, P, T, params)
    @test Keff ≈ 1.2626240745245195e11
end


# EoS = SimpleEoS()
# params = assign_EoS_parameters(:OlivineFo90)
# # Low T, low P
# P = 0.0e9
# T = 300.0
# PT = SA[P, T]
# density_volume(EoS, P, T, params)

# @inline function compute_eff_bulk_modulus(EoS, PT::SVector{2}, params)
#     f, J = value_and_jacobian(PT -> density_volume(EoS, PT, params), AutoForwardDiff(), PT)
#     ρ = f[1]
#     dρdP = J[1, 1] 
#     K = ρ / dρdP
#     return K
# end

# ForwardDiff.jacobian(PT -> density_volume(EoS, PT, params), PT)

# DifferentiationInterface.derivative(P -> density_volume(EoS, P, T, params), AutoForwardDiff(), P)
# DifferentiationInterface.derivative(T -> density_volume(EoS, P, T, params), AutoForwardDiff(), T)
