using Test
using MineralEoS

@testset "Diamond - BM3    " begin
    params = assign_EoS_parameters(:Diamond)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 3.4200170958592153
    @test ρ ≈ 3514.982429343638
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 3.4091933464867825
    @test ρ ≈ 3526.1420454161403
end

@testset "Diamond - exp    " begin
    params = assign_EoS_parameters(:Diamond)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 3.420016933179745
    @test ρ ≈ 3514.9825965403193
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 3.3889675068800655
    @test ρ ≈ 3547.186562159455
end

@testset "OlivineFo90 - BM3" begin
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 43.89235739941078
    @test ρ ≈ 3249.8254468764276
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 43.258393151067864
    @test ρ ≈ 3297.4525776272108
end

@testset "OlivineFo90 - exp" begin
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 43.89235389473845
    @test ρ ≈ 3249.8257063652063
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 43.10351630031558
    @test ρ ≈ 3309.3007773696563
end
