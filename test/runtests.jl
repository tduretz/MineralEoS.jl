using Test
using MineralEoS

@testset "Diamond    " begin
    params = assign_EoS_parameters(:Diamond)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 3.420017095603012
    @test ρ ≈ 3514.982429606956
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 3.4091933464867825
    @test ρ ≈ 3526.1420454161403
end

@testset "OlivineFo90" begin
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 43.89218076613333
    @test ρ ≈ 3249.838524998995
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 43.25983659100763
    @test ρ ≈ 3297.3425523676374
end
