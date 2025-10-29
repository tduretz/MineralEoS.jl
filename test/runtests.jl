using Test
using MineralEoS

@testset "Diamond - BM3    " begin
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


@testset "Diamond - exp    " begin
    params = assign_EoS_parameters(:Diamond)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 3.422747016882794
    @test ρ ≈ 3512.178943025764
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 3.391672804877551
    @test ρ ≈ 3544.357221814621
end

@testset "OlivineFo90 - BM3" begin
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 43.89235734686477
    @test ρ ≈ 3249.8254507669762
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params)
    @test V ≈ 43.26002004420757
    @test ρ ≈ 3297.328569294076
end

@testset "OlivineFo90 - exp" begin
    params = assign_EoS_parameters(:OlivineFo90)
    # Low T, low P
    P = 0.0e9
    T = 300.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 44.24449872017864
    @test ρ ≈ 3223.9601334876206
    # High T, high P
    P = 5.0e9
    T = 1100.0
    ρ, V = density_volume(P, T, params; EoS=:exp)
    @test V ≈ 43.44933234517466
    @test ρ ≈ 3282.9618385572594
end
