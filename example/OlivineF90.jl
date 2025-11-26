using MineralEoS, CairoMakie
import LinearAlgebra: norm

let
    # Get material parameters from database
    EoS    = ComplexEoS()
    params = assign_EoS_parameters(:OlivineFo90)

    # Olivine Fo90-92: data produced by Prof. Ross Angel using EOSFIT7c
    V_EoSfit = [ 
        43.8924  43.9534  44.0179  44.0849  44.1542  44.2253  44.2982  44.3725  44.4483  44.5256  44.6042  44.6841  44.7654  44.8480  44.9320  45.0173  45.1040
        43.5521  43.6101  43.6712  43.7348  43.8005  43.8680  43.9369  44.0074  44.0791  44.1522  44.2265  44.3020  44.3788  44.4567  44.5359  44.6163  44.6979
        43.2260  43.2810  43.3391  43.3996  43.4620  43.5261  43.5916  43.6585  43.7266  43.7959  43.8663  43.9378  44.0105  44.0842  44.1591  44.2350  44.3121
        42.9127  42.9651  43.0205  43.0781  43.1375  43.1985  43.2609  43.3245  43.3892  43.4551  43.5220  43.5899  43.6588  43.7288  43.7997  43.8717  43.9446
        42.6115  42.6615  42.7143  42.7692  42.8259  42.8841  42.9436  43.0042  43.0659  43.1286  43.1923  43.2569  43.3224  43.3889  43.4563  43.5246  43.5939
        42.3213  42.3691  42.4195  42.4721  42.5262  42.5818  42.6386  42.6965  42.7554  42.8152  42.8759  42.9376  43.0000  43.0634  43.1275  43.1925  43.2584
    ]
    ρ_EoSfit = params.V0 ./ V_EoSfit * params.ρ0
    nP, nT = size(V_EoSfit)

    # Define P-T grid
    P = LinRange(0, 5e9, nP)    # Pa
    T = LinRange(300, 1100, nT) # K

    # Allocate density and volume arrays
    V = zeros(nP, nT)
    ρ = zeros(nP, nT)

    # Loop over P-T grid
    for j in eachindex(T), i in eachindex(P)
        ρ[i,j], V[i,j] = density_volume(EoS, P[i], T[j], params)
    end

    # Compute error w.r.t. EOSFIT7c
    @info "Misfit w.r.t EOSFIT7c: $(norm(V_EoSfit .- V))"

    # Check corner values
    @info "Corner values in the P-T space"
    println( "T =  300 K, P = 0 GPa, V = ",  V_EoSfit[1,1],     " ", round(V[1,1]    , digits=5))
    println( "T =  300 K, P = 0 GPa, V = ",  round(ρ_EoSfit[1,1], digits=5),     " ", round(ρ[1,1]    , digits=5))
    println( "T = 1100 K, P = 5 GPa, V = ",  V_EoSfit[end,end], " ", round(V[end,end], digits=5))
    println( "T = 1100 K, P = 5 GPa, V = ",  round(ρ_EoSfit[end,end], digits=5), " ", round(ρ[end,end], digits=5))

    # Show parameters
    display(params)

    # Visualisation
    function Visualisation()

        fig = Figure(size=(600, 600))

        ax = Axis(fig[1,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$MineralEoS.jl - OlivineFo90")
        hm = heatmap!(ax, T, P./1e9, V')
        Colorbar(fig[1, 2], hm, label = L"$V$ (cm³)" )
        
        ax = Axis(fig[2,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$EOSFIT7c - OlivineFo90")
        hm = heatmap!(ax, T, P./1e9, V_EoSfit')
        Colorbar(fig[2, 2], hm, label = L"$V$ (cm³)" )

        ax = Axis(fig[3,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$Absolute difference ($\log_{10}$)")
        hm = heatmap!(ax, T, P./1e9, log10.(abs.(V'.-V_EoSfit')))
        Colorbar(fig[3, 2], hm, label = L"$ΔV$ (cm³)" )
        
        display(fig)

        ######################

        ax = Axis(fig[1,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$MineralEoS.jl - OlivineFo90")
        hm = heatmap!(ax, T, P./1e9, ρ')
        Colorbar(fig[1, 4], hm, label = L"$ρ$ (kg/m³)" )
        
        ax = Axis(fig[2,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$EOSFIT7c - OlivineFo90")
        hm = heatmap!(ax, T, P./1e9, ρ_EoSfit')
        Colorbar(fig[2, 4], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[3,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$Difference")
        hm = heatmap!(ax, T, P./1e9, (ρ'.-ρ_EoSfit'))
        Colorbar(fig[3, 4], hm, label = L"$Δρ$ (kg/m³)" )

        display(fig)
    end
    with_theme(Visualisation, theme_latexfonts())
end
