using MineralEoS, CairoMakie
import LinearAlgebra: norm

#  Diamond
let
    # Get material parameters from database
    params = assign_EoS_parameters(:Diamond)

    # Olivine Fo90-92: data produced by Prof. Ross Angel using EOSFIT7c
    V_EoSfit = [ 
        3.42002  3.42061  3.42144  3.42250  3.42377  3.42522  3.42683  3.42857  3.43043  3.43240  3.43446  3.43660  3.43880  3.44107  3.44340  3.44577  3.44819
        3.41236  3.41294  3.41376  3.41481  3.41607  3.41750  3.41909  3.42081  3.42266  3.42460  3.42664  3.42875  3.43093  3.43317  3.43547  3.43782  3.44021
        3.40478  3.40536  3.40617  3.40721  3.40845  3.40987  3.41144  3.41315  3.41497  3.41689  3.41890  3.42099  3.42315  3.42537  3.42764  3.42996  3.43233
        3.39729  3.39786  3.39866  3.39969  3.40092  3.40232  3.40388  3.40556  3.40736  3.40927  3.41126  3.41332  3.41545  3.41765  3.41989  3.42219  3.42453
        3.38988  3.39044  3.39124  3.39225  3.39347  3.39486  3.39639  3.39806  3.39984  3.40173  3.40369  3.40574  3.40784  3.41001  3.41223  3.41451  3.41682
        3.38255  3.38311  3.38389  3.38490  3.38610  3.38747  3.38899  3.39064  3.39241  3.39427  3.39621  3.39823  3.40032  3.40246  3.40466  3.40691  3.40919
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
        ρ[i,j], V[i,j] = density_volume(ComplexEoS(), P[i], T[j], params)
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

        ax = Axis(fig[1,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$MineralEoS.jl - Diamond")
        hm = heatmap!(ax, T, P./1e9, V')
        Colorbar(fig[1, 2], hm, label = L"$V$ (cm³)" )
        
        ax = Axis(fig[2,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$EOSFIT7c - Diamond")
        hm = heatmap!(ax, T, P./1e9, V_EoSfit')
        Colorbar(fig[2, 2], hm, label = L"$V$ (cm³)" )

        ax = Axis(fig[3,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$Absolute difference ($\log_{10}$)")
        hm = heatmap!(ax, T, P./1e9, log10.(abs.(V'.-V_EoSfit')))
        Colorbar(fig[3, 2], hm, label = L"$ΔV$ (cm³)" )
        
        display(fig)

        ######################

        ax = Axis(fig[1,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$MineralEoS.jl - Diamond")
        hm = heatmap!(ax, T, P./1e9, ρ')
        Colorbar(fig[1, 4], hm, label = L"$ρ$ (kg/m³)" )
        
        ax = Axis(fig[2,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$EOSFIT7c - Diamond")
        hm = heatmap!(ax, T, P./1e9, ρ_EoSfit')
        Colorbar(fig[2, 4], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[3,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$Difference")
        hm = heatmap!(ax, T, P./1e9, (ρ'.-ρ_EoSfit'))
        Colorbar(fig[3, 4], hm, label = L"$Δρ$ (kg/m³)" )

        display(fig)
    end
    with_theme(Visualisation, theme_latexfonts())
end
