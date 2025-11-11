using MineralEoS, CairoMakie
import LinearAlgebra: norm

# Olivine


# !!!!!!!!!!!!!!! check room conditions


let
    nP, nT = 6, 17

    # Define characteristic units
    scales = (σ = 1e9, L = 1e3, t = 1e12, T = 100.0)

    # Define P-T grid
    P = LinRange(0/scales.σ, 5e9/scales.σ, nP)    # Pa
    T = LinRange(300/scales.T, 1100/scales.T, nT) # K

    # Get material parameters from database
    params = assign_EoS_parameters(:OlivineFo90; sc=scales)

    # Allocate density and volume arrays
    V_BM3 = zeros(nP, nT)
    ρ_BM3 = zeros(nP, nT)
    V_exp = zeros(nP, nT)
    ρ_exp = zeros(nP, nT)

    # Loop over P-T grid and compute density using the BM3 model
    for j in eachindex(T), i in eachindex(P)
        ρ_BM3[i,j], V_BM3[i,j] = density_volume(P[i], T[j], params)
    end

    # Loop over P-T grid and compute density using the BM3 model
    for j in eachindex(T), i in eachindex(P)
        ρ_exp[i,j], V_exp[i,j] = density_volume(P[i], T[j], params; EoS=:exp)
    end

    # Visualisation
    function Visualisation()
        ρc = scales.σ * scales.L * scales.t^2.0 / scales.L^3

        fig = Figure(size=(400, 600))
        
        ax = Axis(fig[1,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$OlivineFo90 - BM3 Einstein model ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, ρ_BM3' * ρc)
        Colorbar(fig[1,2], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[2,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$OlivineFo90 - exponential model")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, ρ_exp' * ρc, colorrange=extrema(ρ_BM3 * ρc))
        Colorbar(fig[2,2], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[3,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$\rho_\text{BM3} - \rho_\text{exp}$")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM3 - ρ_exp)' * ρc)
        Colorbar(fig[3,2], hm, label = L"$Δρ$ (kg/m³)" )
        
        display(fig)
    end
    with_theme(Visualisation, theme_latexfonts())

    @info "Room condition volume (cm³)"
    x = density_volume(0.0/scales.σ, 298.0/scales.T, params; EoS=:exp)
    x[2]*scales.L^3

    x = density_volume(0.0/scales.σ, 298.0/scales.T, params; EoS=:BM3)
    x[2]*scales.L^3
end
