using MineralEoS, CairoMakie
import LinearAlgebra: norm

let
    nP, nT = 6, 17

    # Define characteristic units
    scales = (σ = 1e9, L = 1e3, t = 1e12, T = 100.0)

    # Define P-T grid
    P = LinRange(0/scales.σ, 5e9/scales.σ, nP)    # Pa
    T = LinRange(300/scales.T, 1100/scales.T, nT) # K

    # Get material parameters from database
    EoS    = ComplexEoS()
    params = assign_EoS_parameters(:OlivineFo90; sc=scales)

    # Thermal models
    thermal_models = (Einstein(), Debye())

    # Mechanical models
    mechanical_models = (BM2(), BM3(), BM4())

    # Allocate density and volume arrays

    # Reference: simple model
    V_exp    = zeros(nP, nT)
    ρ_exp    = zeros(nP, nT)

    # Einstein 
    V_BM2E   = zeros(nP, nT)
    ρ_BM2E   = zeros(nP, nT)
    V_BM3E   = zeros(nP, nT)
    ρ_BM3E   = zeros(nP, nT)
    V_BM4E   = zeros(nP, nT)
    ρ_BM4E   = zeros(nP, nT)

    # Mie-Grüneisen-Debye 
    V_BM2MGD = zeros(nP, nT)
    ρ_BM2MGD = zeros(nP, nT)
    V_BM3MGD = zeros(nP, nT)
    ρ_BM3MGD = zeros(nP, nT)
    V_BM4MGD = zeros(nP, nT)
    ρ_BM4MGD = zeros(nP, nT)

    # Loop over P-T grid and compute density using the BM3 model
    for j in eachindex(T), i in eachindex(P)

        # Einstein 
        ρ_exp[i,j], V_exp[i,j] = density_volume(SimpleEoS(), P[i], T[j], params)

        # Einstein 
        opts = (thermal_model=Einstein(), mechanical_model=BM2())
        ρ_BM2E[i,j], V_BM2E[i,j] = density_volume(EoS, P[i], T[j], params; options=opts)

        opts = (thermal_model=Einstein(), mechanical_model=BM3())
        ρ_BM3E[i,j], V_BM3E[i,j] = density_volume(EoS, P[i], T[j], params; options=opts)

        opts = (thermal_model=Einstein(), mechanical_model=BM4())
        ρ_BM4E[i,j], V_BM4E[i,j] = density_volume(EoS, P[i], T[j], params; options=opts)

        # Mie-Grüneisen-Debye 
        opts = (thermal_model=Debye(), mechanical_model=BM2())
        ρ_BM2MGD[i,j], V_BM2MGD[i,j] = density_volume(EoS, P[i], T[j], params; options=opts)

        opts = (thermal_model=Debye(), mechanical_model=BM3())
        ρ_BM3MGD[i,j], V_BM3MGD[i,j] = density_volume(EoS, P[i], T[j], params; options=opts)

        opts = (thermal_model=Debye(), mechanical_model=BM4())
        ρ_BM4MGD[i,j], V_BM4MGD[i,j] = density_volume(EoS, P[i], T[j], params; options=opts)
    end

    # Visualisation
    function Visualisation()
        ρc = scales.σ * scales.L * scales.t^2.0 / scales.L^3

        fig = Figure(size=(600, 600))
        
        ax = Axis(fig[1,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$\Delta\rho - BM2-E ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM2E - ρ_exp)' * ρc)
        Colorbar(fig[1,2], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[2,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$\Delta\rho - BM3-E ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM3E - ρ_exp)' * ρc)
        Colorbar(fig[2,2], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[3,1], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$\Delta\rho - BM4-E ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM4E - ρ_exp)' * ρc)
        Colorbar(fig[3,2], hm, label = L"$ρ$ (kg/m³)" )

        #-------------------------------

        ax = Axis(fig[1,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$\Delta\rho - BM2-MGD ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM2MGD - ρ_exp)' * ρc)
        Colorbar(fig[1,4], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[2,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$\Delta\rho - BM3-MGD ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM3MGD - ρ_exp)' * ρc)
        Colorbar(fig[2,4], hm, label = L"$ρ$ (kg/m³)" )

        ax = Axis(fig[3,3], xlabel=L"$T$ (K)", ylabel=L"$P$ (GPa)", title=L"$$\Delta\rho - BM4-MGD ")
        hm = heatmap!(ax, T * scales.T, P * scales.σ./1e9, (ρ_BM4MGD - ρ_exp)' * ρc)
        Colorbar(fig[3,4], hm, label = L"$ρ$ (kg/m³)" )
        
        display(fig)
    end
    with_theme(Visualisation, theme_latexfonts())

end
