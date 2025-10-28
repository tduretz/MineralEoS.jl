# MineralEoS.jl

![CI](https://github.com/tduretz/MineralEoS.jl/actions/workflows/UnitTests.yml/badge.svg)

This package is designed to compute the densities and volumes of the following mineral phases:
- OlivineFo90
- Diamond

`MineralEoS.jl` does aim to provide all the features of [`EoSFit`](http://www.rossangel.com/text_eosfit.htm); rather, it focuses on exposing core functionalities, such as density computation, and making them available for codes that rely on automatic differentiation in Julia (e.g. [StagFDTools](https://github.com/tduretz/StagFDTools)).

`MineralEoS.jl` was initially developped by T. Duretz & R.J. Angel.

# Quickstart

The following example:
```julia
using MineralEoS
params = assign_EoS_parameters(:Diamond)
P      = 5.0e9
T      = 1100.0
ρ, V   = density_volume(P, T, params)
```
should return the density (kg/m^3) and the volume (cm^3):
```julia-repl
julia> (3526.1420454161403, 3.4091933464867825)
```

The following example uses an additional dimensional scaling:
```julia
using MineralEoS
scales = (σ = 1e9, L = 1e3, t = 1e12, T = 100.0)
params = assign_EoS_parameters(:Diamond, sc=scales)
P      = 5.0e9  / scales.σ
T      = 1100.0 / scales.T
ρ, V   = density_volume(P, T, params)
ρc     = scales.σ * scales.L * scales.t^2.0 / scales.L^3 
(ρ * ρc, V * scales.L^3)
```
should also return the density (kg/m^3) and the volume (cm^3):
```julia-repl
julia> (3526.1420454161394, 3.4091933464867834)
```

# Governing equations
The total pressure accounts for mechanical and thermal contributions:

```math
P = P_\text{ref}(V, T_\text{0}) + \Delta P_\text{th}(V, T)  
```

This equation is non-linear and can be solved for $V$ by Newton-Raphson iterations. The latter are greatly facilitated by the use of the automatic differention package [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl).  

The mechanical part accounts for the Birch-Murnaghan model:\
$$ P_\text{ref} = 3K_0 f (1+2f)^\frac{5}{2}  \left(1 + f\frac{3}{2}\left(\frac{\partial K}{\partial P} -4\right) + f^2 \frac{3}{2}\left(K_0\frac{\partial^2 K}{\partial P^2} + \left(\frac{\partial K}{\partial P} - 4\right) \left(\frac{\partial K}{\partial P}-3\right) + \frac{35}{9}\right) \right)$$
where $f = \frac{1}{2}\left(\left(\frac{V_0}{V}\right)^\frac{2}{3} - 1.0\right)$ and $K$ is the bulk modulus.

The thermal is based on the Einstein model:\
$$ \Delta P_\text{th} = 3 N \frac{γ}{V}\left(U(T) - U(T_\text{0})\right)$$
where $N$ is a number of atoms, $U(T) = \frac{R \theta_\text{E}}{\exp{\frac{\theta_\text{E}}{T}} - 1}$ and $\gamma = γ_0  \left(\frac{V}{V_0}\right)^q$. The models can also represent a q-compromised equation of states, provided $q = 1$. 

The parameters are defined within the database of [`EoSFit`](http://www.rossangel.com/text_eosfit.htm) and the results are consistent with those of `EoSFit`.

# Examples

An example for Diamond is available [here](/example/Diamond.jl) and gives the result:
![](/results/Diamond.png)

A similar example using dimensional scaling is available [here](/example/Diamond_scaled.jl).

An example for OlivineF90 is available [here](/example/OlivineF90.jl) and gives the result:
![](/results/OlivineF90.png)

