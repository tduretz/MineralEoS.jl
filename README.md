# MineralEoS.jl

![CI](https://github.com/tduretz/MineralEoS.jl/actions/workflows/UnitTests.yml/badge.svg)

This package is designed to compute the densities and volumes of the following mineral phases:
- OlivineFo90
- Diamond

`MineralEoS.jl` does not aim to provide all the features of [`EoSFit`](http://www.rossangel.com/text_eosfit.htm); rather, it focuses on exposing core functionalities, such as density computation, and making them available for codes that rely on automatic differentiation in Julia (e.g. [StagFDTools](https://github.com/tduretz/StagFDTools)).

`MineralEoS.jl` was initially developped by T. Duretz & R.J. Angel.

# Quickstart

The following example:
```julia
using MineralEoS
EoS    = ComplexEoS()
params = assign_EoS_parameters(:Diamond)
P      = 5.0e9
T      = 1100.0
ρ, V   = density_volume(EoS, P, T, params)
```
should return the density (kg/m<sup>3</sup>) and the volume (cm<sup>3</sup>):
```julia-repl
julia> (3526.1420450764394, 3.4091933468152167)
```

The following example uses an additional dimensional scaling:
```julia
using MineralEoS
scales = (σ = 1e9, L = 1e3, t = 1e12, T = 100.0)
EoS    = ComplexEoS()
params = assign_EoS_parameters(:Diamond, sc=scales)
P      = 5.0e9  / scales.σ
T      = 1100.0 / scales.T
ρ, V   = density_volume(EoS, P, T, params)
ρc     = scales.σ * scales.L * scales.t^2.0 / scales.L^3 
(ρ * ρc, V * scales.L^3)
```
should also return the density (kg/m<sup>3</sup>) and the volume (cm<sup>3</sup>):
```julia-repl
julia> (3526.1420450764394, 3.409193346815217)
```

Using simple and complex EoS, with variable degree of complexity:

```julia
using MineralEoS
params = assign_EoS_parameters(:Diamond)
P      = 5.0e9
T      = 1100.0
# Uses exponential model
ρ, V   = density_volume(SimpleEoS(), P, T, params)
# Complex model. Default Birch-Munaghan O3 + Einstein
ρ, V   = density_volume(ComplexEoS(), P, T, params)
# More flexibility...
opts   = (thermal_model=Debye(), mechanical_model=BM4())
ρ, V   = density_volume(ComplexEoS(), P, T, params; options=opts)
```

# Governing equations

The total pressure accounts for mechanical and thermal contributions:

```math
P = P_\text{ref}(V, T_\text{0}) + \Delta P_\text{th}(V, T)  
```

This equation is non-linear and can be solved for $V$ by Newton-Raphson iterations. The latter are greatly facilitated by the use of the automatic differention package [Enzyme.jl](https://github.com/EnzymeAD/Enzyme.jl). This EoS model is referred to as `:complex` (Birch-Murnaghan, Einstein, Mie-Grüneisen-Debye models) or `:simple` (exponential model).

## Birch-Murnaghan (order 2-3-4)

The mechanical part accounts for the Birch-Murnaghan model. The 4th order model (`mechanical_model = :BM4`) is expressed as:

```math
P_\text{ref} = 3K_0 f (1+2f)^\frac{5}{2}  \left(1 + f\frac{3}{2}\left(\frac{\partial K}{\partial P} -4\right) + f^2 \frac{3}{2}\left(K_0\frac{\partial^2 K}{\partial P^2} + \left(\frac{\partial K}{\partial P} - 4\right) \left(\frac{\partial K}{\partial P}-3\right) + \frac{35}{9}\right) \right)
```

where $f = \frac{1}{2}\left(\left(\frac{V_0}{V}\right)^\frac{2}{3} - 1.0\right)$ and $K$ is the bulk modulus.

The 3rd order model (`mechanical_model = BM3`) sees no contribution from either $\frac{\partial^2 K}{\partial P^2}$ nor $\frac{\partial K}{\partial P}^2$:

```math
P_\text{ref} = 3K_0 f (1+2f)^\frac{5}{2}  \left(1 + f\frac{3}{2}\left(\frac{\partial K}{\partial P} -4\right) \right)
```

The 2rd order model (`mechanical_model = :BM2`) sees no contribution from either $\frac{\partial K}{\partial P}$ nor $\frac{\partial^2 K}{\partial P^2}$:

```math
P_\text{ref} = 3K_0 f (1+2f)^\frac{5}{2}
```

# Thermal pressure
The thermal pressure model is expressed as:

```math
\Delta P_\text{th} = 3 N \frac{γ}{V}\left(U(T) - U(T_\text{0})\right)
```

where $N$ is a number of atoms and  $\gamma = γ_0  \left(\frac{V}{V_0}\right)^q$. The models can also represent a q-compromised equation of states, provided $q = 1$ and $\theta_\text{E} = \theta_\text{E0}$.

The Einstein model (`thermal_model = :Einstein`) relies on the following expression for the vibrational energy:
$U(T) = \frac{R \theta_\text{E}}{\exp{\frac{\theta_\text{E}}{T}} - 1}$
The Einstein temperature is expressed as $\theta_\text{E} = \theta_\text{E0}\left(\frac{V}{V_0} \right)^{-\gamma}$.

The Debye model (`thermal_model = :Debye`) is based on the following expression for the vibrational energy:
$U(T) = R T D\left(\frac{\theta_\text{D}}{T}\right)$ where $D(x)$ is the 3rd Debye function. The Debye temperature is expressed as $\theta_\text{D} = \theta_\text{D0}  \exp{\left(- (\frac{\gamma_0}{q})  \left[\left(\frac{V}{V_0}\right)^q - 1.0\right]\right)}$
Note that the Debye model requires the evaluation of an integral, which is performed using a power series summation.

The parameters are defined within the database of [`EoSFit`](http://www.rossangel.com/text_eosfit.htm) and the results are consistent with those of `EoSFit (see below).

## Exponential model

The exponential model typically used in geodynamic simulation is also available, and follows the expression:

 ```math
\rho    = \rho_0 \exp{\left( \frac{P}{K} - \alpha T \right)}
 ```

 where $\rho_0$ is a reference density. This EoS model is referred to as `exp`.


# Examples

## Accuracy

An example for Diamond is available [here](/example/Diamond.jl) and gives the result:
![](/results/Diamond.png)

A similar example using dimensional scaling is available [here](/example/Diamond_scaled.jl).

An example for OlivineF90 is available [here](/example/OlivineF90.jl) and gives the result:
![](/results/OlivineF90.png)

## Complex (`:complex`) versus simplified (`:simple`) EoS

The possibility to switch between `:simple` and `:complex` models is shown in the example [here](/example/Diamond_exp_BM3.jl) and [there](/example/OlivineF90_exp_BM3.jl). An example of results for diamond, and potential differences between the two approaches can be seen below.

The diamond case:
![](/results/Diamond_exp_BM3.png)

This example clearly highlights the fact that thermal effects are not accurately captured by the exponential model. The resulting density diffences between the two modes reaches 20 kg/m<sup>3</sup> at high temperature. 

The olivine case:
![](/results/OlivineF90_exp_BM3.png)

In this case, the density differences are lower. As in the case of the diamond, the exponential models produces a density excess of about 20 kg/m<sup>3</sup> at high temperature. With increasing pressure, this differences vanishes.

## Olivine: Demo of all complex models (`:BM2`, `:BM3`, `:BM4`, `:Einstein`, `:Debye`) versus simplified (`:simple`) EoS

It is also possible to compare the different Birch-Murnaghan models and thermal pressure models.  This is shown the example [here](/example/Olivine_demo.jl) and the results are depicted below.
The results of the complex EoS are compared to the simple exponential models and the density difference relative to the simple model is computed.


![](/results/OlivineF90_demo.png)

The results for `:BM3` and `:BM4` are similar because of the currently defined olivine parameter (implied values). There is no notable difference between Einstein and Debye thermal pressure models.


