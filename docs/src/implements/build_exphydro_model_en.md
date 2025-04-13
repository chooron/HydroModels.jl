# Build an ExpHydro Model

## Introduction to ExpHydro Model

This section demonstrates how to build a simple hydrological model - the ExpHydro model - using HydroModels.jl, serving as an introduction to conceptual model construction.
The ExpHydro model consists of two computational modules: the Snowpack Bucket and the Soilwater Bucket. Their mathematical formulations are as follows:

```math
\begin{aligned}
& \text{Snowpack Bucket:} \\
& pet = 29.8 \cdot lday \cdot 24 \cdot 0.611 \cdot \frac{\exp(17.3 \cdot temp)}{temp + 237.3} \cdot \frac{1}{temp + 273.2} && (1) \\
& snowfall = H(T_{min} - temp) \cdot prcp && (2) \\
& rainfall = H(temp - T_{min}) \cdot prcp && (3) \\
& melt = H(temp - T_{max}) \cdot H(snowpack) \cdot \min(snowpack, D_f \cdot (temp - T_{max})) && (4) \\
& \frac{d(snowpack)}{dt} = snowfall - melt && (5) \\
\\
& \text{Soilwater Bucket:} \\
& evap = H(soilwater) \cdot pet \cdot \min(1.0, \frac{soilwater}{S_{max}}) && (6) \\
& baseflow = H(soilwater) \cdot Q_{max} \cdot \exp(-f \cdot \max(0.0, S_{max} - soilwater)) && (7) \\
& surfaceflow = \max(0.0, soilwater - S_{max}) && (8) \\
& flow = baseflow + surfaceflow && (9) \\
& \frac{d(soilwater)}{dt} = rainfall + melt - evap - flow && (10)
\end{aligned}
```

Where: $H(x)$ represents the Heaviside step function, equals 1 when $x > 0$, otherwise 0; $T_{min}, T_{max}, D_f, S_{max}, Q_{max}, f$ are model parameters;$temp, lday, prcp$ are input variables;$snowpack, soilwater$ are state variables;Other variables are intermediate calculation variables

## Complete Model Construction Process

```julia
# import packages
using HydroModels

# define variables and parameters
@variables temp lday pet prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# define snowpack bucket
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

# define soilwater bucket
fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)

# define the Exp-Hydro model
exphydro_model = HydroModel(name=:exphydro, components=[snowpack_bucket, soilwater_bucket])
```

## Step-by-Step Analysis

Let's break down the model construction process into detailed steps.

First, we need to import the HydroModels.jl dependency:

```julia
using HydroModels
```

By importing the HydroModels module, we gain direct access to types defined in HydroFlux, StateFlux, HydroBucket, and HydroModel modules. Additionally, the package exports macros from ModelingToolkit.jl such as @variables and @parameters, which are used for defining and manipulating variables and parameters in hydrological models:

```julia
@variables temp lday prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow pet
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f
```

In this code segment, we define all variables used in the ExpHydro model, including:
- Input variables: temperature (`temp`), day length (`lday`), and precipitation (`prcp`)
- Intermediate calculation variables: `pet`, `snowfall`, `rainfall`, `melt`, `evap`, `baseflow`, `surfaceflow`, `flow`
- State variables: `snowpack` and `soilwater`
- Model parameters: `Tmin`, `Tmax`, `Df`, `Smax`, `Qmax`, and `f`

The definition of `HydroFlux` requires determining input/output variables and model parameters based on calculation formulas. For example, in the rain-snow partitioning formula, the input variables are `prcp` and `temp`, output variables are `snowfall` and `rainfall`, and the model parameter is `Tmin`. The formula translation to `HydroFlux` looks like this:

```julia
# define the smooth function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
# define the rain-snow partitioning flux
split_flux = HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp])
```

The definition of `StateFlux` is based on the balance equation of state variables. Generally, the balance equation equates the rate of change of the state variable to the difference between input and output fluxes. For example, the balance equation for `snowpack` is the difference between `snowfall` and `melt`. The formula translation to `StateFlux` looks like this:

```julia
snowpack_dflux = StateFlux([snowfall] => [melt], snowpack)
```

Next, we define the model calculation formulas using `HydroFlux` and `StateFlux`, and integrate them into `HydroBucket` to create the Snowpack Bucket and Soilwater Bucket components:

```julia
# define snowpack bucket
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

# define soilwater bucket
fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)
```

The construction of `HydroBucket` consists of `HydroFlux` and `StateFlux`, where `HydroFlux` defines the model's calculation formulas, and `StateFlux` defines the balance equations for state variables.

Finally, we combine the `HydroBucket` components into a `HydroModel` to create the complete ExpHydro model:

```julia
# define the Exp-Hydro model
exphydro_model = HydroModel(name=:exphydro, components=[snowpack_bucket, soilwater_bucket])
```
