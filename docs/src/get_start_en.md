# Getting Started With HydroModels.jl

`HydroModels.jl` is a modern hydrological modeling framework based on the Julia language. Built upon and extending the SUPERFLEX design philosophy, it offers flexible model construction capabilities, efficient computational performance, and support for deep learning model integration. This tutorial will demonstrate how to build and run hydrological models using `HydroModels.jl`.

## Installation

```julia
] add HydroModels
```

## Build a ExpHydro Model

### Introduction to ExpHydro Model

`ExpHydro` is a simple hydrological model consisting of a snow module and a soil module. Its mathematical expressions are as follows.

```math
\begin{aligned}
& \text{Snowpack Bucket:} \\
& pet = 29.8 \cdot lday \cdot 24 \cdot 0.611 \cdot \frac{\exp(17.3 \cdot temp)}{temp + 237.3} \cdot \frac{1}{temp + 273.2} && (1) \\
& snowfall = H(T_{min} - temp) \cdot prcp && (2) \\
& rainfall = H(temp - T_{m,n}) \cdot prcp,&& (3) \\
& melt = H(temp - T_{max}) \cdot H(snowpack) \cdot \min(snowpack, D_f \cdot (temp - T_{max})) && (4) \\
& \frac{d(snowpack)}{dt} = snowfall - melt && (5) \\
\\
& \text{Soilwater Bucket:} \\
& evap = H(soilwater) \cdot pet \cdot \min(1.0, \frac{soilwater}{S_{max}}) && (6) \\
& baseflow = H(soilwater) \cdot Q_{max} \cdot \exp(-f \cdot \max(0.0, S_{max} - soilwater)) && (7) \\
& surfaceflow = \max(0.0, soilwater - S_{max}) && (8) \\
& flow = baseflow + ,urfaceflow && (9) \\
& \frac{d(soilwater)}{dt} = rainfall + melt - evap - flow && (10)
\end{aligned}
```

where $H(x)$ represents the Heaviside step function that equals 1 when $x > 0$ and 0 otherwise; $T_{min}$, $T_{max}$, $D_f$, $S_{max}$, $Q_{max}$, and $f$ are model parameters; $temp$, $lday$, and $prcp$ are input variables; $snowpack$ and $soilwater$ are state variables; and the remaining variables are intermediate calculation variables.

### Build an ExpHydro Model in HydroModels.jl

After understanding the basic principles of the `ExpHydro` model, we can use `HydroModels.jl` to build this model.

#### Import Packages

First, import the HydroModels.jl package:

```julia
using HydroModels
```

#### Define Variables and Parameters

Define the model parameters, state variables, and other variables (including precipitation, evaporation, snowmelt, surface flow, etc.):

```julia
@variables temp lday pet prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f
```

#### Build Flux Formulas by HydroFlux

Next, we need to construct the various computational formulas of the model, including state equations and intermediate variable calculations. Let's use the snowmelt formula as an example to demonstrate how to use `HydroFlux`:

```julia
# define the step function
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

melt_flux = HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))])
```

When constructing `HydroFlux`, several points need attention:

1. Ensure all variables and parameters used in the expression are declared in the inputs
2. Input and output variables are connected through the `Pair` type (=>)
3. Parameters are passed as vectors
4. Expressions are passed as vectors to the `exprs` parameter

`HydroFlux` first accepts a `Pair` type consisting of input and output variables, then takes a `Vector` of parameters, and finally accepts the snowmelt formula expression through `exprs` to complete the construction of the snowmelt formula.

It's important to note that when building `HydroFlux`, you must ensure that all variables and parameters in `expr` are included in the input variables and parameters. If there are variables or parameters in `exprs` that haven't been provided to `HydroFlux`, errors may occur during subsequent calculations due to undefined variables. In such cases, you'll need to check if the `HydroFlux` construction is correct.

Of course, when building HydroFlux, parameter variables are not always necessary. If parameters are not needed, they can be omitted, as shown in the potential evapotranspiration calculation formula:

```julia
pet_flux = HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp(17.3 * temp) / (temp + 237.3) * 1 / (temp + 273.2)])
```

This formula uses `temp` and `lday` as input variables and calculates `pet` as the output variable. Since this formula doesn't involve parameters, there's no need to pass parameter variables.

It's worth noting that `HydroFlux`'s `exprs` parameter can only accept a `Vector` of expressions. This is because `HydroFlux` sometimes needs to support multiple output variables, thus requiring multiple expressions, and the number of expressions must match the number of outputs. Here's an example using the rain-snow separation formula:

```julia
HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp])
```

This formula uses `prcp` and `temp` as input variables and calculates `snowfall` and `rainfall` as output variables, while using `Tmin` as a parameter. Since this formula outputs two variables (`snowfall` and `rainfall`), two expressions must be provided to represent the calculation formulas for each output variable.

#### Build State Formulas by StateFlux

Generally, state equations are constructed as the sum of input fluxes minus the sum of output fluxes. Therefore, we can use a `Pair` to represent this input-output relationship and include the state variable to build a `StateFlux`. Here's an example with the snowpack state equation:

```julia
snowpack_flux = StateFlux([snowfall] => [melt], snowpack)
```

This formula uses `snowfall` as the input variable and calculates `snowpack` as the output variable, while incorporating the state variable `snowpack` to construct the snowpack state equation.

Typically, before constructing state equations, it's preferable to build all involved fluxes using `HydroFlux` first, then combine them using `StateFlux` with input-output flux `Pair`s and state variables to get an intuitive `StateFlux`. However, to avoid introducing unnecessary variables, you can also directly construct more complex `StateFlux` equations. In fact, this construction method is highly consistent with `HydroFlux` - it similarly takes input and output variables as `StateFlux` inputs, names the state variables and required parameters, and finally constructs the state equation. Let's rebuild the snowpack state equation:

```julia
snowpack_flux = StateFlux([prcp, temp, melt], snowpack, [Tmin], expr=step_func(Tmin - temp) * prcp-melt)
```

In this state equation, instead of directly calculating `snowfall`, we perform the calculation using `prcp` and `temp`, combining with the state equation `dsnowpack/dt=snowfall-melt` to construct a complex state equation that omits the `snowfall` intermediate variable.

#### Build Snowfall Bucket

After completing the necessary `HydroFlux` and `StateFlux` constructions, we can build a `HydroBucket`. Using the snowpack `Bucket` as an example, we need to pass the constructed `HydroFlux` and `StateFlux` to `HydroBucket`.

First, store the relevant `HydroFlux` and `StateFlux` in a `Vector`:

```julia
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)
```

Then use these `Vector` variables as input to build a hydrological computation module, i.e., `HydroBucket`:

```julia
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)
```

Through `HydroBucket`, we can complete the construction of a hydrological computation module. During the construction of `HydroBucket`, the program executes some automatic construction methods to express the ordinary differential equations and other variable calculation functions in the computation module. For implementation details, please refer to the implementation section.

#### Build Exphydro Model

Following the same logic, we can build the soil computation module, then concatenate the two modules into a `Vector` type and input it to `HydroModel` to complete the construction of the `ExpHydro` model. The complete construction code is as follows:

```julia
#* import packages
using HydroModels

#* define variables and parameters
@variables temp lday pet prcp 
@variables snowfall rainfall melt evap baseflow surfaceflow flow
@variables snowpack soilwater
@parameters Tmin Tmax Df Smax Qmax f

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

#* define snowpack bucket
fluxes_1 = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
]
dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
snowpack_bucket = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

#* define soilwater bucket
fluxes_2 = [
    HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
    HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
    HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
    HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
]
dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)

#* define the Exp-Hydro model
exphydro_model = HydroModel(name=:exphydro, components=[snowpack_bucket, soilwater_bucket])
```

This completes the construction of an `ExpHydro` model. When we print this type, we can see the basic information contained in the model:

```txt
HydroModel: exphydro
   Components: surface, soil
   Inputs: [temp, lday, prcp]
   States: [snowpack, soilwater]
   Outputs: [pet, snowfall, rainfall, melt, evap, baseflow, surfaceflow, flow]
   Parameters: [Tmin, Tmax, Df, Smax, Qmax, f]
   Components:
       Fluxes: 0 fluxes
       Buckets: 2 buckets
       Routes: 0 routes
```

This information shows the model's name, component names, input variables, state variables, output variables, parameter variables, and the number of `Flux`es, `Bucket`s, and `Route`s contained in the model.

#### Define Neural Network by Lux.jl

`HydroModels.jl` supports neural network definition through `Lux.jl`. Here's an example using a simple fully connected neural network:

```julia
using Lux

# define the ET NN and Q NN
ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)

q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
```

After constructing two fully connected neural networks `ep_nn` and `q_nn` using `Lux.jl`, we can build neural network formulas using `NeuralFlux` provided by `HydroModel.jl`:

```julia
@variables norm_snw norm_slw norm_temp norm_prcp
@variables log_evap_div_lday log_flow

ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)
```

`NeuralFlux` is a derivative type of `HydroFlux`. It similarly accepts a `Pair` of input and output variables, but unlike `HydroFlux`, it doesn't require model parameters and calculation formulas. Instead, it needs the neural network to be passed as a parameter to `NeuralFlux` to construct the neural network formula.
It's important to note that the input and output variable dimensions must match those of the `NeuralFlux`.

After completing the construction of `NeuralFlux`, we can pass it along with other `HydroFlux`es to `HydroBucket` to build a neural network-coupled hydrological model:

```julia
soil_fluxes = [
    #* normalize
    HydroFlux([snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [snowpack_mean, soilwater_mean, prcp_mean, temp_mean, snowpack_std, soilwater_std, prcp_std, temp_std],
        exprs=[(var - mean) / std for (var, mean, std) in zip([snowpack, soilwater, prcp, temp],
            [snowpack_mean, soilwater_mean, prcp_mean, temp_mean],
            [snowpack_std, soilwater_std, prcp_std, temp_std]
        )]),
    ep_nn_flux,
    q_nn_flux
]
```

This soil module includes normalization formulas, evaporation formulas, and runoff formulas. Next, we can pass these formulas to `HydroBucket` and combine them with state equations to complete the construction of the soil module:

```julia
state_expr = rainfall + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
soil_dfluxes = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr)]
soilwater_bucket = HydroBucket(name=:soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)
m50_model = HydroModel(name=:m50, components=[snowpack_bucket, soilwater_bucket])
```

## Run a ExpHydro Model

After completing the model construction, since the model is `callable`, we can input the data and parameters into the model to obtain output results.

### Prepare Input Data

The model accepts input observation data of type `AbstractMatrix`, with dimensions being the number of input variables multiplied by the number of time steps. For example, if there are 3 input variables and 1000 time steps, the input observation data dimensions would be (3, 1000). It's important to note that the index of input variables in the dimension will affect the model calculation results, as the model cannot automatically identify which row in the Matrix corresponds to which input variable. Therefore, users can first prepare a `Dict` or `NamedTuple` type, then use `HydroModels.get_input_names(model)` to get the names of the model's input variables, and combine these names with the `Dict` or `NamedTuple` to construct the `Matrix` type input data:

```julia
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))')
```

### Prepare Parameters and Initial States

For parameter preparation, `HydroModels.jl` accepts parameters of type `ComponentVector` (`using ComponentArrays`), storing parameter names and values in the `ComponentVector`, for example:

```julia
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
```

For initial state preparation, `HydroModels.jl` similarly accepts initial states of type `ComponentVector`, storing state names and values in the `ComponentVector`, for example:

```julia
init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248)
```

Then store the parameters and initial states in a single `ComponentVector` (sometimes initial states are also treated as parameters to be optimized, so they need to be stored together with parameters for unified management):

```julia
pas = ComponentVector(params=params, initstates=init_states)
```

### Run the ExpHydro Model

Finally, we can input the data and parameters into the model to obtain the output results:

```julia
result = exphydro_model(input_arr, pas)
```

The model output is also of type `AbstractMatrix`, with dimensions being the number of output variables (including state variables) multiplied by the number of time steps. For example, if there are 8 output variables and 1000 time steps, the output dimensions would be (8, 1000). If you want to export the results as a `DataFrame` type, you can use `HydroModels.get_output_names(model)` to get the names of the output variables, then combine them with a `Dict` or `NamedTuple` to construct a `DataFrame`:

```julia
states_and_output_names = vcat(HydroModels.get_state_names(exphydro_model), HydroModels.get_output_names(exphydro_model))
output = NamedTuple{Tuple(states_and_output_names)}(eachslice(result, dims=1))
df = DataFrame(output)
```

The calculation results look like this:

```txt
10000×10 DataFrame
   Row │ snowpack  soilwater  pet      snowfall  rainfall  melt      evap      baseflow   surfaceflow  flow      
       │ Float64   Float64    Float64  Float64   Float64   Float64   Float64   Float64    Float64      Float64   
───────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────
     1 │    0.0         0.0      3.1        0.0   1305.22  1.13779   0.868735  0.0212194          0.0  0.0212194
     2 │    0.0         0.0      4.24       0.0   1308.28  1.51019   1.15578   0.0223371          0.0  0.0223371
     3 │    0.0         0.0      8.02       0.0   1315.03  1.63204   1.25547   0.0250095          0.0  0.0250095
     4 │    0.0         0.0     15.27       0.0   1329.34  1.21771   0.946937  0.0317802          0.0  0.0317802
     5 │    0.0         0.0      8.48       0.0   1336.99  1.02779   0.803845  0.0361228          0.0  0.0361228
   ⋮   │    ⋮          ⋮         ⋮        ⋮         ⋮         ⋮         ⋮          ⋮           ⋮           ⋮
  9996 │  254.264       0.0      0.0       -0.0   1521.37  0.221762  0.197362  0.791897           0.0  0.791897
  9997 │  266.324      12.06     0.0       -0.0   1520.37  0.238907  0.21248   0.778688           0.0  0.778688
  9998 │  277.874      11.55     0.0       -0.0   1519.33  0.288362  0.256291  0.765307           0.0  0.765307
  9999 │  279.714       1.84     0.0       -0.0   1518.29  0.311022  0.27624   0.752073           0.0  0.752073
 10000 │  279.854       0.14     0.0       -0.0   1517.38  0.176836  0.156967  0.740711           0.0  0.740711
```

### Plot the Results

The model results are compared with the actual values ​​and plotted into a line graph as follows.

![exphydro model predict](assets/exphydro_predict.png)

### Running with Config

The model computation is not limited to these basic features. To support necessary settings during the computation process, `HydroModel` can accept a `config` parameter.

For example, by default, the model uses only the Euler method to solve ordinary differential equations. However, sometimes a more accurate solving method is needed, which can be set using the `config` parameter.

To demonstrate this functionality, we need to additionally import `OrdinaryDiffEq.jl` and `HydroModelTools.jl`, then set the solving method:

```julia
using OrdinaryDiffEq
using HydroModelTools
using DataInterpolations

config = (solver=ODESolver(alg=Tsit5(), abstol=1e-3, reltol=1e-3), interp=LinearInterpolation)
output = exphydro_model(input_arr, pas, config=config)
```

The `config` parameter is set as a `NamedTuple` type, where the `solver` parameter is used to set the solving method, and the `interp` parameter is used to set the interpolation method. `ODESolver` is a solver wrapper class provided by `HydroModelTools.jl` that performs some data conversion on the solution results. It accepts solving methods from `OrdinaryDiffEq.jl` as parameters and sets the solving method parameters.

[Would you like me to continue with the neural network implementation section?]

### Run the Neural Network Enhanced Model

Earlier we built an `ExpHydro` model coupled with neural networks, and now we can use this model to calculate results.

First, we need to prepare the model parameters. In addition to the `ExpHydro` model parameters, we also need to prepare neural network parameters:

```julia
using StableRNGs

ep_nn_params = Vector(ComponentVector(LuxCore.initialparameters(ep_nn)))
q_nn_params = Vector(ComponentVector(LuxCore.initialparameters(q_nn)))
```

Models based on `Lux.jl` can decouple parameters from model structure, which is the primary reason we chose to use `Lux.jl` to build neural networks. After using `LuxCore.initialparameters`, we convert the parameters (`NamedTuple`) to `Vector` type, then combine them with the model parameters and initial states to build the complete model parameters:

```julia
nn_pas = ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
pas = ComponentVector(params=params, initstates=init_states, nns=nn_pas)
```

It's important to note that we need to use the key name `nns` to store the neural network parameters, keeping them independent from model parameters and initial states. Additionally, the neural network parameters `nn_pas` need to use model names and parameters as key-value pairs.

Finally, we can input the model parameters into the model to obtain the output results:

```julia
output = m50_model(input_arr, pas)
```

## Conclusion

This tutorial has provided a detailed introduction to the basic usage of `HydroModels.jl`, including:
- Model building process
- Parameter and initial state configuration
- Model execution and computation

Through this tutorial, users can quickly master the core functionality of `HydroModels.jl` and begin building their own hydrological models.

`HydroModels.jl` also offers more advanced features, including:
- Semi-distributed and distributed hydrological model construction and computation
- Neural network-coupled parameter-adaptive models
- `GR4J` model with unit hydrograph routing
- Model parameter optimization
- Model extension functionality based on `Wrapper`

These advanced features will be covered in subsequent tutorials.
