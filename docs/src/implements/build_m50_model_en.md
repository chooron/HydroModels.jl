# Embedding Neural Networks in ExpHydro Model: Implementation of [M50](https://hess.copernicus.org/articles/26/5085/2022/) Model

[marv-in](https://github.com/marv-in) proposed the M50 and M100 models in their paper, with code stored in [HydroNODE](https://github.com/marv-in/HydroNODE), pioneering the integration of neural networks into hydrological models (notably by replacing ExpHydro model's calculation formulas while participating in balance equation solving). This model inspired our repository, and recognizing the current complexity in model development, our primary goal is to simplify the process of embedding neural networks into hydrological models. Below, we'll compare the source code with our repository's implementation of M50 and M100 models, showcasing a new approach to model construction.

## M50 Implementation in HydroNODE

First, let's examine the implementation of M50 and M100 models in HydroNODE:

```julia
# define the ODE problem
function NeuralODE_M50_core!(dS,S,p,t)

    Tmin, Tmax, Df = (p_bucket_precal...,)

    Lday = itp_Lday(t)
    P    = itp_P(t)
    T    = itp_T(t)

    g_ET = ann_ET([norm_S0(S[1]), norm_S1(S[2]), norm_T(T)],p[:p1]) #p[idcs_params_ann_ET])
    g_Q = ann_Q([norm_S1(S[2]), norm_P(P)],p[:p2]) #p[idcs_params_ann_Q])

    melting = M(S[1], T, Df, Tmax)
    dS[1] = Ps(P, T, Tmin) - melting
    dS[2] = Pr(P, T, Tmin) + melting - step_fct(S[2])*Lday*exp(g_ET[1])- step_fct(S[2])*exp(g_Q[1])

end
# build and solve ODE problem
prob = ODEProblem(NeuralODE_M50_core!, S_init, Float64.((t_out[1], maximum(t_out))), p)
sol = solve(prob, BS3(), dt=1.0, saveat=t_out, reltol=1e-3, abstol=1e-3, sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()))
# calculate Qout
P_interp = norm_P.(itp_P.(t_out))
S1_ = norm_S1.(sol[2,:])
Qout_ =  exp.(ann_Q(permutedims([S1_ P_interp]),p[:p2])[1,:])
```

In the HydroNODE code, we can see that model construction consists of three steps:

- ODE function definition
- ODE problem construction and solution
- Calculation of Qout from intermediate states

While this code follows the standard style of the DifferentialEquations.jl library, it requires manual construction of model calculation formulas (`M`, `Ps`, `Pr`) and their combination with state equation expressions (`dS`). The model has high coupling, and obtaining intermediate calculation fluxes like `Qout` requires additional calls to corresponding calculation formulas, making the code relatively complex to write.

## M50 Implementation in HydroModels.jl

First, we need to define the parameters and variables designed in the M50 model:

```julia
using HydroModels

#! parameters in the Exp-Hydro model
@parameters Tmin Tmax Df Smax f Qmax
#! parameters in normalize flux
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

#! hydrological flux in the Exp-Hydro model
@variables prcp temp lday pet rainfall snowfall
@variables snowpack soilwater lday pet
@variables melt log_evap_div_lday log_flow flow
@variables norm_snw norm_slw norm_temp norm_prcp

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
```

Compared to conventional conceptual hydrological models, the M50 model uses neural networks to replace hydrological flux calculation formulas (`ET` and `Q`). These neural networks differ significantly from standard calculation formulas, so HydroModels.jl uses `NeuralFlux` to represent neural networks. The NeuralFlux for ET and Q predictions is shown below:

```julia
using Lux
using StableRNGs

# define the ET NN and Q NN
ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)
ep_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), ep_nn))))
q_nn = Lux.Chain(
    Lux.Dense(2 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:qnn
)
q_nn_params = Vector(ComponentVector(first(Lux.setup(StableRNGs.LehmerRNG(1234), q_nn))))


ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)
```

A special aspect of building the M50 model is the need to wrap embedded neural networks using the NeuralFlux model, providing input-output information for the neural networks and converting them into expressions during construction.

Then we can build the remaining parts of the M50 model:

```julia
# define the snow pack reservoir
snow_fluxes = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfluxes = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(name=:exphydro_snow, fluxes=snow_fluxes, dfluxes=snow_dfluxes)

# define the soil water reservoir
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
state_expr = rainfall + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
soil_dfluxes = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr)]
soil_ele = HydroBucket(name=:m50_soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)
convert_flux = HydroFlux([log_flow] => [flow], exprs=[exp(log_flow)])
# define the Exp-Hydro model
m50_model = HydroModel(name=:m50, components=[snow_ele, soil_ele, convert_flux]);
```

In the code above, we present a complete implementation of the M50 model, starting with the Snowpack Bucket implementation from ExpHydro, followed by the Soilwater Bucket implementation with neural network embedding. In the Soilwater Bucket module implementation, we can use defined NeuralFlux (`ep_nn_flux` and `q_nn_flux`) to express embedded neural networks, combining them with other `HydroFlux` components and using `StateFlux` to construct the Soilwater Bucket.
