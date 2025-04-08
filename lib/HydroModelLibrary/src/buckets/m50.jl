@reexport module M50
using ..HydroModels
using ..Lux

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
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

#! define the snow pack reservoir
snow_funcs = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(name=:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

#! define the ET NN and Q NN
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

#! get init parameters for each NN
ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

#! define the soil water reservoir
soil_funcs = [
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
soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr)]
soil_ele = HydroBucket(name=:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)
convert_flux = HydroFlux([log_flow] => [flow], exprs=[exp(log_flow)])
#! define the Exp-Hydro model
m50_model = HydroModel(name=:m50, components=[snow_ele, soil_ele, convert_flux])

export m50_model, eq_nn_flux, q_nn_flux
end