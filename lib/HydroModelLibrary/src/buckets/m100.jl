@reexport module M100
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
@variables melt log_evap_div_lday log_flow asinh_melt asinh_ps asinh_pr
@variables norm_snw norm_slw norm_temp norm_prcp

#! define the m100 NN
m100_nn = Lux.Chain(
    Lux.Dense(4, 32, tanh),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 5),
    name=:nn
)

#! get init parameters for each NN
#! define the soil water reservoir
m100_funcs = [
    #* normalize
    HydroFlux([snowpack] => [norm_snw], [snowpack_mean, snowpack_std], exprs=[(snowpack - snowpack_mean) / snowpack_std]),
    HydroFlux([soilwater] => [norm_slw], [soilwater_mean, soilwater_std], exprs=[(soilwater - soilwater_mean) / soilwater_std]),
    HydroFlux([prcp] => [norm_prcp], [prcp_mean, prcp_std], exprs=[(prcp - prcp_mean) / prcp_std]),
    HydroFlux([temp] => [norm_temp], [temp_mean, temp_std], exprs=[(temp - temp_mean) / temp_std]),
    NeuralFlux([norm_snw, norm_slw, norm_prcp, norm_temp] => [log_evap_div_lday, log_flow, asinh_melt, asinh_ps, asinh_pr], m100_nn),
    HydroFlux([asinh_melt, snowpack] => [melt], exprs=[relu(sinh(asinh_melt) * step_func(snowpack))]),
]

state_expr1 = relu(sinh(asinh_ps)) * step_func(-temp) - melt
state_expr2 = relu(sinh(asinh_pr)) + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
m100_dfuncs = [
    StateFlux([asinh_ps, temp, melt], snowpack, expr=state_expr1),
    StateFlux([asinh_pr, melt, soilwater, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr2),
]
m100_bucket = HydroBucket(name=:m100_bucket, funcs=m100_funcs, dfuncs=m100_dfuncs)
#! define the Exp-Hydro model
m100_model = HydroModel(name=:m100, components=[m100_bucket])

export m100_bucket, m100_model

end