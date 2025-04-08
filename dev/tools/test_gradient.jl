using Zygote
using ComponentArrays
using HydroModels
using Lux
using StableRNGs
using CSV, DataFrames
include("../src/HydroModelTools.jl")

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
@variables melt log_evap_div_lday log_flow
@variables norm_snw norm_slw norm_temp norm_prcp

#! define the snow pack reservoir
snow_fluxes = [
    HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
    HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
    HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))]),
]
snow_dfluxes = [StateFlux([snowfall] => [melt], snowpack)]
snow_ele = HydroBucket(name=:exphydro_snow, fluxes=snow_fluxes, dfluxes=snow_dfluxes)

#! define the ET NN and Q NN
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

#! get init parameters for each NN
ep_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], ep_nn)
q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

#! define the soil water reservoir
soil_fluxes = [
    #* normalize
    HydroFlux([snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
        [snowpack_mean, soilwater_mean, prcp_mean, temp_mean, snowpack_std, soilwater_std, prcp_std, temp_std],
        exprs=[(var - mean) / std for (var, mean, std) in zip([snowpack, soilwater, prcp, temp],
            [snowpack_mean, soilwater_mean, prcp_mean, temp_mean],
            [snowpack_std, soilwater_std, prcp_std, temp_std]
        )]),
    ep_nn_flux,
    q_nn_flux,
]

state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
soil_dfluxes = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, expr=state_expr)]
soil_ele = HydroBucket(name=:m50_soil, fluxes=soil_fluxes, dfluxes=soil_dfluxes)
#! define the Exp-Hydro model
m50_model = HydroModel(name=:m50, components=[snow_ele, soil_ele])
ioadapter = NamedTupleIOAdapter(m50_model)

# predefine the parameters
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
# load data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:100)
# cols: Baseflow,Evap,Flow,Infiltration,Lday,Melt,Pet,Prcp,Rainfall,Snowfall,Surfaceflow,Temp,SoilWater,SnowWater
lday_vec = df[ts, "Lday"]
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
flow_vec = df[ts, "Flow"]

log_flow_vec = log.(flow_vec)
log_evap_div_lday_vec = log.(df[ts, "Evap"] ./ lday_vec)
norm_prcp_vec = (prcp_vec .- mean(prcp_vec)) ./ std(prcp_vec)
norm_temp_vec = (temp_vec .- mean(temp_vec)) ./ std(temp_vec)
norm_snw_vec = (df[ts, "SnowWater"] .- mean(df[ts, "SnowWater"])) ./ std(df[ts, "SnowWater"])
norm_slw_vec = (df[ts, "SoilWater"] .- mean(df[ts, "SoilWater"])) ./ std(df[ts, "SoilWater"])
nn_input = (norm_snw=norm_snw_vec, norm_slw=norm_slw_vec, norm_temp=norm_temp_vec, norm_prcp=norm_prcp_vec)

m50_opt = HydroModelTools.GradOptimizer(component=NamedTupleIOAdapter(m50_model), solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote(), maxiters=10)
config = (solver=HydroModelTools.ODESolver(sensealg=BacksolveAdjoint(autodiff=ZygoteVJP())),)
norm_pas = ComponentVector(
    snowpack_mean=mean(norm_snw_vec), soilwater_mean=mean(norm_slw_vec), prcp_mean=mean(norm_prcp_vec), temp_mean=mean(norm_temp_vec),
    snowpack_std=std(norm_snw_vec), soilwater_std=std(norm_slw_vec), prcp_std=std(norm_prcp_vec), temp_std=std(norm_temp_vec)
)
m50_const_pas = ComponentVector(
    initstates=ComponentVector(snowpack=0.0, soilwater=1300.0),
    params=ComponentVector(Tmin=Tmin, Tmax=Tmax, Df=Df; norm_pas...)
)
m50_tunable_pas = ComponentVector(
    nn=ComponentVector(epnn=ep_nn_params, qnn=q_nn_params)
)
m50_input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
q_output = (log_flow=flow_vec,)
m50_opt_params, m50_loss_df = m50_opt(
    [m50_input], [q_output],
    tunable_pas=m50_tunable_pas,
    const_pas=m50_const_pas,
    config=[config],
    return_loss_df=true
)