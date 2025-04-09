# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using BenchmarkTools
using DataInterpolations
using Zygote
using ForwardDiff

include("../src/HydroModels.jl")

HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/exphydro.jl")

ele = bucket_1
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
# 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)[HydroModels.get_param_names(ele)]
init_states = ComponentVector(snowpack=100.0)[HydroModels.get_state_names(ele)]
pas = ComponentVector(params=params, initstates=init_states)
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)])) |> permutedims
single_ele = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)
# result = single_ele(input_arr, pas)

# multi node input
node_num = 10
mul_ele = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)
node_names = [Symbol(:node, i) for i in 1:node_num]
node_params = ComponentVector(
    f=fill(f, node_num), Smax=fill(Smax, node_num), Qmax=fill(Qmax, node_num),
    Df=fill(Df, node_num), Tmax=fill(Tmax, node_num), Tmin=fill(Tmin, node_num)
)[HydroModels.get_param_names(ele)]
node_states = ComponentVector(
    snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num)
)
node_pas = ComponentVector(params=node_params)
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
node_input = permutedims(node_input, (2, 3, 1))
config = (ptyidx=1:10, styidx=1:10, timeidx=ts, interp=LinearInterpolation, solver=HydroModels.ManualSolver(mutable=true))
result = mul_ele(node_input, node_pas; initstates=node_states, config...)

# ForwardDiff.gradient(p -> mul_ele(node_input, p; initstates=node_states, config...)[:,:,:] |> sum, node_pas)

