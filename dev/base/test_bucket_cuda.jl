# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using BenchmarkTools
using CUDA, cuDNN
using Lux
using DataInterpolations

include("../src/HydroModels.jl")

HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/exphydro.jl")

gpu = gpu_device()
ele = bucket_1

HydroModels.get_param_names(ele)
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
# 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
# params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
# init_states = ComponentVector(snowpack=0.0, soilwater=1303.004248) |> gpu
# pas = ComponentVector(params=params, initstates=init_states) |> gpu

pas = ComponentVector(params=[-2.09, 0.17, 18.46], initstates=[0.0, 1303.004248]) |> gpu

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = CuArray(Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))'))
results = ele(input_arr, pas)
ele.ode_func([2,3,4], [10], [1,2,3], [])
ele.flux_func([2,3,4], [10], [1,2,3], [])

# # multi node input
# node_num = 10
# node_names = [Symbol(:node, i) for i in 1:node_num]
# node_params = ComponentVector(
#     f=fill(f, node_num), Smax=fill(Smax, node_num), Qmax=fill(Qmax, node_num),
#     Df=fill(Df, node_num), Tmax=fill(Tmax, node_num), Tmin=fill(Tmin, node_num)
# )
# node_states = ComponentVector(
#     snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num)
# )

# node_pas = ComponentVector(params=node_params[HydroModels.get_param_names(ele)], initstates=node_states[HydroModels.get_state_names(ele)])
# input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)]))
# node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
# node_input = permutedims(node_input, (2, 3, 1))
# config = (ptyidx=1:10, styidx=1:10, timeidx=ts)
# result = ele(node_input, node_pas, config=config)
