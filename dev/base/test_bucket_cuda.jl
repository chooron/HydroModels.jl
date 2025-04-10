# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using BenchmarkTools
using CUDA, cuDNN
using Lux
using DataInterpolations
using HydroModels

include("../models/exphydro.jl")

gpu = gpu_device()
ele = bucket_1

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])


# multi node input
node_num = 10
node_names = [Symbol(:node, i) for i in 1:node_num]
node_params = ComponentVector(
    f=fill(0.0167, node_num), Smax=fill(1709.46, node_num), Qmax=fill(18.47, node_num),
    Df=fill(2.6745, node_num), Tmax=fill(0.1757, node_num), Tmin=fill(-2.093, node_num)
) |> gpu
node_states = ComponentVector(
    snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num)
) |> gpu

node_pas = ComponentVector(params=node_params)
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(ele)])) |> gpu
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
node_input = permutedims(node_input, (2, 3, 1))
config = (ptyidx=1:10, styidx=1:10, timeidx=ts, solver=HydroModels.ManualSolver(mutable=false, dev=gpu))
result = ele(node_input, node_pas, config=config)
