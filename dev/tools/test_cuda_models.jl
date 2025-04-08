# 导入模块
using CSV
using DataFrames
using ComponentArrays
using CUDA
using LuxCUDA
using cuDNN
using HydroModels
using Lux
using BenchmarkTools
using Pipe: @pipe

# dev = gpu_device()
dev = identity

include("../models/exphydro.jl")
include("../src/HydroModelTools.jl")

# define parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))')

# run model with multi node input
node_num = 10
node_params = ComponentVector(
    f=fill(f, node_num), Smax=fill(Smax, node_num), Qmax=fill(Qmax, node_num),
    Df=fill(Df, node_num), Tmax=fill(Tmax, node_num), Tmin=fill(Tmin, node_num)
) |> dev
node_states = ComponentVector(snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num)) |> dev

node_pas = ComponentVector(params=node_params) |> dev
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], node_num))
node_input = permutedims(node_input, (2, 3, 1)) |> dev
config = Dict(:ptyidx => 1:node_num, :styidx => 1:node_num, :timeidx => ts, :solver => HydroModelTools.DiscreteSolver(dev=dev), :device => dev)
result = exphydro_model(node_input, node_pas, initstates=node_states, config=config)
ntp_pre, ntp_post = HydroModelTools.NamedTuplePreprocessor(exphydro_model), HydroModelTools.NamedTuplePostprocessor(exphydro_model)
node_input_ntp = fill(input, node_num)
input_tuple = (node_input_ntp, node_pas, Dict(:initstates=>node_states, :config=>config))
@pipe input_tuple |> ntp_pre(_[1], _[2], _[3]) |> exphydro_model(_[1], _[2]; _[3]...) |> ntp_post(_)

# exphydro_model.components[2].flux_funcs[2](eachslice(view(node_input, exphydro_model.varindices[1], :, :), dims=1), eachslice(node_input[1:1, :, :], dims=1), node_pas)