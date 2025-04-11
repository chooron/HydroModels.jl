# 导入模块
using CSV
using DataFrames
using ComponentArrays
using HydroModels
using DataInterpolations

include("../models/exphydro.jl")

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

bucket_input_names = HydroModels.get_input_names(bucket_1)
bucket_param_names = HydroModels.get_param_names(bucket_1)
bucket_state_names = HydroModels.get_state_names(bucket_1)
# 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)[bucket_param_names]
init_states = ComponentVector(snowpack=100.0)[bucket_state_names]
pas = ComponentVector(params=params)
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = reduce(hcat, input[bucket_input_names]) |> permutedims
result = bucket_1(input_arr, pas)

# multi node input
node_num = 10
node_names = [Symbol(:node, i) for i in 1:node_num]
node_params = ComponentVector(
    f=fill(f, node_num), Smax=fill(Smax, node_num), Qmax=fill(Qmax, node_num),
    Df=fill(Df, node_num), Tmax=fill(Tmax, node_num), Tmin=fill(Tmin, node_num)
)
node_states = ComponentVector(
    snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num)
)
node_pas = ComponentVector(params=node_params)
input_arr = reduce(hcat, collect(input[bucket_input_names]))
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
node_input = permutedims(node_input, (2, 3, 1))
config = (ptyidx=1:10, styidx=1:10, timeidx=ts, interp=LinearInterpolation, solver=HydroModelTools.ODESolver())
result = bucket_1(node_input, node_pas; initstates=node_states, config...)