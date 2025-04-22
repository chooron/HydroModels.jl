# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using DataInterpolations
using HydroModels

include("../models/exphydro.jl")

# define parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=100.0, soilwater=1303.004248)
pas = ComponentVector(params=params)

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))')
config = (timeidx=ts, interp=LinearInterpolation, solver=HydroModels.ManualSolver(mutable=true))
# run model with single node input
single_result = exphydro_model(input_arr, pas, initstates=init_states, config=config)
exphydro_model.infos
# # run model with multi node input
# node_num = 10
# node_params = ComponentVector(
#     f=fill(f, node_num), Smax=fill(Smax, node_num), Qmax=fill(Qmax, node_num),
#     Df=fill(Df, node_num), Tmax=fill(Tmax, node_num), Tmin=fill(Tmin, node_num)
# )
# node_states = ComponentVector(
#     snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num)
# )

# node_pas = ComponentVector(params=node_params, initstates=node_states)
# node_input = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), fill(input_arr, node_num)), (1, 3, 2))
# config = (ptyidx=1:node_num, styidx=1:node_num, timeidx=ts)
# result = exphydro_model(node_input, node_pas, config=config)


# tmp_input = ones(3, node_num)
# tmp_state = ones(1, node_num)
# exphydro_model.components[2].ode_funcs[2](eachslice(tmp_input, dims=1), eachslice(tmp_state, dims=1), node_pas)