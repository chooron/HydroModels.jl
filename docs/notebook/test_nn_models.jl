# 导入模块
using CSV
using DataFrames
using ComponentArrays
using ModelingToolkit
using Lux
using StableRNGs
using Statistics

include("../src/HydroModels.jl")
HydroFlux = HydroModels.HydroFlux
NeuralFlux = HydroModels.NeuralFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/m50.jl")

#! load data
df = DataFrame(CSV.File("data/m50/01013500.csv"))
ts = collect(1:10000)
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
dayl_vec = df[ts, "Lday"]
snowpack_vec = df[ts, "SnowWater"]
soilwater_vec = df[ts, "SoilWater"]
qobs_vec = df[ts, "Flow"]

inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
means, stds = mean.(inputs), std.(inputs)
(prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec) = [@.((inp - mean) / std) for (inp, mean, std) in zip(inputs, means, stds)]

et_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), ep_nn))
q_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn))

base_params = (Df=2.674, Tmax=0.17, Tmin=-2.09)
var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
nn_params = (epnn=et_nn_p, qnn=q_nn_p)
params = reduce(merge, [base_params, var_means, var_stds])
initstates = (snowpack=0.0, soilwater=1303.00)
pas = ComponentVector(params=params, nns=nn_params)
input_ntp = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
input_mat = Matrix(reduce(hcat, collect(input_ntp[HydroModels.get_input_names(m50_model)]))')

# Run the model and get results as a matrix
result_mat = m50_model(input_mat, pas, config=(timeidx=ts,), initstates=ComponentVector(initstates))
# # run model with multi node input
# node_num = 10
# node_params = (Df=fill(2.674, node_num), Tmax=fill(0.17, node_num), Tmin=fill(-2.09, node_num))
# node_states = (snowpack=fill(0.0, node_num), soilwater=fill(1303.004248, node_num))
# node_var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(eachslice(repeat(reshape(stds, 1, :), 10), dims=2))
# node_var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(eachslice(repeat(reshape(means, 1, :), 10), dims=2))

# node_pas = ComponentVector(params=reduce(merge, [node_params, node_var_means, node_var_stds]), initstates=node_states, nns=nn_params)
# input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(m50_model)]))
# node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], 10))
# node_input = permutedims(node_input, (2, 3, 1))
# config = (ptyidx=1:10, styidx=1:10, timeidx=ts)
# result = m50_model(node_input, node_pas, config=config, initstates=ComponentVector(node_states))
# m50_model.components[2].flux_funcs[2]