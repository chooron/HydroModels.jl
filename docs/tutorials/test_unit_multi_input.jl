# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
include("../../src/LumpedHydro.jl")

unit = LumpedHydro.ExpHydro.Unit(name=:exphydro)

node_num = 10
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params_i = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states_i = ComponentVector(snowwater=0.0, soilwater=1303.004248)
params = ComponentVector(NamedTuple{Tuple([Symbol(:node, i) for i in 1:node_num])}(repeat([params_i], node_num)))
init_states = ComponentVector(NamedTuple{Tuple([Symbol(:node, i) for i in 1:node_num])}(repeat([init_states_i], node_num)))
pas = ComponentVector(params=params, initstates=init_states)

file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
inputs = repeat([input], node_num)

solver = LumpedHydro.ODESolver()
results = unit(inputs, pas, timeidx=ts, solver=solver)
# input_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), [reduce(hcat, [input[nm] for nm in unit.input_names]) for input in inputs])
# results[1]
# plot(results[1].flow)