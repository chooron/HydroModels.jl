# import lib
using CSV
using Random
using DataFrames
using CairoMakie
using CairoMakie: Axis
using BenchmarkTools
using ModelingToolkit
using OrdinaryDiffEq
using ComponentArrays
using StructArrays
include("../../src/HydroModels.jl")
# using HydroModels

# test gr4j model
# load data
file_path = "data/symhyd/sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);
# time = 1:10000
prcp_vec = df[!, "precip"]
et_vec = df[!, "pet"]

model = HydroModels.SIMHYD.Unit(name=:simhyd, mtk=true);
HydroModels.get_input_names(model)
HydroModels.get_output_names(model)
HydroModels.get_param_names(model)
HydroModels.get_state_names(model)

# build model
input = (prcp=prcp_vec, pet=et_vec)
unit_params = (INSC=25.0, COEFF=200.0, SQ=5.0, SMSC=200.0, SUB=0.2, CRAK=0.2, K=0.5)
init_states = (SMS=100.0, GW=0.0)
pas = ComponentVector(params=unit_params, initstates=init_states)
result = model(input, pas, timeidx=collect(1:length(prcp_vec)));

# plot result
fig = Figure(size=(400, 400))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, 1:1:length(prcp_vec), result.U, color=:blue)
lines!(ax, 1:1:length(prcp_vec), df[!, "qsim"], color=:green)
fig