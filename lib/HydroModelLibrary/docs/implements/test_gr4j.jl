# import lib
# tested in https://github.com/ckrapu/gr4j_theano
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
file_path = "data/gr4j/sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);
# time = 1:10000
prcp_vec = df[!, "prec"]
et_vec = df[!, "pet"]

# build model
unit_params = (x1=320.11, x2=2.42, x3=69.63, x4=1.39, ω=3.5, γ=5.0)
init_states = (soilwater=235.966719473926, routingstore=45.4697)
pas = ComponentVector((params=unit_params, initstates=init_states))

model = HydroModels.GR4J.Unit(name=:gr4j, mtk=false)
# model.elements[2].dfuncs[1].inner_func
input = (prcp=prcp_vec, pet=et_vec)
solver = HydroModels.DiscreteSolver()
result = model(input, pas, timeidx=collect(1:length(prcp_vec)), solver=solver);
result_df = DataFrame(result)

# plot result
fig = Figure(size=(400, 400))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
x = range(1, 100, length=100)
lines!(ax, 1:1:length(prcp_vec), result.flow, color=:blue)
lines!(ax, 1:1:length(prcp_vec), df[!, "qsim"], color=:green)
fig