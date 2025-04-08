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

# load data
file_path = "data/cemaneige/sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:length(df[!, :precipitation]))
model = HydroModels.Cemaneige.SurfaceStorage(name=:cemaneige, mtk=false)
HydroModels.get_input_names(model)
HydroModels.get_state_names(model)
HydroModels.get_param_names(model)
input = (
    prcp=df[!, :precipitation], min_temp=df[!, :min_temp],
    mean_temp=df[!, :mean_temp], max_temp=df[!, :max_temp],
)
params = ComponentVector(CTG=0.25, Kf=3.74, altitude=550.0, zthresh=1500.0, snwthresh=129.38, height=495.0)
initstates = ComponentVector(snowwater=0.0, thermal=0.0)
pas = ComponentVector(params=params, initstates=initstates)
timeidx = collect(1:length(input[:prcp]))
result = model(input, pas, timeidx=timeidx, solved=false)
result_df = DataFrame(result)
CSV.write("result.csv", result_df)

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, ts[1:200], df[1:200, :liquid_outflow], color=:red)
lines!(ax, ts[1:200], result.infiltration[1:200], color=:blue)

fig