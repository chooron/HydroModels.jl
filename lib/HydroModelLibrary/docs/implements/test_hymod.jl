using CSV
using Random
using DataFrames
using ComponentArrays
using CairoMakie
using CairoMakie: Axis
# using HydroModels

include("../../src/HydroModels.jl")
file_path = "E:\\JlCode\\hymod\\sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);
prcp_vec = df[!, "precip"]
et_vec = df[!, "pet"]
q_vec = df[!, "q"]

# build model
unit_params = ComponentVector(alpha=0.56377, bexp=0.45505, kf=0.96141, cmax=85.8419, ks=0.349047)
init_states = ComponentVector(soilwater=0.0, fastwater1=0.0, fastwater2=0.0, fastwater3=0.0, slowwater=0.0) # 
ps = ComponentVector(params=unit_params, initstates=init_states)

model = HydroModels.HyMOD.Unit(name=:hymod, mtk=true);
HydroModels.get_state_names(model)
input = (prcp=prcp_vec, pet=et_vec)
output = model(input, ps, timeidx=collect(1:length(q_vec)))
output_df = DataFrame(output)
# # plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, 1:length(et_vec), output[:flow], color=:blue)
lines!(ax, 1:length(et_vec), df[!, "q"] /1e3, color=:red)
fig