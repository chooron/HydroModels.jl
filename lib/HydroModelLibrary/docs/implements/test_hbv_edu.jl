# import lib
using CSV
using DataFrames
using CairoMakie
using BenchmarkTools
using ComponentArrays

include("../../src/HydroModels.jl")
# load data
file_path = "data/hbv_edu/hbv_sample.csv"
data = CSV.File(file_path);
df = DataFrame(data);
pet_vec = df[!, "pet"]
prcp_vec = df[!, "prec"]
temp_vec = df[!, "temp"]
flow_vec = df[!, "qsim"]
ts = collect(1:length(pet_vec))

# build model
model = HydroModels.HBV_EDU.Unit(name=:hbv_edu, mtk=true)
HydroModels.get_state_names(model.elements[3])

# model parameters
tt, dd, FC, Beta, PWP, L = 0.0, 4.25, 177.1, 2.35, 105.89, 4.87
k0, k1, k2, kp = 0.05, 0.03, 0.02, 0.05
params = ComponentVector(tt=tt, dd=dd, FC=FC, Beta=Beta, PWP=PWP, L=L, k0=k0, k1=k1, k2=k2, kp=kp)
initstates = ComponentVector(snowwater=0.0, soilwater=100.0, s1=3, s2=10)
pas = ComponentVector(params=params, initstates=initstates)

input = (prcp=prcp_vec, pet=pet_vec, temp=temp_vec)
solver = HydroModels.ODESolver(reltol=1e-3, abstol=1e-3)
# solver = HydroModels.DiscreteSolver()
result = model(input, pas, timeidx=ts, solver=solver);

# plot result
fig = Figure(size=(400, 300))
ax = CairoMakie.Axis(fig[1, 1], title="predict results", xlabel="time", ylabel="flow(mm)")
lines!(ax, ts, df[!, "qsim"], color=:red)
lines!(ax, ts, result.flow, color=:blue)
fig