# 导入模块
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
using Lux, LuxCore
using StableRNGs
using OrdinaryDiffEq
# using HydroModels

include("../../src/HydroModels.jl")
model = HydroModels.M50.Unit(name=:m50, mtk=false, step=false)
# unit_sys = model.units[1]
# unknowns(unit_sys.system)
# base param names
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# input data
file_path = "data/m50/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)

mean_snowwater, std_snowwater = mean(df[ts, "SnowWater"]), std(df[ts, "SnowWater"])
mean_soilwater, std_soilwater = mean(df[ts, "SoilWater"]), std(df[ts, "SoilWater"])
mean_temp, std_temp = mean(df[ts, "Temp"]), std(df[ts, "Temp"])
mean_prcp, std_prcp = mean(df[ts, "Prcp"]), std(df[ts, "Prcp"])

input = StructArray((lday=df[ts, "Lday"], temp=df[ts, "Temp"], prcp=df[ts, "Prcp"])) # , snowwater=df[ts, "SnowWater"], infiltration=df[ts, "Infiltration"]),

et_ann = Lux.Chain(Lux.Dense(3, 16, Lux.tanh), Lux.Dense(16, 1, Lux.leakyrelu))
q_ann = Lux.Chain(Lux.Dense(2, 16, Lux.tanh), Lux.Dense(16, 1, Lux.leakyrelu))
et_ann_p = LuxCore.initialparameters(StableRNG(42), et_ann)
q_ann_p = LuxCore.initialparameters(StableRNG(42), q_ann)

params = ComponentVector(
    f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin,
    mean_snowwater=mean_snowwater, std_snowwater=std_snowwater,
    mean_soilwater=mean_soilwater, std_soilwater=std_soilwater,
    mean_temp=mean_temp, std_temp=std_temp,
    mean_prcp=mean_prcp, std_prcp=std_prcp)
initstates = ComponentVector(snowwater=0.0, soilwater=1303.004248)
pas = ComponentVector(params=params, initstates=initstates, nn=(etnn=et_ann_p, qnn=q_ann_p)) # , weight=1.0

solver = HydroModels.ODESolver(alg=Rosenbrock23())
results = model(input, pas, timeidx=ts, solver=solver)

# q_ann = Lux.Chain(
#     Lux.Dense(2 => 16, Lux.tanh),
#     # Lux.Dense(16 => 16, Lux.leakyrelu),
#     Lux.Dense(16 => 1, Lux.leakyrelu)
# )
# func = (x, p) -> LuxCore.stateless_apply(q_ann, x, p)

