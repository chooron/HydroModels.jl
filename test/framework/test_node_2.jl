# 导入模块
using ModelingToolkit
using CSV
using DataFrames
using ComponentArrays
using BenchmarkTools
using NamedTupleTools
using DataFrames
include("../../src/DeepFlex.jl")

model = DeepFlex.HydroNode(
    :exphydro_node,
    units=[DeepFlex.ExpHydro.Unit(name=:exphydro1), DeepFlex.ExpHydro.Unit(name=:exphydro2)],
    routes=namedtuple([:exphydro1, :exphydro2],
        [DeepFlex.ExpHydro.Route(name=:exphydro1), DeepFlex.ExpHydro.Route(name=:exphydro2)])
)
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = 1:100

unit_params = (f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
unit_input = (time=ts, lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
unit_init_states = (snowwater=0.0, soilwater=1303.004248)

input = (exphydro1=unit_input, exphydro2=unit_input)
ps = ComponentVector(
    exphydro1=(unit=(params=unit_params, initstates=unit_init_states), route=(params=ComponentVector(),), weight=0.5),
    exphydro2=(unit=(params=unit_params, initstates=unit_init_states), route=(params=ComponentVector(),), weight=0.5)
)

results = model(input, ps, step=false) # 分开计算和同时计算消耗的资源差不多