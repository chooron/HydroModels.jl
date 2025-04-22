# 测试HydroRec
# 导入模块
using CSV
using DataFrames
using ComponentArrays
# using HydroModels
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: isparameter, get_variables
using HydroModelTools
using Zygote
using SciMLSensitivity
using BenchmarkTools
using OrdinaryDiffEq

include("../../src/HydroModels.jl")
HydroFlux = HydroModels.HydroFlux
NeuralFlux = HydroModels.NeuralFlux
StateFlux = HydroModels.StateFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/exphydro.jl")

f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

bucket_input_names = HydroModels.get_input_names(bucket_1)
bucket_param_names = HydroModels.get_param_names(bucket_1)
bucket_state_names = HydroModels.get_state_names(bucket_1)
# 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)[bucket_param_names]
init_states = ComponentVector(snowpack=100.0)[bucket_state_names]
pas = ComponentVector(params=params)
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = reduce(hcat, input[bucket_input_names]) |> permutedims
result = bucket_1(input_arr, pas; solver=HydroModelTools.ODESolver(sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())), interp=LinearInterpolation)
tmp_func1(p, k) = bucket_1(input_arr, p; k...)