using CSV, DataFrames, ComponentArrays
using HydroModels
using DataInterpolations
using OptimizationOptimisers
using SciMLSensitivity
using Statistics
using Zygote
using Lux
using StableRNGs
using Pipe
using BenchmarkTools: @btime
using HydroModelTools
include("../models/m50.jl")

# define parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

df = DataFrame(CSV.File("data/m50/01013500.csv"))
ts = collect(1:10000)
prcp_vec = df[ts, "Prcp"]
temp_vec = df[ts, "Temp"]
dayl_vec = df[ts, "Lday"]
snowpack_vec = df[ts, "SnowWater"]
soilwater_vec = df[ts, "SoilWater"]
flow_vec = df[ts, "Flow"]

inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
means, stds = mean.(inputs), std.(inputs)
(prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec) = [@.((inp - mean) / std) for (inp, mean, std) in zip(inputs, means, stds)]

et_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), ep_nn)) |> Vector
q_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn)) |> Vector

base_params = (Df=2.674, Tmax=0.17, Tmin=-2.09)
var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
nn_params = (epnn=et_nn_p, qnn=q_nn_p)
params = reduce(merge, [base_params, var_means, var_stds])
initstates = ComponentVector(snowpack=0.0, soilwater=1303.00)
tunable_pas = ComponentVector(params=params, nns=nn_params)
const_pas = ComponentVector()
input_ntp = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
input_mat = Matrix(reduce(hcat, collect(input_ntp[HydroModels.get_input_names(m50_model)]))')

# parameters optimization
output = (flow=flow_vec,)
config = (solver=HydroModelTools.ODESolver(sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())), interp=LinearInterpolation)
initstates = ComponentVector(snowpack=0.0, soilwater=1300.0)
run_kwargs = (config=config, initstates=initstates)
# ntp_pre, ntp_post = HydroModelTools.NamedTuplePreprocessor(m50_model), HydroModelTools.NamedTuplePostprocessor(m50_model)
# opt_func(i, p, c) = @pipe (i, p, c) |> ntp_pre(_[1], _[2], _[3]) |> m50_model(_[1], _[2]; _[3]...) |> ntp_post(_)
# opt_func(input_ntp, tunable_pas, run_kwargs)
# # opt_func(input_ntp, tunable_pas, run_kwargs)
Zygote.gradient(p -> m50_model(input_mat, p; run_kwargs...)[end,:] |> sum, tunable_pas)