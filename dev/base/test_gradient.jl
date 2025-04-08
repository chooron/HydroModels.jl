using CSV
using DataFrames
using Lux
using StableRNGs
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Statistics
using BenchmarkTools
using Plots
using JLD2
using SciMLSensitivity
using HydroModels
using HydroModelTools
using Zygote

# include("E:\\JlCode\\HydroModels\\src\\HydroModels.jl")
HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
NeuralFlux = HydroModels.NeuralFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel
include("../models/m50.jl")

# load data
data = CSV.read("data/exphydro/01013500.csv", DataFrame)
ts = collect(1:100)
input = (lday=data[ts, "dayl(day)"], temp=data[ts, "tmean(C)"], prcp=data[ts, "prcp(mm/day)"])
flow_vec = data[ts, "flow(mm)"]
config = (solver=ODESolver(sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP())), interp=LinearInterpolation, alg=BS3())

m50_opt_params = load("E:\\JlCode\\HydroModelsPaperCode\\src\\implements\\save\\m50_opt.jld2")["opt_params"]
input_mat = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(m50_model)]))')
m50_output = m50_model(input_mat, m50_opt_params, config=config)

function loss1(p)
    m50_output = m50_ele(input_mat, p, config=config)
    sum(m50_output[end,:] .- flow_vec)
end

gradient(loss1, m50_opt_params)
