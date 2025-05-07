# 导入模块
using CSV
using DataFrames
using ComponentArrays
using DataInterpolations
using ModelingToolkit
using ModelingToolkit: isparameter, get_variables
using HydroModelTools
using Zygote
using SciMLSensitivity
using BenchmarkTools
using OrdinaryDiffEq
using Lux

include("../../src/HydroModels.jl")
NeuralFlux = HydroModels.NeuralFlux
@variables a b c d

ep_nn = Lux.Chain(
    Lux.Dense(3 => 16, tanh),
    Lux.Dense(16 => 16, leakyrelu),
    Lux.Dense(16 => 1, leakyrelu),
    name=:epnn
)
nnflux_1 = HydroModels.@neuralflux d ~ ep_nn([a, b, c])
