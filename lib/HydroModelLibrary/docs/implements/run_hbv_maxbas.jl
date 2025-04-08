using CSV
using Plots
using Random
using DataFrames
using BenchmarkTools
using ComponentArrays
include("../../src/HydroModels.jl")

FC, Beta, PWP, L = 177.1, 2.35, 105.89, 4.87
k0, k1, k2, kp = 0.05, 0.03, 0.02, 0.05
resolution = 0.025
area_conv = (resolution * 110.94)^2.0 * 1000 / 3600
params_i = ComponentVector(FC=FC, Beta=Beta, PWP=PWP, L=L, k0=k0, k1=k1, k2=k2, kp=kp, conv=area_conv, lag=1.5)
initstates_i = ComponentVector(soilwater=100.0, s1=3, s2=10)

node_num = 10
params = ComponentVector(NamedTuple{Tuple([Symbol(:node, i) for i in 1:node_num])}(repeat([params_i], node_num)))
init_states = ComponentVector(NamedTuple{Tuple([Symbol(:node, i) for i in 1:node_num])}(repeat([initstates_i], node_num)))
pas = ComponentVector(params=params, initstates=init_states)
# pas = ComponentVector(param=params_i, state=initstates_i)
unit = HydroModels.HBV_MAXBAS.Unit(name=:hbv_maxbas)
input = (prcp=ones(100) .* 5, pet=zeros(100))
inputs = repeat([input], node_num)

output = unit(inputs, pas, timeidx=collect(1:100))