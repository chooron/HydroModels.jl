using NamedTupleTools
using CSV
using DataFrames
using Plots
using Lux
using LuxCore
using StableRNGs
using Optimization
using OptimizationOptimisers
include("../../src/HydroModels.jl")

#* num of basin grid is 62
#* set the model parameters
FC, Beta, PWP, L = 177.1, 2.35, 105.89, 4.87
k0, k1, k2, kp = 0.5, 0.03, 0.02, 0.05

wnn = Lux.Chain(Lux.Dense(2 => 16), Lux.Dense(16 => 1, Lux.sigmoid_fast), name=:wnn)
ksnn = Lux.Chain(Lux.Dense(3 => 16), Lux.Dense(16 => 5, Lux.sigmoid_fast), name=:ksnn)
qnn = Lux.Chain(Lux.Dense(5 => 32), Lux.Dense(32 => 1, Lux.leakyrelu), name=:qnn)
wnn_ps = collect(ComponentVector(LuxCore.initialparameters(StableRNG(42), wnn)))
ksnn_ps = collect(ComponentVector(LuxCore.initialparameters(StableRNG(42), ksnn)))
q_ps = collect(ComponentVector(LuxCore.initialparameters(StableRNG(42), qnn)))
unit = LumpedHydro.HBV_NN.Unit(name=:hbvnn)
LumpedHydro.get_param_names(unit.components)

params = ComponentVector(PWP=PWP, L=L, kp=kp, k2=k2)
initstates = ComponentVector(soilwater=100.0, s1=3, s2=10)
# pas = ComponentVector(params=params, initstates=initstates, nn=(wnn=wnn_ps, ksnn=ksnn_ps, qnn=q_ps))
input = (prcp=ones(100), pet=zeros(100))
inputs = repeat([input], 5)
pas = ComponentVector(
    params=NamedTuple{Tuple([Symbol(:node_, i) for i in 1:5])}(repeat([params], 5)),
    initstates=NamedTuple{Tuple([Symbol(:node_, i) for i in 1:5])}(repeat([initstates], 5)),
    nn=(wnn=wnn_ps, ksnn=ksnn_ps, qnn=q_ps)
)

@btime unit_result = unit(inputs, pas, timeidx=collect(1:100))