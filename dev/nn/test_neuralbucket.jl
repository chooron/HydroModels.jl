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
NeuralBucket = HydroModels.NeuralBucket
HydroModel = HydroModels.HydroModel

@variables lday temp prcp snowfall rainfall melt snowpack snowpack2
@parameters f Smax Qmax Df Tmax Tmin

#* build a neural bucket
@neuralbucket :bucket1 begin
    fluxes = begin
        HydroModels.@hydroflux begin
            snowfall ~ step_func(Tmin - temp) * prcp
            rainfall ~ step_func(temp - Tmin) * prcp
        end
        HydroModels.@hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
    end
    dfluxes = begin
        HydroModels.@stateflux snowpack ~ snowfall - melt
        HydroModels.@stateflux snowpack2 ~ snowfall - melt
    end
end

bucket_1 = NeuralBucket(fluxes=fluxes, dfluxes=dfluxes)
f_, Smax_, Qmax_, Df_, Tmax_, Tmin_ = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
node_num = 10
node_params = ComponentVector(
    f=fill(f_, node_num), Smax=fill(Smax_, node_num), Qmax=fill(Qmax_, node_num),
    Df=fill(Df_, node_num), Tmax=fill(Tmax_, node_num), Tmin=fill(Tmin_, node_num)
)
init_states = ComponentVector(snowpack=fill(100.0, node_num), snowpack2=fill(100.0, node_num))
pas = ComponentVector(params=node_params)
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:10000)

# single node input
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_mat = reduce(hcat, input[bucket_input_names]) |> permutedims
input_arr = permutedims(repeat(input_mat, 1, 1, node_num), (1, 3, 2))
result = bucket_1(input_arr, pas; initstates=init_states, solver=HydroModels.ManualSolver(mutable=false))
# Zygote.gradient(pas) do p
#     bucket_1(input_arr, p; initstates=init_states) |> sum
# end

