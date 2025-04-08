include("../../src/HydroModels.jl")
using ModelingToolkit
using Lux

HydroFlux = HydroModels.HydroFlux
StateFlux = HydroModels.StateFlux
NeuralFlux = HydroModels.NeuralFlux
HydroBucket = HydroModels.HydroBucket
HydroModel = HydroModels.HydroModel

@variables a b c d e
@parameters k1 k2
eq1 = a ~ k1 * b + k2 * c
flux1 = HydroModels.@hydroflux :flux1 a ~ k1 * b -k2 * c
flux1 = HydroModels.@hydroflux :flux1 begin
    a ~ k1 * b + k2 * c
    d ~ k1 * b + k2 * c
end
dflux1 = HydroModels.@stateflux :dflux1 a ~ b + c
dflux1 = HydroModels.@stateflux :dflux1 a ~ c - b

chain = Lux.Chain(Lux.Dense(2, 64), Lux.relu, Lux.Dense(64, 1), name=:m)
HydroModels.@neuralflux :neuralflux [d, e] ~ chain([a,b])
HydroModels.@neuralflux :neuralflux d ~ chain([a,b])

bucket1 = HydroModels.@hydrobucket :bucket1 begin
    fluxes = [
        HydroModels.@hydroflux :flux1 a ~ k1 * b - k2 * c
        HydroModels.@hydroflux :flux2 d ~ k1 * b - k2 * c
    ]
    dfluxes = [
        HydroModels.@stateflux :dflux1 c ~ b - a - d
    ]
end

bucket2 = HydroModels.@hydrobucket :route1 begin
    fluxes = [
        HydroModels.@hydroflux :flux1 a ~ k1 * b - k2 * c
        HydroModels.@hydroflux :flux2 d ~ k1 * b - k2 * c
    ]
    dfluxes = [
        HydroModels.@stateflux :dflux1 c ~ b - a - d
    ]

    rfunc = ...
end

HydroModels.@hydromodel :model1 begin
    bucket1
    HydroModels.@neuralflux :neuralflux1 [e] ~ chain([a,b])
end
