include("../../src/HydroModels.jl")
using ModelingToolkit
using ModelingToolkit: isparameter, get_variables
using ComponentArrays

@parameters lag
@variables q t

UnitHydrograph = HydroModels.UnitHydrograph

uh1 =  HydroModels.@unithydro :maxbas_uh begin
    uh_func = begin
        2lag => (1 - 0.5 * (2 - t / lag)^2.5)
        lag => (0.5 * (t / lag)^2.5)
    end
    uh_vars = [q]
    configs = (solvetype=:SPARSE, suffix=:_lag)
end

input_flow = Float32[2 3 4 2 3 1]
params = ComponentVector(params=(lag=3.5,))

uh1(input_flow, params)