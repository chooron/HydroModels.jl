using CSV, DataFrames
using ComponentArrays
using Distributions
using ModelingToolkit
using DataInterpolations
using Plots
include("../../src/HydroModels.jl")

df = CSV.read("data/marrmot/3604000.csv", DataFrame)

input = (P=df[!, "prec"], Ep=df[!, "pet"], T=df[!, "temp"])
model = HydroModels.load_model(:classic) 
param_bounds = getbounds.(HydroModels.get_params(model))
random_param_values = map(param_bounds) do param_bound
    rand(Uniform(param_bound[1], param_bound[2]))
end
init_params = ComponentVector(params=NamedTuple{Tuple(HydroModels.get_param_names(model))}(random_param_values))
init_states = NamedTuple{Tuple(HydroModels.get_state_names(model))}(zeros(length(HydroModels.get_state_names(model)))) |> ComponentVector
@info "Input variables: $(HydroModels.get_input_names(model))"
input_arr = stack(input[HydroModels.get_input_names(model)], dims=1)
config= (;solver=HydroModels.ODESolver(), interp=LinearInterpolation)
result = model(input_arr, init_params, initstates=init_states, config=config)
model_output_names = vcat(HydroModels.get_state_names(model), HydroModels.get_output_names(model)) |> Tuple
output_df = DataFrame(NamedTuple{model_output_names}(eachslice(result, dims=1)))

plot(df[!, "flow"])
plot!(output_df[!, "Qt"])