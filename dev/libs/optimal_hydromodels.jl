using CSV, DataFrames, Dates, ComponentArrays
using HydroModels, HydroModelTools
using DataInterpolations
using OptimizationEvolutionary
using OptimizationGCMAES
using OptimizationBBO
using Optimization
using ModelingToolkit
using Distributions
using ProgressMeter
using Plots
using JLD2
include("../src/HydroModelLibrary.jl")

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:10000)
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
pet_vec = @. 29.8 * lday_vec * 24 * 0.611 * exp((17.3 * temp_vec) / (temp_vec + 237.3)) / (temp_vec + 273.2)
flow_vec = df[ts, "flow(mm)"]
warm_up = 365
max_iter = 1000

model_nm = "hbv"
model = HydroModelLibrary.load_model(Symbol(model_nm))
param_bounds = getbounds.(get_params(model))
random_param_values = map(param_bounds) do param_bound
    rand(Uniform(param_bound[1], param_bound[2]))
end
init_params = ComponentVector(params=NamedTuple{Tuple(get_param_names(model))}(random_param_values))
ps_axes = getaxes(init_params)
init_states = NamedTuple{Tuple(get_state_names(model))}(zeros(length(get_state_names(model)))) |> ComponentVector

input = (P=prcp_vec, Ep=pet_vec, T=temp_vec)
input_matrix = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(model)]))')
config = (solver=HydroModelTools.DiscreteSolver(), interp=LinearInterpolation)
run_kwargs = (config=config, initstates=init_states)
y_mean = mean(flow_vec)
r2_func(y, y_hat) = sum((y .- y_hat) .^ 2) ./ sum((y .- y_mean) .^ 2)
kge_func(y, y_hat) = sqrt((r2_func(y, y_hat))^2 + (std(y_hat) / std(y) - 1)^2 + (mean(y_hat) / mean(y) - 1)^2)
function obj_func(p, _)
    return kge_func(
        flow_vec[warm_up:end],
        model(input_matrix, ComponentVector(p, ps_axes); run_kwargs...)[end, warm_up:end]
    )
end

progress = Progress(max_iter, desc="Optimization")
recorder = []
callback_func!(state, l) = begin
    push!(recorder, (iter=state.iter, loss=l, time=now(), params=state.u))
    next!(progress)
    false
end

lb_list = first.(param_bounds) .|> eltype(input_matrix)
ub_list = last.(param_bounds) .|> eltype(input_matrix)

optf = Optimization.OptimizationFunction(obj_func, Optimization.AutoForwardDiff())
optprob = Optimization.OptimizationProblem(optf, Vector(init_params), lb=lb_list, ub=ub_list)
sol = Optimization.solve(
    optprob,
    # CMAES(),
    BBO_adaptive_de_rand_1_bin_radiuslimited(),
    # GCMAESOpt(),
    maxiters=max_iter,
    callback=callback_func!
)
recorder_df = DataFrame(recorder)
CSV.write("cache/recorder_$(model_nm).csv", recorder_df)

params = ComponentVector(sol.u, ps_axes)
output = model(input_matrix, params; run_kwargs...)
model_output_names = vcat(get_state_names(model), get_output_names(model)) |> Tuple
output_df = DataFrame(NamedTuple{model_output_names}(eachslice(output, dims=1)))
@info loss = kge_func(flow_vec, output_df[!, :q_routed])
plot(output_df[!, :q_routed])
plot!(flow_vec)