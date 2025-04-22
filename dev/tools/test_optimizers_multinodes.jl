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
include("../models/exphydro.jl")

# define parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
init_states = ComponentVector(snowpack=100.0, soilwater=1303.004248)
pas = ComponentVector(params=params)

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
ts = collect(1:100)
flow_vec = df[ts, "flow(mm)"]
input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
input_arr = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))')

# run model with multi node input
node_num_1 = 10
node_num_2 = 10
node_params = ComponentVector(
    f=fill(f, node_num_1), Smax=fill(Smax, node_num_1), Qmax=fill(Qmax, node_num_1),
    Df=fill(Df, node_num_1), Tmax=fill(Tmax, node_num_1), Tmin=fill(Tmin, node_num_1)
)
node_states = ComponentVector(snowpack=fill(0.0, node_num_2), soilwater=fill(1303.004248, node_num_2))

node_pas = ComponentVector(params=node_params)
input_arr = reduce(hcat, collect(input[HydroModels.get_input_names(exphydro_model)]))
node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], node_num_2))
node_input = permutedims(node_input, (2, 3, 1))
config = (ptyidx=[1, 2, 3, 4, 5, 2, 3, 4, 1, 5], styidx=collect(1:10), timeidx=ts,
    # solver=HydroModelTools.ODESolver(sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP()))
    solver=HydroModelTools.DiscreteSolver()
)
run_kwargs = (config=config,)
result = exphydro_model(node_input, node_pas, config=config, initstates=node_states)

tmp_input = ones(3, node_num_2)
tmp_state = ones(1, node_num_2)
exphydro_model.components[2].ode_funcs[2](eachslice(tmp_input, dims=1), eachslice(tmp_state, dims=1), node_pas)

ntp_pre = HydroModelTools.NamedTuplePreprocessor(exphydro_model.infos)
ntp_post = HydroModelTools.NamedTuplePostprocessor(exphydro_model.infos)
select_output = HydroModelTools.SelectComponentOutlet(exphydro_model.infos, 1)
opt_func(i, p, c) = @pipe (i, p, c) |> ntp_pre(_[1], _[2], _[3]) |> exphydro_model(_[1], _[2]; _[3]...) |> select_output(_) |> ntp_post(_)
output = (flow=flow_vec,)
input_ntp = fill(input, node_num_2)
opt_func(input_ntp, node_pas, run_kwargs)

test_input = rand(3, node_num_2)
test_state = rand(1, node_num_2)

exphydro_model(node_input, node_pas; initstates=node_states, config=config)

# Zygote.gradient(node_pas) do p
#     exphydro_model(node_input, p; initstates=node_states, config=config)[end, end, :] |> sum
# end

# Zygote.gradient(p -> opt_func(input_ntp, p, run_kwargs)[:flow] |> sum, node_pas)
# solve_alg = Adam()
# adtype = AutoZygote()
# grad_opt = HydroModelTools.HydroOptimizer(component=opt_func, maxiters=100, solve_alg=solve_alg, adtype=adtype)
# opt_params, loss_df = grad_opt([input_ntp], [output], run_kwargs=[run_kwargs], tunable_pas=tunable_pas, const_pas=const_pas, return_loss_df=true)