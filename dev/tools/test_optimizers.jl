using CSV, DataFrames, ComponentArrays
using HydroModels
using DataInterpolations
using OptimizationOptimisers
using SciMLSensitivity
using BenchmarkTools
using Zygote
using Pipe
include("../models/exphydro.jl")
include("../src/HydroModelTools.jl")

# define parameters and initial states
f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

# load data
file_path = "data/exphydro/01013500.csv"
data = CSV.File(file_path)
df = DataFrame(data)
ts = collect(1:10000)
lday_vec = df[ts, "dayl(day)"]
prcp_vec = df[ts, "prcp(mm/day)"]
temp_vec = df[ts, "tmean(C)"]
flow_vec = df[ts, "flow(mm)"]

tunable_pas = ComponentVector(params=ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
const_pas = ComponentVector()

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
input_matrix = Matrix(reduce(hcat, collect(input))')
output = (flow=flow_vec,)
config = (solver = HydroModelTools.ODESolver(), interp = LinearInterpolation)
initstates = ComponentVector(snowpack=0.0, soilwater=1300.0)
run_kwargs = (config = config, initstates = initstates)
ntp_pre, ntp_post = HydroModelTools.NamedTuplePreprocessor(exphydro_model), HydroModelTools.NamedTuplePostprocessor(exphydro_model)
opt_func(i, p, c) = @pipe (i, p, c) |> ntp_pre(_[1], _[2], _[3]) |> exphydro_model(_[1], _[2]; _[3]...) |> ntp_post(_)

opt_func(input, tunable_pas, run_kwargs)
# build optimizer
hydro_opt = HydroModelTools.HydroOptimizer(component=opt_func, maxiters=100)
lb_list = [0.0, 100.0, 10.0, 0.0, 0.0, -3.0]
ub_list = [0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]
# opt_params, loss_df = hydro_opt([input], [output], run_kwargs=[run_kwargs], tunable_pas=tunable_pas, const_pas=const_pas, lb=lb_list, ub=ub_list, return_loss_df=true)

Zygote.gradient(p -> opt_func(input, p, run_kwargs)[:flow] |> sum, tunable_pas)
# grad_opt = HydroModelTools.GradOptimizer(component=opt_func, maxiters=100, solve_alg=Adam(1e-3), adtype=AutoForwardDiff())
# opt_params, loss_df = grad_opt([input], [output], run_kwargs=[run_kwargs], tunable_pas=tunable_pas, const_pas=const_pas, return_loss_df=true)