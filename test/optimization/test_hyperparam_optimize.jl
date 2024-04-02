# import lib
using CSV
using DataFrames
using CairoMakie
using ComponentArrays
using OptimizationOptimisers
using BenchmarkTools

# test exphydro model
include("../../src/DeepFlex.jl")

# predefine the parameters
# init_parameter = [0.0, 100.0, 0.01, 20, 1.0, 1.0, -1.0]
init_parameter = [0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084]
search_params = [
    (name=:f, default=init_parameter[1], lb=0.0, ub=0.1),
    (name=:Smax, default=init_parameter[2], lb=100.0, ub=1500.0),
    (name=:Qmax, default=init_parameter[3], lb=10.0, ub=50.0),
    (name=:Df, default=init_parameter[4], lb=0.0, ub=5.0),
    (name=:Tmax, default=init_parameter[5], lb=0.0, ub=3.0),
    (name=:Tmin, default=init_parameter[6], lb=-3.0, ub=0.0),
]
const_params = [(name=:snowwater, value=0.0), (name=:soilwater, value=1300.0)]
model = DeepFlex.ExpHydro(name=:exphydro)

# load data
file_path = "data/camels/01013500.csv"
data = CSV.File(file_path);
df = DataFrame(data);
lday_vec = df[1:1000, "dayl(day)"]
prcp_vec = df[1:1000, "prcp(mm/day)"]
temp_vec = df[1:1000, "tmean(C)"]
flow_vec = df[1:1000, "flow(mm)"]

# parameters optimization
input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec, time=1:1:length(lday_vec))
output = (flow=flow_vec,)

# best_params = DeepFlex.param_box_optim(
#     model,
#     search_params=search_params,
#     const_params=const_params,
#     input=input,
#     target=output,
# )

best_params = DeepFlex.param_grad_optim(
    model,
    search_params=search_params,
    const_params=const_params,
    input=input,
    target=output,
)

total_params = merge(best_params, const_params)
reulst = model(input, total_params[model.nameinfo.param_names], total_params[model.nameinfo.state_names])

