using Aqua
using CSV
using DataFrames
using Lux
using Test
using StableRNGs
using Statistics
using ComponentArrays
using DataInterpolations
using Graphs
using HydroModels
using Test
using DifferentialEquations
using SciMLSensitivity
using Zygote, ForwardDiff

@testset "HydroModels.jl" begin
    include("base/run_hydro_flux.jl")
    include("base/run_neural_flux.jl")
    include("base/run_single_bucket.jl")
    include("base/run_single_lumped.jl")
    include("base/run_unithydro.jl")
    include("base/run_hydro_route.jl")
    include("base/run_multi_bucket.jl")
    include("base/run_multi_lumped.jl")
    include("base/run_spatial_model.jl")
#     # include("base/run_neural.jl")
#     include("base/run_uh.jl")
#     include("base/run_route.jl")
#     include("base/run_lumped_model_single.jl")
#     include("base/run_spatial_model.jl")
end

# @testset "test cuda support" begin
#     include("miscellaneous/run_cuda.jl")
# end

# # cost a lot of time
# @testset "test gradient" begin
#     # include("gradient/run_hydro_gradient.jl")
#     include("gradient/run_hybrid_gradient.jl")
# end