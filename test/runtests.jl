using Aqua
using CSV
using CUDA, cuDNN
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
    include("base/run_neural.jl")
    include("base/run_bucket.jl")
    include("base/run_uh.jl")
    include("base/run_route.jl")
    include("base/run_lumped_model.jl")
    include("base/run_spatial_model.jl")
end

# @testset "test cuda support" begin
#     include("miscellaneous/run_cuda.jl")
# end

# # cost a lot of time
# @testset "test gradient" begin
#     include("gradient/run_hydro_gradient.jl")
#     include("gradient/run_hybrid_gradient.jl")
# end