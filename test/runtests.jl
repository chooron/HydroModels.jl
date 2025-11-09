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
using DifferentialEquations
using SciMLSensitivity

# Smooth step function for tests
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Include test helpers
include("test_helpers.jl")

@testset "HydroModels.jl" begin
    @testset "Basic Components" begin
        include("base/run_hydro_flux.jl")
        include("base/run_neural_flux.jl")
    end
    
    @testset "Single Node Models" begin
        include("base/run_single_bucket.jl")
        include("base/run_single_lumped.jl")
    end
    
    @testset "Routing Components" begin
        include("base/run_unithydro.jl")
        include("base/run_hydro_route.jl")
    end
    
    @testset "Multi Node Models" begin
        include("base/run_multi_bucket.jl")
        include("base/run_multi_lumped.jl")
        include("base/run_spatial_model.jl")
    end
end

# @testset "test cuda support" begin
#     include("miscellaneous/run_cuda.jl")
# end

# # cost a lot of time
# @testset "test gradient" begin
#     # include("gradient/run_hydro_gradient.jl")
#     include("gradient/run_hybrid_gradient.jl")
# end