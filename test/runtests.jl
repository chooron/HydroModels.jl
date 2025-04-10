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

# @testset "HydroModels.jl" begin
#     include("run_flux.jl")
#     include("run_bucket.jl")
#     include("run_uh.jl")
#     include("run_route.jl")
#     include("run_lumped_model.jl")
#     include("run_spatial_model.jl")
# end

@testset "HydroModelSolvers.jl" begin
    include("solvers/run_solvers.jl")
end