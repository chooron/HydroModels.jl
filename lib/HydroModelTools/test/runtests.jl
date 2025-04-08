using Aqua
using CSV
using DataFrames
using Lux
using Test
using StableRNGs
using ComponentArrays
using DataInterpolations
using OrdinaryDiffEq
using Statistics
using HydroModels
using HydroModelTools
using SciMLSensitivity
using OptimizationOptimisers

@testset "HydroModelTools.jl" begin
    include("run_exphydro_optimize.jl")
    include("run_m50_optimize.jl")
    include("run_solver.jl")
end