#!/usr/bin/env julia

"""
Test runner script for HydroModels.jl

Usage:
    julia test/run_tests.jl [test_group]

Test groups:
    all         - Run all tests (default)
    basic       - Run basic component tests
    single      - Run single-node model tests
    multi       - Run multi-node model tests
    routing     - Run routing component tests
    gradient    - Run gradient computation tests
    misc        - Run miscellaneous tests

Examples:
    julia test/run_tests.jl
    julia test/run_tests.jl basic
    julia test/run_tests.jl multi
"""

using Pkg
Pkg.activate(".")

using Test
using HydroModels
using CSV
using DataFrames
using Lux
using StableRNGs
using Statistics
using ComponentArrays
using Graphs

# Smooth step function for tests
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Include test helpers
include("test_helpers.jl")

# Parse command line arguments
test_group = length(ARGS) > 0 ? ARGS[1] : "all"

println("=" ^ 70)
println("Running HydroModels.jl tests: $test_group")
println("=" ^ 70)

if test_group == "all" || test_group == "basic"
    @testset "Basic Components" begin
        include("base/run_hydro_flux.jl")
        include("base/run_neural_flux.jl")
    end
end

if test_group == "all" || test_group == "single"
    @testset "Single Node Models" begin
        include("base/run_single_bucket.jl")
        include("base/run_single_lumped.jl")
    end
end

if test_group == "all" || test_group == "routing"
    @testset "Routing Components" begin
        include("base/run_unithydro.jl")
        include("base/run_hydro_route.jl")
    end
end

if test_group == "all" || test_group == "multi"
    @testset "Multi Node Models" begin
        include("base/run_multi_bucket.jl")
        include("base/run_multi_lumped.jl")
        include("base/run_spatial_model.jl")
    end
end

# if test_group == "gradient"
#     @testset "Gradient Computation" begin
#         include("gradient/run_hybrid_gradient.jl")
#         # include("gradient/run_hydro_gradient.jl")
#     end
# end

# if test_group == "misc"
#     @testset "Miscellaneous Tests" begin
#         include("miscellaneous/run_sort.jl")
#         include("miscellaneous/run_d8_grid.jl")
#         # include("miscellaneous/run_cuda.jl")
#     end
# end

println("\n" * "=" ^ 70)
println("Test suite completed!")
println("=" ^ 70)
