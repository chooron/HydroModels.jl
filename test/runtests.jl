using CSV
using DataFrames
using Test
using StableRNGs
using Statistics
using ComponentArrays
using Graphs
using HydroModels

# Smooth step function for tests
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Include test helpers
include("test_helpers.jl")

@testset "HydroModels.jl" begin
    @testset "Core Package Loading" begin
        @test isnothing(Base.get_extension(HydroModels, :HydroModelsLuxExt))
        @test_throws ArgumentError create_neural_bucket(
            name = :stub_bucket,
            n_inputs = 1,
            n_states = 1,
            n_outputs = 1,
            inputs = [:prcp],
            states = [:soilwater],
            outputs = [:flow],
        )
    end

    @testset "Basic Components" begin
        include("base/run_hydro_flux.jl")
    end
    
    @testset "Single Node Models" begin
        include("base/run_single_bucket.jl")
        include("base/run_single_lumped.jl")
    end
    
    @testset "Routing Components" begin
        include("base/run_unithydro.jl")
        include("base/run_hydro_route.jl")
        include("base/run_channel_route.jl")
    end
    
    @testset "Multi Node Models" begin
        include("base/run_multi_bucket.jl")
        include("base/run_spatial_model.jl")
        include("base/run_multi_lumped.jl")
    end

    @testset "YAML Model Loading" begin
        include("base/run_yaml_model.jl")
    end

    Base.eval(Main, :(using Lux, LuxCore))

    @testset "Lux Extension Loading" begin
        @test !isnothing(Base.get_extension(HydroModels, :HydroModelsLuxExt))
    end

    @testset "NN Components" begin
        include("nn/run_neural_flux.jl")
    end

    @testset "NN Models" begin
        include("nn/run_single_lumped_nn.jl")
        include("nn/run_multi_lumped_nn.jl")
    end
end

# @testset "test cuda support" begin
#     include("miscellaneous/run_cuda.jl")
# end

@testset "AD and ODE gradients" begin
    include("gradient/test_ad.jl")
end
