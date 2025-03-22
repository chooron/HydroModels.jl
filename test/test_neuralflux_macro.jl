using Test
using HydroModels
using ModelingToolkit
using Lux
using Random
using ComponentArrays

@testset "NeuralFlux Macro Tests" begin
    # Define variables
    @variables x, y, z₁, z₂
    
    # Set random seed for reproducibility
    rng = StableRNG(42)
    
    # Test 1: Single output neural flux
    @testset "Single output neural flux" begin
        # Create a neural network chain
        chain = Chain(
            Dense(2 => 10, relu),
            Dense(10 => 1),
            name=:test_net1
        )
        
        # Create a neural flux
        flux = @neuralflux z₁ ~ chain([x, y])
        
        # Check type
        @test flux isa NeuralFlux
        
        # Check inputs, outputs, and neural network parameters
        @test length(flux.meta.inputs) == 2
        @test length(flux.meta.outputs) == 1
        @test :test_net1 in keys(flux.meta.nns)
        
        # Check that inputs contain x and y
        @test x in flux.meta.inputs
        @test y in flux.meta.inputs
        
        # Check that outputs contain z₁
        @test z₁ in flux.meta.outputs
        
        # Test the function with values
        params = ComponentVector(Lux.initialparameters(rng, chain))
        params_cv = ComponentVector(test_net1=params)
        result = flux([1.0, 2.0], params_cv)
        
        # Check that the result has the correct shape
        @test size(result) == (1,)
    end
    
    # Test 2: Multiple output neural flux
    @testset "Multiple output neural flux" begin
        # Create a neural network chain
        chain = Chain(
            Dense(2 => 16, relu),
            Dense(16 => 2),
            name=:test_net2
        )
        
        # Create a neural flux
        flux = @neuralflux [z₁, z₂] ~ chain([x, y])
        
        # Check type
        @test flux isa NeuralFlux
        
        # Check inputs, outputs, and neural network parameters
        @test length(flux.meta.inputs) == 2
        @test length(flux.meta.outputs) == 2
        @test :test_net2 in keys(flux.meta.nns)
        
        # Check that inputs contain x and y
        @test x in flux.meta.inputs
        @test y in flux.meta.inputs
        
        # Check that outputs contain z₁ and z₂
        @test z₁ in flux.meta.outputs
        @test z₂ in flux.meta.outputs
        
        # Test the function with values
        params = ComponentVector(Lux.initialparameters(rng, chain))
        params_cv = ComponentVector(test_net2=params)
        result = flux([1.0, 2.0], params_cv)
        
        # Check that the result has the correct shape
        @test size(result) == (2,)
    end
    
    # Test 3: Error handling
    @testset "Error handling" begin
        # Create a neural network chain
        chain = Chain(
            Dense(2 => 10, relu),
            Dense(10 => 1),
            name=:test_net3
        )
        
        # Test that an error is raised when the syntax is incorrect
        @test_throws ErrorException eval(:(
            @neuralflux z₁ = chain([x, y])  # Using = instead of ~
        ))
        
        # Test that an error is raised when inputs are not provided as an array
        @test_throws ErrorException eval(:(
            @neuralflux z₁ ~ chain(x, y)  # Not using an array for inputs
        ))
    end
end
