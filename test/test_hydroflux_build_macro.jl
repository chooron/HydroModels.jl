using Test
using HydroModels
using ModelingToolkit
using ComponentArrays

# Define variables and parameters
@variables x y z w g
@parameters a b c k

@testset "HydroFlux Build Macro Tests" begin
    # Test 1: Single equation with name
    @testset "Single equation with name" begin
        flux = @hydroflux :test_flux c ~ a * x + b * y

        # Check type and name
        @test flux isa HydroFlux{:test_flux}

        # Check inputs, outputs, and parameters
        @test length(flux.meta.inputs) == 2
        @test length(flux.meta.outputs) == 1
        @test length(flux.meta.params) == 2

        # Check that inputs contain x and y
        @test [:y, :x] == HydroModels.get_input_names(flux)

        # Check that outputs contain c
        @test [:c] == HydroModels.get_output_names(flux)

        # Check that parameters contain a and b
        @test [:a, :b] == HydroModels.get_param_names(flux)

        # Test the function with values
        params_cv = ComponentVector(params=(a=3.0, b=2.0))
        result = flux([1.0, 2.0], params_cv)
        @test result ≈ [3.0*2.0 + 2.0*1.0]
    end

    # Test 2: Multiple equations with name
    @testset "Multiple equations with name" begin
        flux = @hydroflux :multi_flux begin
            g ~ a*x + b*y
            w ~ c*x + k*y
        end

        # Check type and name
        @test flux isa HydroFlux{:multi_flux}

        # Check inputs, outputs, and parameters
        @test length(flux.meta.inputs) == 2
        @test length(flux.meta.outputs) == 2
        @test length(flux.meta.params) == 4

        # Check that inputs contain x and y
        @test [:y, :x] == HydroModels.get_input_names(flux)

        # Check that outputs contain c and w
        @test [:g, :w] == HydroModels.get_output_names(flux)

        # Check that parameters contain a, b, c, and k
        @test Set([:a, :b, :c, :k]) == Set(HydroModels.get_param_names(flux))

        # Test the function with values
        params_cv = ComponentVector(params=(a=2.0, b=3.0, c=1.5, k=0.5)[HydroModels.get_param_names(flux)])
        result = flux([2.0, 1.0], params_cv)
        @test result ≈ [2.0*1.0 + 3.0*2.0, 1.5*1.0 + 0.5*2.0]
    end

    # Test 3: Single equation without name
    @testset "Single equation without name" begin
        flux = @hydroflux z ~ a*x + b*y

        # Check that it's a HydroFlux
        @test flux isa HydroFlux

        # Check inputs, outputs, and parameters
        @test length(flux.meta.inputs) == 2
        @test length(flux.meta.outputs) == 1
        @test length(flux.meta.params) == 2

        # Test the function with values
        params_cv = ComponentVector(params=(a=2.0, b=3.0)[HydroModels.get_param_names(flux)])
        result = flux([2.0, 1.0], params_cv)
        @test result ≈ [2.0*1.0 + 3.0*2.0]
    end
end
