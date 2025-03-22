using Test
using HydroModels
using ModelingToolkit
using ComponentArrays

@testset "StateFlux Build Macro Tests" begin
    # Define variables and parameters
    @variables S, P, ET, Q
    @parameters a, b, c
    
    # Test 1: State flux with a name
    @testset "State flux with a name" begin
        flux = @stateflux_build :test_flux S = P - ET - Q
        
        # Check type and name
        @test flux isa StateFlux{:test_flux}
        
        # Check inputs, state, and parameters
        @test length(flux.meta.inputs) == 3
        @test length(flux.meta.states) == 1
        @test length(flux.meta.params) == 0
        
        # Check that inputs contain P, ET, and Q
        @test P in flux.meta.inputs
        @test ET in flux.meta.inputs
        @test Q in flux.meta.inputs
        
        # Check that state is S
        @test S in flux.meta.states
        
        # Check expression
        @test length(flux.exprs) == 1
        # The expression should be P - ET - Q
        @test string(flux.exprs[1]) == "P - ET - Q"
    end
    
    # Test 2: State flux with parameters
    @testset "State flux with parameters" begin
        flux = @stateflux_build :param_flux S = P - a*S - b*S^c
        
        # Check type and name
        @test flux isa StateFlux{:param_flux}
        
        # Check inputs, state, and parameters
        @test length(flux.meta.inputs) == 1  # P
        @test length(flux.meta.states) == 1  # S
        @test length(flux.meta.params) == 3  # a, b, c
        
        # Check that inputs contain P
        @test P in flux.meta.inputs
        
        # Check that state is S
        @test S in flux.meta.states
        
        # Check that parameters contain a, b, and c
        @test a in flux.meta.params
        @test b in flux.meta.params
        @test c in flux.meta.params
        
        # Check expression
        @test length(flux.exprs) == 1
        # The expression should be P - a*S - b*S^c
        @test string(flux.exprs[1]) == "P - a*S - b*S^c"
    end
    
    # Test 3: State flux without a name
    @testset "State flux without a name" begin
        flux = @stateflux_build S = P - ET - Q
        
        # Check that it's a StateFlux
        @test flux isa StateFlux
        
        # Check inputs, state, and parameters
        @test length(flux.meta.inputs) == 3
        @test length(flux.meta.states) == 1
        @test length(flux.meta.params) == 0
        
        # Check expression
        @test length(flux.exprs) == 1
        # The expression should be P - ET - Q
        @test string(flux.exprs[1]) == "P - ET - Q"
    end
    
    # Test 4: Error handling
    @testset "Error handling" begin
        # Test that an error is raised when the syntax is incorrect
        @test_throws ErrorException eval(:(
            @stateflux_build S + P  # Not an equation
        ))
    end
    
    # Test 5: State variable appearing on both sides
    @testset "State variable on both sides" begin
        flux = @stateflux_build S = P - a*S
        
        # Check inputs, state, and parameters
        @test length(flux.meta.inputs) == 1  # P
        @test length(flux.meta.states) == 1  # S
        @test length(flux.meta.params) == 1  # a
        
        # Check that inputs contain P but not S
        @test P in flux.meta.inputs
        @test !(S in flux.meta.inputs)  # S should not be in inputs
        
        # Check that state is S
        @test S in flux.meta.states
        
        # Check expression
        @test length(flux.exprs) == 1
        # The expression should be P - a*S
        @test string(flux.exprs[1]) == "P - a*S"
    end
end
