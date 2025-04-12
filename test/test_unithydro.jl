using HydroModels
using Test

@testset "UnitHydrograph Macro" begin
    # Test the @unithydro macro
    @variables q
    @parameters lag
    
    # Create a unit hydrograph using the macro
    uh = @unithydro :maxbas_uh begin
        uh_func = begin
            2lag => (1 - 0.5 * (2 - t / lag)^2.5)
            lag => (0.5 * (t / lag)^2.5)
        end
        uh_vars = [q]
        configs = (solvetype=:DISCRETE, suffix=:_lag)
    end
    
    # Test that the unit hydrograph was created correctly
    @test uh isa UnitHydrograph
    @test get_name(uh) == :maxbas_uh
    @test get_input_names(uh) == [:q]
    @test get_output_names(uh) == [:q_lag]
    
    # Test with a single segment
    uh_single = @unithydro :simple_uh begin
        uh_func = begin
            lag => (t / lag)^2
        end
        uh_vars = [q]
        configs = (solvetype=:SPARSE, suffix=:_routed)
    end
    
    @test uh_single isa UnitHydrograph
    @test get_name(uh_single) == :simple_uh
    @test get_input_names(uh_single) == [:q]
    @test get_output_names(uh_single) == [:q_routed]
    
    # Test with multiple variables
    @variables q1 q2
    uh_multi = @unithydro :multi_uh begin
        uh_func = begin
            2lag => (1 - 0.5 * (2 - t / lag)^2.5)
            lag => (0.5 * (t / lag)^2.5)
        end
        uh_vars = [q1, q2]
        configs = (solvetype=:DISCRETE, suffix=:_lag)
    end
    
    @test uh_multi isa UnitHydrograph
    @test get_name(uh_multi) == :multi_uh
    @test get_input_names(uh_multi) == [:q1, :q2]
    @test get_output_names(uh_multi) == [:q1_lag, :q2_lag]
end
