# Test routing components

@testset "Grid-based routing" begin
    # Define grid structure
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    # Define variables and parameters
    @variables q q_routed s_river
    @parameters lag

    # Create routing component
    grid_route = @hydroroute begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag)
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(flwdir, positions)
        hru_types = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    end

    @testset "Route interface" begin
        @test HydroModels.get_input_names(grid_route) == [:q]
        @test HydroModels.get_output_names(grid_route) == [:q_routed]
        @test HydroModels.get_param_names(grid_route) == [:lag]
        @test HydroModels.get_state_names(grid_route) == [:s_river]
    end

    @testset "Run with MutableSolver" begin
        input_arr = ones(1, 9, 20)
        timeidx = collect(1:20)
        params = ComponentVector(params = (lag = fill(0.2, 3),))
        initstates = ComponentVector(s_river = fill(0.0, length(positions)))
        
        config = create_test_config(solver = MutableSolver, timeidx = timeidx)
        output_arr = grid_route(input_arr, params, config; initstates = initstates)
        
        @test size(output_arr) == (2, 9, 20)  # 2 = state + output
        
        # Check that outputs are reasonable
        @test all(output_arr[1, :, :] .>= 0)  # River storage should be non-negative
        @test all(output_arr[2, :, :] .>= 0)  # Routed flow should be non-negative
    end
    
    @testset "Run with ImmutableSolver" begin
        input_arr = ones(1, 9, 20)
        timeidx = collect(1:20)
        params = ComponentVector(params = (lag = fill(0.2, 3),))
        initstates = ComponentVector(s_river = fill(0.0, length(positions)))
        
        config = create_test_config(solver = ImmutableSolver, timeidx = timeidx)
        output_arr = grid_route(input_arr, params, config; initstates = initstates)
        
        @test size(output_arr) == (2, 9, 20)
    end
end

@testset "Vector-based routing" begin
    # Define variables and parameters
    @variables q q_routed s_river
    @parameters lag

    # Define network structure
    network = DiGraph(9)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)

    # Create routing component
    vector_route = @hydroroute begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag)
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(network)
        hru_types = [1, 1, 1, 2, 2, 2, 3, 3, 3]
    end

    @testset "Route interface" begin
        @test HydroModels.get_input_names(vector_route) == [:q]
        @test HydroModels.get_output_names(vector_route) == [:q_routed]
        @test HydroModels.get_param_names(vector_route) == [:lag]
        @test HydroModels.get_state_names(vector_route) == [:s_river]
    end

    @testset "Run with different solvers" begin
        input_arr = rand(1, 9, 20)
        timeidx = collect(1:20)
        params = ComponentVector(params = (lag = fill(0.2, 3),))
        initstates = ComponentVector(s_river = fill(0.0, 9))
        
        # Test with MutableSolver
        config_mut = create_test_config(solver = MutableSolver, timeidx = timeidx)
        output_mut = vector_route(input_arr, params, config_mut; initstates = initstates)
        @test size(output_mut) == (2, 9, 20)
        
        # Test with ImmutableSolver
        config_immut = create_test_config(solver = ImmutableSolver, timeidx = timeidx)
        output_immut = vector_route(input_arr, params, config_immut; initstates = initstates)
        @test size(output_immut) == (2, 9, 20)
        
        # Results should be similar (allowing for numerical differences)
        @test output_mut ≈ output_immut atol = 1e-8
    end
    
    @testset "Flow accumulation correctness" begin
        # Create simple linear network: 1 -> 2 -> 3
        simple_network = DiGraph(3)
        add_edge!(simple_network, 1, 2)
        add_edge!(simple_network, 2, 3)
        
        simple_route = @hydroroute begin
            fluxes = begin
                @hydroflux q_routed ~ s_river / (1 + lag)
            end
            dfluxes = begin
                @stateflux s_river ~ q - q_routed
            end
            aggr_func = HydroModels.build_aggr_func(simple_network)
            hru_types = [1, 1, 1]
        end
        
        # Input: constant flow at each node
        input_arr = ones(1, 3, 10)
        params = ComponentVector(params = (lag = [0.1],))
        initstates = ComponentVector(s_river = fill(0.0, 3))
        
        config = create_test_config(solver = MutableSolver, timeidx = collect(1:10))
        output = simple_route(input_arr, params, config; initstates = initstates)
        
        # At steady state, outlet flow should be greater than upstream flows
        # (due to accumulation)
        @test output[2, 3, end] > output[2, 1, end]  # Outlet > Source
    end
end

@testset "Routing with flux input" begin
    # Test routing with q_gen input (generated flow)
    @variables q_gen q_routed s_river
    @parameters lag

    network = DiGraph(3)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 3)
    
    route_with_gen = @hydroroute begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag) + q_gen
        end
        dfluxes = begin
            @stateflux s_river ~ q_gen - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(network)
        hru_types = collect(1:3)
    end
    
    @testset "Route with generated flow" begin
        input_arr = rand(1, 3, 15)
        params = ComponentVector(params = (lag = fill(0.3, 3),))
        initstates = ComponentVector(s_river = fill(0.0, 3))
        
        config = create_test_config(solver = MutableSolver, timeidx = collect(1:15))
        output = route_with_gen(input_arr, params, config; initstates = initstates)
        
        @test size(output) == (2, 3, 15)
        @test all(output .>= -1e-10)  # Should be non-negative (allowing for numerical error)
    end
end
