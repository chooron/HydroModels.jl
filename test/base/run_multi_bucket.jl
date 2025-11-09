# Test multi-node bucket components

# Load test data
ts = collect(1:10)
input_ntp, input, df = load_test_data(:exphydro, ts)

# Define variables and parameters
@variables temp lday prcp pet snowfall rainfall melt snowpack
@parameters Tmin Tmax Df

# Number of nodes for testing
const NUM_NODES = 10

@testset "Multi hydro bucket with state" begin
    # Define single-node bucket for comparison
    snow_single = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
    end

    # Define multi-node bucket with independent parameters
    snow_multi_independent = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
        hru_types = collect(1:NUM_NODES)
    end

    # Define multi-node bucket with shared parameters
    snow_multi_shared = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
        hru_types = [1, 2, 2, 2, 1, 3, 3, 2, 3, 2]  # 3 parameter types
    end

    @testset "Multi-node with independent parameters" begin
        # Prepare multi-node inputs
        input_arr = create_multinode_input(input, NUM_NODES)
        
        # Prepare parameters (all nodes have same parameter values)
        single_params = ComponentVector(params = (Df = EXPHYDRO_PARAMS.Df, Tmax = EXPHYDRO_PARAMS.Tmax, Tmin = EXPHYDRO_PARAMS.Tmin))
        node_params = create_multinode_params(single_params, NUM_NODES)
        node_states = create_multinode_states(ComponentVector(snowpack = 0.0), NUM_NODES)
        
        # Run multi-node model
        config = create_test_config()
        node_output = snow_multi_independent(input_arr, node_params; initstates = node_states)
        
        # Run single-node model for comparison
        single_output = snow_single(input, single_params; initstates = ComponentVector(snowpack = 0.0), timeidx = ts)
        
        # Expected output: replicate single-node output for all nodes
        target_output = permutedims(
            reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], NUM_NODES)),
            (1, 3, 2)
        )
        
        @test node_output ≈ target_output atol = 1e-10
        @test size(node_output) == (
            length(HydroModels.get_state_names(snow_single)) + 
            length(HydroModels.get_output_names(snow_single)),
            NUM_NODES,
            size(input, 2)
        )
    end

    @testset "Multi-node with shared parameters" begin
        # Prepare multi-node inputs
        input_arr = create_multinode_input(input, NUM_NODES)
        
        # Prepare parameters (3 parameter types)
        single_params = ComponentVector(params = (Df = EXPHYDRO_PARAMS.Df, Tmax = EXPHYDRO_PARAMS.Tmax, Tmin = EXPHYDRO_PARAMS.Tmin))
        node_params = ComponentVector(params = (
            Df = fill(EXPHYDRO_PARAMS.Df, 3),
            Tmax = fill(EXPHYDRO_PARAMS.Tmax, 3),
            Tmin = fill(EXPHYDRO_PARAMS.Tmin, 3)
        ))
        node_states = create_multinode_states(ComponentVector(snowpack = 0.0), NUM_NODES)
        
        # Run multi-node model
        config = create_test_config()
        node_output = snow_multi_shared(input_arr, node_params; initstates = node_states)
        
        # Run single-node model for comparison
        single_output = snow_single(input, single_params; initstates = ComponentVector(snowpack = 0.0), timeidx = ts)
        
        # Expected output: same for all nodes (since all parameter types have same values)
        target_output = permutedims(
            reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], NUM_NODES)),
            (1, 3, 2)
        )
        
        @test node_output ≈ target_output atol = 1e-10
        @test size(node_output) == (
            length(HydroModels.get_state_names(snow_single)) + 
            length(HydroModels.get_output_names(snow_single)),
            NUM_NODES,
            size(input, 2)
        )
    end
end

@testset "Multi hydro bucket without state" begin
    # Define single-node bucket without state
    snow_single = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
    end

    # Define multi-node buckets
    snow_multi_independent = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
        hru_types = collect(1:NUM_NODES)
    end
    
    snow_multi_shared = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
        hru_types = [1, 2, 3, 3, 2, 1, 2, 1, 3, 2]  # 3 parameter types
    end

    @testset "Multi-node with independent parameters" begin
        input_arr = create_multinode_input(input, NUM_NODES)
        single_params = ComponentVector(params = (Tmin = EXPHYDRO_PARAMS.Tmin,))
        node_params = create_multinode_params(single_params, NUM_NODES)
        node_states = ComponentVector()
        
        config = create_test_config()
        node_output = snow_multi_independent(input_arr, node_params; initstates = node_states)
        
        # Compare with single-node
        single_output = snow_single(input, single_params; initstates = node_states, timeidx = ts)
        target_output = permutedims(
            reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], NUM_NODES)),
            (1, 3, 2)
        )
        
        @test node_output ≈ target_output atol = 1e-10
    end

    @testset "Multi-node with shared parameters" begin
        input_arr = create_multinode_input(input, NUM_NODES)
        node_params = ComponentVector(params = (Tmin = fill(EXPHYDRO_PARAMS.Tmin, 3),))
        node_states = ComponentVector()
        
        config = create_test_config()
        node_output = snow_multi_shared(input_arr, node_params; initstates = node_states)
        
        # Compare with single-node
        single_params = ComponentVector(params = (Tmin = EXPHYDRO_PARAMS.Tmin,))
        single_output = snow_single(input, single_params; initstates = node_states, timeidx = ts)
        target_output = permutedims(
            reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], NUM_NODES)),
            (1, 3, 2)
        )
        
        @test node_output ≈ target_output atol = 1e-10
    end
end
