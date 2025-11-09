# Test spatial hydrological models with routing

# Common setup for spatial models
@parameters Tmin Tmax Df Smax f Qmax area_coef lag
@variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow q q_routed s_river

# Define base buckets with hru_types for spatial use
const NUM_SPATIAL_NODES = 9

snow_bucket = @hydrobucket :surface begin
    fluxes = begin
        @hydroflux begin
            snowfall ~ step_func(Tmin - temp) * prcp
            rainfall ~ step_func(temp - Tmin) * prcp
        end
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
    hru_types = collect(1:NUM_SPATIAL_NODES)
end

soil_bucket = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux evap ~ step_func(soilwater) * pet * min(1.0, soilwater / Smax)
        @hydroflux baseflow ~ step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))
        @hydroflux surfaceflow ~ max(0.0, soilwater - Smax)
        @hydroflux flow ~ baseflow + surfaceflow
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + melt) - (evap + flow)
    end
    hru_types = collect(1:NUM_SPATIAL_NODES)
end

# Convert flow to discharge
convert_flux = @hydroflux begin
    q ~ flow * area_coef
    hru_types = collect(1:NUM_SPATIAL_NODES)
end

@testset "Grid-based spatial routing model" begin
    # Define grid structure (3x3)
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]

    # Define routing component
    grid_route = @hydroroute begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag) + q
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(flwdir, positions)
        hru_types = collect(1:NUM_SPATIAL_NODES)
    end

    # Define complete model
    grid_spatial_model = @hydromodel :exphydro_grid begin
        snow_bucket
        soil_bucket
        convert_flux
        grid_route
    end

    @testset "Model interface" begin
        @test Set(HydroModels.get_input_names(grid_spatial_model)) == Set([:temp, :lday, :prcp])
        @test Set(HydroModels.get_param_names(grid_spatial_model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax, :lag, :area_coef])
        @test Set(HydroModels.get_state_names(grid_spatial_model)) == Set([:snowpack, :soilwater, :s_river])
        @test Set(HydroModels.get_output_names(grid_spatial_model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :q, :q_routed])
    end

    @testset "Run grid spatial model" begin
        # Load data and prepare inputs
        ts = collect(1:100)
        input_ntp, input_mat, df = load_test_data(:exphydro, ts)
        input_arr = create_multinode_input(input_mat, NUM_SPATIAL_NODES)
        
        # Prepare parameters
        single_params = ComponentVector(params = (
            Df = 2.674, Tmax = 0.17, Tmin = -2.09,
            Smax = 1709.46, f = 0.0167, Qmax = 18.47,
            lag = 0.2, area_coef = 1.0
        ))
        node_params = create_multinode_params(single_params, NUM_SPATIAL_NODES)
        
        # Prepare states
        node_states = ComponentVector(
            snowpack = fill(0.0, NUM_SPATIAL_NODES),
            soilwater = fill(0.0, NUM_SPATIAL_NODES),
            s_river = fill(0.0, NUM_SPATIAL_NODES)
        )
        
        # Run model
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result = grid_spatial_model(input_arr, node_params, config; initstates = node_states)
        
        expected_n_outputs = length(HydroModels.get_state_names(grid_spatial_model)) + 
                            length(HydroModels.get_output_names(grid_spatial_model))
        @test size(result) == (expected_n_outputs, NUM_SPATIAL_NODES, length(ts))
        
        # Sanity checks
        @test all(result[1, :, :] .>= 0)  # snowpack >= 0
        @test all(result[2, :, :] .>= 0)  # soilwater >= 0
        @test all(result[3, :, :] .>= 0)  # s_river >= 0
    end
    
    @testset "Flow routing correctness" begin
        # Test that outlet receives more flow than sources (due to accumulation)
        ts = collect(1:50)
        input_ntp, input_mat, df = load_test_data(:exphydro, ts)
        input_arr = create_multinode_input(input_mat, NUM_SPATIAL_NODES)
        
        single_params = ComponentVector(params = (
            Df = 2.674, Tmax = 0.17, Tmin = -2.09,
            Smax = 1709.46, f = 0.0167, Qmax = 18.47,
            lag = 0.2, area_coef = 1.0
        ))
        node_params = create_multinode_params(single_params, NUM_SPATIAL_NODES)
        
        node_states = ComponentVector(
            snowpack = fill(0.0, NUM_SPATIAL_NODES),
            soilwater = fill(1303.0, NUM_SPATIAL_NODES),
            s_river = fill(0.0, NUM_SPATIAL_NODES)
        )
        
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result = grid_spatial_model(input_arr, node_params, config; initstates = node_states)
        
        # Output index for q_routed (last output)
        q_routed_idx = expected_n_outputs
        
        # Node 9 (3,3) is the outlet - should accumulate flow
        # Node 1 (1,1) is a source - should have local flow only
        @test mean(result[q_routed_idx, 9, :]) > mean(result[q_routed_idx, 1, :])
    end
end

@testset "Vector-based spatial routing model" begin
    # Define network structure
    network = DiGraph(NUM_SPATIAL_NODES)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)

    # Define routing component
    vector_route = @hydroroute :exphydro_routed begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag) + q
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(network)
        hru_types = collect(1:NUM_SPATIAL_NODES)
    end

    # Define complete model
    vector_spatial_model = @hydromodel :exphydro_vector begin
        snow_bucket
        soil_bucket
        convert_flux
        vector_route
    end

    @testset "Model interface" begin
        @test Set(HydroModels.get_input_names(vector_spatial_model)) == Set([:temp, :lday, :prcp])
        @test Set(HydroModels.get_param_names(vector_spatial_model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax, :lag, :area_coef])
        @test Set(HydroModels.get_state_names(vector_spatial_model)) == Set([:snowpack, :soilwater, :s_river])
    end

    @testset "Run vector spatial model" begin
        # Load data and prepare inputs
        ts = collect(1:100)
        input_ntp, input_mat, df = load_test_data(:exphydro, ts)
        input_arr = create_multinode_input(input_mat, NUM_SPATIAL_NODES)
        
        # Prepare parameters
        single_params = ComponentVector(params = (
            Df = 2.674, Tmax = 0.17, Tmin = -2.09,
            Smax = 1709.46, f = 0.0167, Qmax = 18.47,
            lag = 0.2, area_coef = 1.0
        ))
        node_params = create_multinode_params(single_params, NUM_SPATIAL_NODES)
        
        # Prepare states
        node_states = ComponentVector(
            snowpack = fill(0.0, NUM_SPATIAL_NODES),
            soilwater = fill(0.0, NUM_SPATIAL_NODES),
            s_river = fill(0.0, NUM_SPATIAL_NODES)
        )
        
        # Run model
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result = vector_spatial_model(input_arr, node_params, config; initstates = node_states)
        
        expected_n_outputs = length(HydroModels.get_state_names(vector_spatial_model)) + 
                            length(HydroModels.get_output_names(vector_spatial_model))
        @test size(result) == (expected_n_outputs, NUM_SPATIAL_NODES, length(ts))
        
        # Sanity checks
        @test all(result[1, :, :] .>= 0)  # snowpack >= 0
        @test all(result[2, :, :] .>= 0)  # soilwater >= 0
        @test all(result[3, :, :] .>= 0)  # s_river >= 0
    end
    
    @testset "Network flow accumulation" begin
        # Test flow accumulation along network path
        # Node 9 is outlet, should have highest accumulated flow at steady state
        ts = collect(1:100)
        input_ntp, input_mat, df = load_test_data(:exphydro, ts)
        input_arr = create_multinode_input(input_mat, NUM_SPATIAL_NODES)
        
        single_params = ComponentVector(params = (
            Df = 2.674, Tmax = 0.17, Tmin = -2.09,
            Smax = 1709.46, f = 0.0167, Qmax = 18.47,
            lag = 0.2, area_coef = 1.0
        ))
        node_params = create_multinode_params(single_params, NUM_SPATIAL_NODES)
        
        node_states = ComponentVector(
            snowpack = fill(0.0, NUM_SPATIAL_NODES),
            soilwater = fill(1303.0, NUM_SPATIAL_NODES),
            s_river = fill(0.0, NUM_SPATIAL_NODES)
        )
        
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result = vector_spatial_model(input_arr, node_params, config; initstates = node_states)
        
        # Outlet (node 9) should have higher flow than upstream nodes
        q_routed_idx = expected_n_outputs
        @test mean(result[q_routed_idx, 9, end-10:end]) > mean(result[q_routed_idx, 1, end-10:end])
    end
end

@testset "Spatial model solver consistency" begin
    # Test that different solvers give consistent results for spatial models
    network = DiGraph(NUM_SPATIAL_NODES)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)

    route = @hydroroute begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag) + q
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(network)
        hru_types = collect(1:NUM_SPATIAL_NODES)
    end

    spatial_model = @hydromodel :spatial_test begin
        snow_bucket
        soil_bucket
        convert_flux
        route
    end

    # Prepare test inputs
    ts = collect(1:50)
    input_ntp, input_mat, df = load_test_data(:exphydro, ts)
    input_arr = create_multinode_input(input_mat, NUM_SPATIAL_NODES)
    
    single_params = ComponentVector(params = (
        Df = 2.674, Tmax = 0.17, Tmin = -2.09,
        Smax = 1709.46, f = 0.0167, Qmax = 18.47,
        lag = 0.2, area_coef = 1.0
    ))
    node_params = create_multinode_params(single_params, NUM_SPATIAL_NODES)
    
    node_states = ComponentVector(
        snowpack = fill(0.0, NUM_SPATIAL_NODES),
        soilwater = fill(1303.0, NUM_SPATIAL_NODES),
        s_river = fill(0.0, NUM_SPATIAL_NODES)
    )

    # Run with both solvers
    config_mut = create_test_config(solver = MutableSolver, timeidx = ts)
    config_immut = create_test_config(solver = ImmutableSolver, timeidx = ts)
    
    result_mut = spatial_model(input_arr, node_params, config_mut; initstates = node_states)
    result_immut = spatial_model(input_arr, node_params, config_immut; initstates = node_states)
    
    # Results should be very similar
    @test result_mut ≈ result_immut atol = 1e-6
end
