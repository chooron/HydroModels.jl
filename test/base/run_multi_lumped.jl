# Test complete multi-node lumped hydrological models

const NUM_TEST_NODES = 10

@testset "Multi-node ExpHydro model" begin
    # Define variables and parameters
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow

    # Load data
    ts = collect(1:100)
    input_ntp, input_mat, df = load_test_data(:exphydro, ts)

    # Define snow bucket with hru_types
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
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define soil bucket with hru_types
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
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define complete model
    multi_exphydro = @hydromodel :multi_exphydro begin
        snow_bucket
        soil_bucket
    end

    @testset "Multi-node with independent parameters" begin
        # Prepare multi-node inputs
        input_arr = create_multinode_input(input_mat, NUM_TEST_NODES)
        
        # Prepare parameters (all nodes use same values)
        single_params = ComponentVector(params = ComponentVector(EXPHYDRO_PARAMS))
        node_params = create_multinode_params(single_params, NUM_TEST_NODES)
        node_states = create_multinode_states(ComponentVector(EXPHYDRO_STATES), NUM_TEST_NODES)
        
        # Run model
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result_arr = multi_exphydro(input_arr, node_params, config; initstates = node_states)
        
        expected_n_outputs = length(HydroModels.get_state_names(multi_exphydro)) + 
                            length(HydroModels.get_output_names(multi_exphydro))
        @test size(result_arr) == (expected_n_outputs, NUM_TEST_NODES, length(ts))
        
        # Sanity checks
        @test all(result_arr[1, :, :] .>= 0)  # snowpack >= 0
        @test all(result_arr[2, :, :] .>= 0)  # soilwater >= 0
    end
end

@testset "Multi-node GR4J model" begin
    # Define variables and parameters
    @variables prcp ep soilwater pn en ps es perc pr slowflow fastflow t
    @variables slowflow_routed fastflow_routed routingstore exch routedflow flow
    @parameters x1 x2 x3 x4

    # Load data
    ts = collect(1:100)
    input_ntp, input_mat, df = load_test_data(:gr4j, ts)

    # Define production bucket with hru_types
    prod_bucket = @hydrobucket begin
        fluxes = begin
            @hydroflux pn ~ prcp - min(prcp, ep)
            @hydroflux en ~ ep - min(prcp, ep)
            @hydroflux ps ~ max(0.0, pn * (1 - (soilwater / x1)^2))
            @hydroflux es ~ en * (2 * soilwater / x1 - (soilwater / x1)^2)
            @hydroflux perc ~ ((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)
            @hydroflux pr ~ pn - ps + perc
            @hydroflux slowflow ~ 0.9 * pr
            @hydroflux fastflow ~ 0.1 * pr
        end
        dfluxes = begin
            @stateflux soilwater ~ ps - es - perc
        end
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define unit hydrographs with hru_types
    uh_slow = @unithydro :uh_slow begin
        uh_func = begin
            x4 => (t / x4)^2.5
        end
        uh_vars = slowflow => slowflow_routed
        hru_types = collect(1:NUM_TEST_NODES)
    end

    uh_fast = @unithydro :uh_fast begin
        uh_func = begin
            2x4 => (1 - 0.5 * (2 - t / x4)^2.5)
            x4 => (0.5 * (t / x4)^2.5)
        end
        uh_vars = fastflow => fastflow_routed
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define routing bucket with hru_types
    routing_bucket = @hydrobucket begin
        fluxes = begin
            @hydroflux exch ~ x2 * abs(routingstore / x3)^3.5
            @hydroflux routedflow ~ x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5
            @hydroflux flow ~ routedflow + max(fastflow_routed + exch, 0.0)
        end
        dfluxes = begin
            @stateflux routingstore ~ slowflow_routed + exch - routedflow
        end
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define complete model
    multi_gr4j = @hydromodel :multi_gr4j begin
        prod_bucket
        uh_slow
        uh_fast
        routing_bucket
    end

    @testset "Multi-node GR4J execution" begin
        # Prepare multi-node inputs
        input_arr = create_multinode_input(input_mat, NUM_TEST_NODES)
        
        # Prepare parameters
        single_params = ComponentVector(params = ComponentVector(GR4J_PARAMS))
        node_params = create_multinode_params(single_params, NUM_TEST_NODES)
        node_states = create_multinode_states(ComponentVector(GR4J_STATES), NUM_TEST_NODES)
        
        # Run model
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result_arr = multi_gr4j(input_arr, node_params, config; initstates = node_states)
        
        expected_n_outputs = length(HydroModels.get_state_names(multi_gr4j)) + 
                            length(HydroModels.get_output_names(multi_gr4j))
        @test size(result_arr) == (expected_n_outputs, NUM_TEST_NODES, length(ts))
    end
end

@testset "Multi-node M50 model" begin
    # Define parameters
    @parameters Tmin Tmax Df
    @parameters snowpack_std snowpack_mean soilwater_std soilwater_mean
    @parameters prcp_std prcp_mean temp_std temp_mean

    # Define variables
    @variables prcp temp lday pet rainfall snowfall snowpack soilwater
    @variables melt log_evap_div_lday log_flow
    @variables norm_snw norm_slw norm_temp norm_prcp

    # Load data
    df = DataFrame(CSV.File("../data/m50/01013500.csv"))
    ts = collect(1:10000)
    prcp_vec = df[ts, "Prcp"]
    temp_vec = df[ts, "Temp"]
    dayl_vec = df[ts, "Lday"]
    snowpack_vec = df[ts, "SnowWater"]
    soilwater_vec = df[ts, "SoilWater"]

    # Calculate normalization parameters
    inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
    means, stds = mean.(inputs), std.(inputs)

    # Define snow bucket with hru_types
    snow_bucket = @hydrobucket :m50_snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define neural networks
    et_nn = Lux.Chain(
        Lux.Dense(3 => 16, Lux.tanh),
        Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 1, Lux.leakyrelu),
        name = :etnn
    )
    et_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), et_nn))
    
    q_nn = Lux.Chain(
        Lux.Dense(2 => 16, Lux.tanh),
        Lux.Dense(16 => 16, Lux.leakyrelu),
        Lux.Dense(16 => 1, Lux.leakyrelu),
        name = :qnn
    )
    q_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn))

    # Define soil bucket with neural networks and hru_types
    soil_bucket = @hydrobucket :m50_soil begin
        fluxes = begin
            @hydroflux norm_snw ~ (snowpack - snowpack_mean) / snowpack_std
            @hydroflux norm_slw ~ (soilwater - soilwater_mean) / soilwater_std
            @hydroflux norm_prcp ~ (prcp - prcp_mean) / prcp_std
            @hydroflux norm_temp ~ (temp - temp_mean) / temp_std
            @neuralflux log_evap_div_lday ~ et_nn([norm_snw, norm_slw, norm_temp])
            @neuralflux log_flow ~ q_nn([norm_slw, norm_prcp])
        end
        dfluxes = begin
            @stateflux soilwater ~ rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
        end
        hru_types = collect(1:NUM_TEST_NODES)
    end

    # Define complete model
    multi_m50 = @hydromodel :multi_m50 begin
        snow_bucket
        soil_bucket
    end

    @testset "Multi-node M50 execution" begin
        # Prepare multi-node inputs
        input_ntp = (prcp = prcp_vec, lday = dayl_vec, temp = temp_vec)
        input_mat = Matrix(reduce(hcat, collect(input_ntp[[:prcp, :temp, :lday]]))')
        input_arr = create_multinode_input(input_mat, NUM_TEST_NODES)
        
        # Prepare parameters
        base_params = (Df = 2.674, Tmax = 0.17, Tmin = -2.09)
        var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
        var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
        nn_params = (etnn = et_nn_p, qnn = q_nn_p)
        
        # Create single-node params first
        single_params_dict = reduce(merge, [base_params, var_means, var_stds])
        
        # Expand to multi-node (all nodes use same parameter values)
        node_params_dict = NamedTuple{keys(single_params_dict)}(
            fill(v, NUM_TEST_NODES) for (k, v) in pairs(single_params_dict)
        )
        
        pas = ComponentVector(params = node_params_dict, nns = nn_params)
        node_states = create_multinode_states(ComponentVector(snowpack = 0.0, soilwater = 1303.0), NUM_TEST_NODES)
        
        # Run model
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result_arr = multi_m50(input_arr, pas, config; initstates = node_states)
        
        expected_n_outputs = length(HydroModels.get_state_names(multi_m50)) + 
                            length(HydroModels.get_output_names(multi_m50))
        @test size(result_arr) == (expected_n_outputs, NUM_TEST_NODES, length(ts))
    end
end
