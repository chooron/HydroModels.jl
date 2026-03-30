# Test complete single-node lumped hydrological models

@testset "ExpHydro model (no neural network, no unit hydrograph)" begin
    # Define variables and parameters
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow

    # Load data
    ts = collect(1:100)
    input_ntp, input_mat, df = load_test_data(:exphydro, ts)

    # Setup parameters and states
    params = ComponentVector(params = ComponentVector(EXPHYDRO_PARAMS))
    initstates = ComponentVector(EXPHYDRO_STATES)

    # Define snow bucket
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
    end

    # Define soil bucket
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
    end

    # Define complete model
    exphydro_model = @hydromodel :exphydro begin
        snow_bucket
        soil_bucket
    end

    @testset "Model interface" begin
        @test Set(HydroModels.get_input_names(exphydro_model)) == Set([:temp, :lday, :prcp])
        @test Set(HydroModels.get_param_names(exphydro_model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax])
        @test Set(HydroModels.get_state_names(exphydro_model)) == Set([:snowpack, :soilwater])
        @test Set(HydroModels.get_output_names(exphydro_model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow])
    end

    @testset "Run with MutableSolver" begin
        config = create_test_config(solver = MutableSolver)
        result_mat = exphydro_model(input_mat, params, config; initstates = initstates)
        
        expected_n_outputs = length(HydroModels.get_state_names(exphydro_model)) + 
                            length(HydroModels.get_output_names(exphydro_model))
        @test size(result_mat) == (expected_n_outputs, length(ts))
        
        # Sanity checks
        @test all(result_mat[1, :] .>= 0)  # snowpack >= 0
        @test all(result_mat[2, :] .>= 0)  # soilwater >= 0
    end
    
    @testset "Run with ImmutableSolver" begin
        config = create_test_config(solver = ImmutableSolver)
        result_mat = exphydro_model(input_mat, params, config; initstates = initstates)
        
        expected_n_outputs = length(HydroModels.get_state_names(exphydro_model)) + 
                            length(HydroModels.get_output_names(exphydro_model))
        @test size(result_mat) == (expected_n_outputs, length(ts))
    end
    
    @testset "Solver consistency" begin
        # Both solvers should give similar results
        config_mut = create_test_config(solver = MutableSolver)
        config_immut = create_test_config(solver = ImmutableSolver)
        
        result_mut = exphydro_model(input_mat, params, config_mut; initstates = initstates)
        result_immut = exphydro_model(input_mat, params, config_immut; initstates = initstates)
        
        @test result_mut ≈ result_immut atol = 1e-8
    end
end

@testset "GR4J model (with unit hydrograph)" begin
    # Define variables and parameters
    @variables prcp ep soilwater pn en ps es perc pr slowflow fastflow t
    @variables slowflow_routed fastflow_routed routingstore exch routedflow flow
    @parameters x1 x2 x3 x4

    # Load data
    ts = collect(1:100)
    input_ntp, input_mat, df = load_test_data(:gr4j, ts)

    # Setup parameters and states
    params = ComponentVector(params = ComponentVector(GR4J_PARAMS))
    initstates = ComponentVector(GR4J_STATES)

    # Define production store
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
    end

    # Define unit hydrographs
    uh_slow = @unithydro :uh_slow begin
        uh_func = begin
            x4 => (t / x4)^2.5
        end
        uh_vars = slowflow => slowflow_routed
    end

    uh_fast = @unithydro :uh_fast begin
        uh_func = begin
            2x4 => (1 - 0.5 * (2 - t / x4)^2.5)
            x4 => (0.5 * (t / x4)^2.5)
        end
        uh_vars = fastflow => fastflow_routed
    end

    # Define routing store
    routing_bucket = @hydrobucket begin
        fluxes = begin
            @hydroflux exch ~ x2 * abs(routingstore / x3)^3.5
            @hydroflux routedflow ~ x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5
            @hydroflux flow ~ routedflow + max(fastflow_routed + exch, 0.0)
        end
        dfluxes = begin
            @stateflux routingstore ~ slowflow_routed + exch - routedflow
        end
    end

    # Define complete model
    gr4j_model = @hydromodel :gr4j begin
        prod_bucket
        uh_slow
        uh_fast
        routing_bucket
    end

    @testset "Model interface" begin
        @test Set(HydroModels.get_input_names(gr4j_model)) == Set([:prcp, :ep])
        @test Set(HydroModels.get_param_names(gr4j_model)) == Set([:x1, :x2, :x3, :x4])
        @test Set(HydroModels.get_state_names(gr4j_model)) == Set([:soilwater, :routingstore])
    end

    @testset "Run model" begin
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result_mat = gr4j_model(input_mat, params, config; initstates = initstates)
        
        expected_n_outputs = length(HydroModels.get_state_names(gr4j_model)) + 
                            length(HydroModels.get_output_names(gr4j_model))
        @test size(result_mat) == (expected_n_outputs, length(ts))
        
        # Sanity checks
        @test all(result_mat[1, :] .>= 0)  # soilwater >= 0
        @test all(result_mat[2, :] .>= 0)  # routingstore >= 0
    end
end
