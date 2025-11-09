# Test single-node bucket components

# Load test data
ts = collect(1:10)
input_ntp, input, df = load_test_data(:exphydro, ts)

# Define variables and parameters
@variables temp lday prcp pet snowfall rainfall melt snowpack
@parameters Tmin Tmax Df

# Setup test parameters and states
params = ComponentVector(params = ComponentVector(Df = EXPHYDRO_PARAMS.Df, Tmax = EXPHYDRO_PARAMS.Tmax, Tmin = EXPHYDRO_PARAMS.Tmin))
init_states = ComponentVector(snowpack = 0.0)

@testset "Single hydro bucket with state" begin
    # Define snow bucket
    snow_bucket = @hydrobucket :snow begin
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

    @testset "Bucket interface" begin
        @test Set(HydroModels.get_input_names(snow_bucket)) == Set([:temp, :lday, :prcp])
        @test Set(HydroModels.get_param_names(snow_bucket)) == Set([:Tmin, :Tmax, :Df])
        @test Set(HydroModels.get_output_names(snow_bucket)) == Set([:pet, :snowfall, :rainfall, :melt])
        @test Set(HydroModels.get_state_names(snow_bucket)) == Set([:snowpack])
    end

    @testset "Run with MutableSolver" begin
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result = snow_bucket(input, params, config; initstates = init_states)
        
        @test size(result) == (
            length(HydroModels.get_state_names(snow_bucket)) + 
            length(HydroModels.get_output_names(snow_bucket)), 
            length(ts)
        )
        
        # Convert to NamedTuple for easier inspection
        output_names = vcat(HydroModels.get_state_names(snow_bucket), HydroModels.get_output_names(snow_bucket))
        result_nt = NamedTuple{Tuple(output_names)}(eachslice(result, dims=1))
        
        # Basic sanity checks
        @test all(result_nt.snowpack .>= 0)  # Snowpack should be non-negative
        @test all(result_nt.pet .>= 0)  # PET should be non-negative
    end
    
    @testset "Run with ImmutableSolver" begin
        config = create_test_config(solver = ImmutableSolver, timeidx = ts)
        result = snow_bucket(input, params, config; initstates = init_states)
        
        @test size(result) == (
            length(HydroModels.get_state_names(snow_bucket)) + 
            length(HydroModels.get_output_names(snow_bucket)), 
            length(ts)
        )
    end
end

@testset "Single hydro bucket without state" begin
    # Define bucket without state
    snow_bucket_no_state = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
    end

    @testset "Bucket interface" begin
        @test Set(HydroModels.get_input_names(snow_bucket_no_state)) == Set([:temp, :lday, :prcp])
        @test Set(HydroModels.get_param_names(snow_bucket_no_state)) == Set([:Tmin])
        @test Set(HydroModels.get_output_names(snow_bucket_no_state)) == Set([:pet, :snowfall, :rainfall])
        @test isempty(HydroModels.get_state_names(snow_bucket_no_state))
    end

    @testset "Run without config timeidx" begin
        config = create_test_config(solver = MutableSolver)
        result = snow_bucket_no_state(input, params, config; initstates = init_states)
        
        @test size(result) == (
            length(HydroModels.get_output_names(snow_bucket_no_state)), 
            size(input, 2)
        )
    end
    
    @testset "Run with config timeidx" begin
        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result = snow_bucket_no_state(input, params, config; initstates = init_states)
        
        @test size(result) == (
            length(HydroModels.get_output_names(snow_bucket_no_state)), 
            length(ts)
        )
    end
end
