@testset "M50 model (with neural network)" begin
    # Define parameters
    @parameters Tmin Tmax Df
    @parameters snowpack_std snowpack_mean soilwater_std soilwater_mean
    @parameters prcp_std prcp_mean temp_std temp_mean

    # Define variables
    @variables prcp temp lday pet rainfall snowfall snowpack soilwater
    @variables melt log_evap_div_lday log_flow
    @variables norm_snw norm_slw norm_temp norm_prcp

    # Load data
    df = DataFrame(CSV.File(joinpath(dirname(dirname(@__DIR__)), "data", "m50", "01013500.csv")))
    ts = collect(1:10000)
    prcp_vec = df[ts, "Prcp"]
    temp_vec = df[ts, "Temp"]
    dayl_vec = df[ts, "Lday"]
    snowpack_vec = df[ts, "SnowWater"]
    soilwater_vec = df[ts, "SoilWater"]

    # Calculate normalization parameters
    inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
    means, stds = mean.(inputs), std.(inputs)

    # Define snow bucket
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

    # Define soil bucket with neural networks
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
    end

    # Define complete model
    m50_model = @hydromodel :m50 begin
        snow_bucket
        soil_bucket
    end

    @testset "Model interface" begin
        @test Set(HydroModels.get_input_names(m50_model)) == Set([:prcp, :temp, :lday])
        @test Set(HydroModels.get_state_names(m50_model)) == Set([:snowpack, :soilwater])
        @test Set(HydroModels.get_nn_names(m50_model)) == Set([:etnn, :qnn])
    end

    @testset "Run model" begin
        # Prepare parameters
        base_params = (Df = 2.674, Tmax = 0.17, Tmin = -2.09)
        var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
        var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
        nn_params = (etnn = et_nn_p, qnn = q_nn_p)
        params = reduce(merge, [base_params, var_means, var_stds])

        initstates = ComponentVector(snowpack = 0.0, soilwater = 1303.0)
        pas = ComponentVector(params = params, nns = nn_params)

        input_ntp = (prcp = prcp_vec, lday = dayl_vec, temp = temp_vec)
        input_mat = Matrix(reduce(hcat, collect(input_ntp[[:prcp, :temp, :lday]]))')

        config = create_test_config(solver = MutableSolver, timeidx = ts)
        result_mat = m50_model(input_mat, pas, config; initstates = initstates)

        expected_n_outputs = length(HydroModels.get_state_names(m50_model)) +
                            length(HydroModels.get_output_names(m50_model))
        @test size(result_mat) == (expected_n_outputs, length(ts))
    end
end
