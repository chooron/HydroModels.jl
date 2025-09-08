step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

@testset "test lumped hydro model (exp-hydro with no neural network and no unit hydrograph)" begin
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow

    #! load data
    ts = collect(1:100)
    df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input_mat = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

    f_value, Smax_value, Qmax_value, Df_value, Tmax_value, Tmin_value = 0.0167, 1709.46, 18.47, 2.674, 0.17, -2.09
    snowpack_value, soilwater_value = 0.0, 1303.00

    initstates = ComponentVector(snowpack=snowpack_value, soilwater=soilwater_value)
    params = ComponentVector(f=f_value, Smax=Smax_value, Qmax=Qmax_value, Df=Df_value, Tmax=Tmax_value, Tmin=Tmin_value)
    pas = ComponentVector(params=params)

    bucket_1 = @hydrobucket :surface begin
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

    bucket_2 = @hydrobucket :soil begin
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

    model = @hydromodel :exphydro begin
        bucket_1
        bucket_2
    end

    @test Set(HydroModels.get_input_names(model)) == Set([:temp, :lday, :prcp])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow])
    @test Set(reduce(union, HydroModels.get_var_names(model))) == Set([:temp, :lday, :prcp, :pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :snowpack, :soilwater])

    result_mat = model(input_mat, pas, initstates=initstates)
    @test size(result_mat) == (length(HydroModels.get_state_names(model)) + length(HydroModels.get_output_names(model)), length(ts))
end

@testset "test lumped hydro model (gr4j with unit hydrograph)" begin
    @variables prcp ep soilwater new_soilwater pn en ps es perc pr slowflow fastflow t
    @variables slowflow_routed fastflow_routed routingstore new_routingstore exch routedflow flow
    @parameters x1 x2 x3 x4

    #! load data
    # load data
    df = DataFrame(CSV.File("../data/gr4j/sample.csv"))
    for col in names(df)[3:end]
        df[ismissing.(df[:, col]), col] .= 0.0
    end
    prcp_vec = df[!, "prec"]
    et_vec = df[!, "pet"]
    qobs_vec = df[!, "qobs"]
    ts = collect(1:length(qobs_vec))
    input_ntp = (prcp=prcp_vec, ep=et_vec)
    input_mat = Matrix(reduce(hcat, collect(input_ntp[[:prcp, :ep]]))')

    x1_value, x2_value, x3_value, x4_value = 320.11, 2.42, 69.63, 1.39
    soilwater_value, routingstore_value = 235.97, 45.47

    params = ComponentVector(x1=x1_value, x2=x2_value, x3=x3_value, x4=x4_value)
    initstates = ComponentVector(soilwater=soilwater_value, routingstore=routingstore_value)
    pas = ComponentVector(params=params)

    #* define the production store
    prod_ele = @hydrobucket begin
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


    uh_1 = @unithydro :maxbas_uh begin
        uh_func = begin
            x4 => (t / x4)^2.5
        end
        uh_vars = slowflow => slowflow_routed
    end

    uh_2 = @unithydro begin
        uh_func = begin
            2x4 => (1 - 0.5 * (2 - t / x4)^2.5)
            x4 => (0.5 * (t / x4)^2.5)
        end
        uh_vars = fastflow => fastflow_routed
    end

    rst_ele = @hydrobucket begin
        fluxes = begin
            @hydroflux exch ~ x2 * abs(routingstore / x3)^3.5
            @hydroflux routedflow ~ x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5
            @hydroflux flow ~ routedflow + max(fastflow_routed + exch, 0.0)
        end
        dfluxes = begin
            @stateflux routingstore ~ slowflow_routed + exch - routedflow
        end
    end

    #* define the gr4j model
    model = @hydromodel :gr4j begin
        prod_ele
        uh_1
        uh_2
        rst_ele
    end

    @test Set(HydroModels.get_input_names(model)) == Set([:prcp, :ep])
    @test Set(HydroModels.get_param_names(model)) == Set([:x1, :x2, :x3, :x4])
    @test Set(HydroModels.get_state_names(model)) == Set([:soilwater, :routingstore])
    @test Set(HydroModels.get_output_names(model)) == Set([:en, :routedflow, :pr, :exch, :pn, :fastflow, :ps, :flow, :slowflow_routed, :perc,
        :es, :slowflow, :fastflow_routed])
    @test Set(reduce(union, HydroModels.get_var_names(model))) == Set([:prcp, :ep, :soilwater, :pn, :en, :ps, :es, :perc, :pr, :slowflow,
        :fastflow, :slowflow_routed, :fastflow_routed, :exch, :routedflow, :flow, :routingstore])

    # Test single-node model run
    result_mat = model(input_mat, pas, initstates=initstates, config=(timeidx=ts,))
    @test size(result_mat) == (length(HydroModels.get_state_names(model)) + length(HydroModels.get_output_names(model)), length(ts))
end


@testset "test lumped hydro model (m50 with neural network)" begin
    #! parameters in the Exp-Hydro model
    @parameters Tmin Tmax Df Smax f Qmax
    #! parameters in normalize flux
    @parameters snowpack_std snowpack_mean
    @parameters soilwater_std soilwater_mean
    @parameters prcp_std prcp_mean
    @parameters temp_std temp_mean

    #! hydrological flux in the Exp-Hydro model
    @variables prcp temp lday pet rainfall snowfall
    @variables snowpack soilwater lday pet
    @variables melt log_evap_div_lday log_flow
    @variables norm_snw norm_slw norm_temp norm_prcp

    #! load data
    df = DataFrame(CSV.File("../data/m50/01013500.csv"))
    ts = collect(1:10000)
    prcp_vec = df[ts, "Prcp"]
    temp_vec = df[ts, "Temp"]
    dayl_vec = df[ts, "Lday"]
    snowpack_vec = df[ts, "SnowWater"]
    soilwater_vec = df[ts, "SoilWater"]
    qobs_vec = df[ts, "Flow"]

    inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
    means, stds = mean.(inputs), std.(inputs)
    (prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec) = [@.((vec - mean) / std) for (vec, mean, std) in zip(inputs, means, stds)]

    #! define the snow pack reservoir
    snow_ele = @hydrobucket :m50_snow begin
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

    #! define the ET NN and Q NN
    et_nn = Lux.Chain(Lux.Dense(3 => 16, Lux.tanh), Lux.Dense(16 => 16, Lux.leakyrelu), Lux.Dense(16 => 1, Lux.leakyrelu), name=:etnn)
    et_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), et_nn))
    q_nn = Lux.Chain(Lux.Dense(2 => 16, Lux.tanh), Lux.Dense(16 => 16, Lux.leakyrelu), Lux.Dense(16 => 1, Lux.leakyrelu), name=:qnn)
    q_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn))

    #! define the soil water reservoir
    soil_ele = @hydrobucket :m50_soil begin
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

    #! define the Exp-Hydro model
    model = @hydromodel :m50 begin
        snow_ele
        soil_ele
    end

    @test Set(HydroModels.get_input_names(model)) == Set([:prcp, :temp, :lday])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :snowpack_std, :snowpack_mean, :soilwater_std, :soilwater_mean, :prcp_std, :prcp_mean, :temp_std, :temp_mean])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater])
    @test Set(HydroModels.get_nn_names(model)) == Set([:etnn, :qnn])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :rainfall, :snowfall, :melt, :log_evap_div_lday, :log_flow, :norm_snw, :norm_slw, :norm_temp, :norm_prcp])
    @test Set(reduce(union, HydroModels.get_var_names(model))) == Set([:prcp, :temp, :lday, :pet, :rainfall, :snowfall, :snowpack, :soilwater, :melt,
        :log_evap_div_lday, :log_flow, :norm_snw, :norm_slw, :norm_temp, :norm_prcp])

    base_params = (Df=2.674, Tmax=0.17, Tmin=-2.09)
    var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
    var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
    nn_params = (etnn=et_nn_p, qnn=q_nn_p)
    params = reduce(merge, [base_params, var_means, var_stds])
    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00)
    pas = ComponentVector(params=params, nns=nn_params)
    input_ntp = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
    input_mat = Matrix(reduce(hcat, collect(input_ntp[[:prcp, :temp, :lday]]))')

    # Run the model and get results as a matrix
    result_mat = model(input_mat, pas, initstates=initstates, config=(timeidx=ts,))
    @test size(result_mat) == (length(HydroModels.get_state_names(model)) + length(HydroModels.get_output_names(model)), length(ts))
end