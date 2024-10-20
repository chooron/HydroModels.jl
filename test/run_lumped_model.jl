@testset "test lumped hydro model (exp-hydro with no neural network and no unit hydrograph)" begin
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow

    #! load data
    ts = collect(1:100)
    df = DataFrame(CSV.File("data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00)
    params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09)
    pas = ComponentVector(initstates=initstates, params=params)

    #! define the snow pack reservoir
    snow_funcs = [
        SimpleFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroBucket(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

    #! define the soil water reservoir
    soil_funcs = [
        SimpleFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        SimpleFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        SimpleFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
        SimpleFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
    ]
    soil_dfuncs = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
    soil_ele = HydroBucket(:exphydro_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

    #! define the Exp-Hydro model
    model = HydroModel(:exphydro, components=[snow_ele, soil_ele])

    @test Set(HydroModels.get_input_names(model)) == Set([:temp, :lday, :prcp])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow])
    @test Set(HydroModels.get_var_names(model)) == Set([:temp, :lday, :prcp, :pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :snowpack, :soilwater])

    result_ntp = model(input_ntp, pas, timeidx=ts)
    @test result_ntp isa NamedTuple

    result_mat = model(input_ntp, pas, timeidx=ts, output_type=:matrix)
    @test size(result_mat) == (length(HydroModels.get_var_names(model)), length(ts))

    input_ntp_vec = repeat([input_ntp], 10)
    node_names = [Symbol(:node_, i) for i in 1:10]
    node_params = NamedTuple{Tuple(node_names)}(repeat([params], 10))
    node_initstates = NamedTuple{Tuple(node_names)}(repeat([initstates], 10))
    node_pas = ComponentVector(params=node_params, initstates=node_initstates)
    result_ntp_vec = model(input_ntp_vec, node_pas, timeidx=ts, output_type=:namedtuple)
    @test result_ntp_vec isa Vector
    @test result_ntp_vec[1] isa NamedTuple

    result_mat_vec = model(input_ntp_vec, node_pas, timeidx=ts, output_type=:array)
    @test size(result_mat_vec) == (length(HydroModels.get_var_names(model)), length(input_ntp_vec), length(ts))
end

@testset "test lumped hydro model (gr4j with unit hydrograph)" begin
    @variables prcp ep soilwater new_soilwater pn en ps es perc pr slowflow fastflow
    @variables slowflow_routed fastflow_routed routingstore new_routingstore exch routedflow flow
    @parameters x1 x2 x3 x4

    #! load data
    # load data
    df = DataFrame(CSV.File("data/gr4j/sample.csv"))
    for col in names(df)[3:end]
        df[ismissing.(df[:, col]), col] .= 0.0
    end
    prcp_vec = df[!, "prec"]
    et_vec = df[!, "pet"]
    qobs_vec = df[!, "qobs"]
    ts = collect(1:length(qobs_vec))
    input_ntp = (prcp=prcp_vec, ep=et_vec)

    params = ComponentVector(x1=320.11, x2=2.42, x3=69.63, x4=1.39)
    initstates = ComponentVector(soilwater=235.97, routingstore=45.47)
    pas = ComponentVector(initstates=initstates, params=params)

    #* define the production store
    prod_funcs = [
        SimpleFlux([prcp, ep] => [pn, en], exprs=[prcp - min(prcp, ep), ep - min(prcp, ep)]),
        SimpleFlux([pn, soilwater] => [ps], [x1], exprs=[max(0.0, pn * (1 - (soilwater / x1)^2))]),
        SimpleFlux([en, soilwater] => [es], [x1], exprs=[en * (2 * soilwater / x1 - (soilwater / x1)^2)]),
        SimpleFlux([soilwater] => [perc], [x1], exprs=[((x1)^(-4)) / 4 * ((4 / 9)^(4)) * (soilwater^5)]),
        SimpleFlux([pn, ps, perc] => [pr], [x1], exprs=[pn - ps + perc]),
        SimpleFlux([pr] => [slowflow, fastflow], exprs=[0.9 * pr, 0.1 * pr]),
        SimpleFlux([ps, es, perc, soilwater] => [new_soilwater], exprs=[soilwater + ps - es - perc])
    ]
    prod_dfuncs = [StateFlux(soilwater => new_soilwater)]

    uh_flux_1 = UnitHydroFlux(slowflow, x4, HydroModels.uh_1_half)
    uh_flux_2 = UnitHydroFlux(fastflow, x4, HydroModels.uh_2_full)

    prod_ele = HydroBucket(:gr4j_prod, funcs=prod_funcs, dfuncs=prod_dfuncs)
    #* define the routing store
    rst_funcs = [
        SimpleFlux([routingstore] => [exch], [x2, x3], exprs=[x2 * abs(routingstore / x3)^3.5]),
        SimpleFlux([routingstore, slowflow_routed, exch] => [routedflow], [x3], exprs=[x3^(-4) / 4 * (routingstore + slowflow_routed + exch)^5]),
        SimpleFlux([routedflow, fastflow_routed, exch] => [flow], exprs=[routedflow + max(fastflow_routed + exch, 0.0)]),
        SimpleFlux([slowflow_routed, exch, routedflow, routingstore] => [new_routingstore], exprs=[routingstore + slowflow_routed + exch - routedflow])
    ]
    rst_dfuncs = [StateFlux(routingstore => new_routingstore)]
    rst_ele = HydroBucket(:gr4j_rst, funcs=rst_funcs, dfuncs=rst_dfuncs)
    #* define the gr4j model
    model = HydroModel(:gr4j, components=[prod_ele, uh_flux_1, uh_flux_2, rst_ele])

    @test Set(HydroModels.get_input_names(model)) == Set([:prcp, :ep])
    @test Set(HydroModels.get_param_names(model)) == Set([:x1, :x2, :x3, :x4])
    @test Set(HydroModels.get_state_names(model)) == Set([:soilwater, :routingstore])
    @test Set(HydroModels.get_output_names(model)) == Set([:en, :routedflow, :pr, :exch, :pn, :fastflow, :ps, :flow, :slowflow_routed, :perc,
        :new_soilwater, :es, :new_routingstore, :slowflow, :fastflow_routed])
    @test Set(HydroModels.get_var_names(model)) == Set([:prcp, :ep, :soilwater, :new_soilwater, :pn, :en, :ps, :es, :perc, :pr, :slowflow,
        :fastflow, :slowflow_routed, :fastflow_routed, :exch, :routedflow, :flow, :new_routingstore, :routingstore])

    # Test single-node model run
    result_ntp = model(input_ntp, pas, timeidx=ts)
    @test result_ntp isa NamedTuple

    result_mat = model(input_ntp, pas, timeidx=ts, output_type=:matrix)
    @test size(result_mat) == (length(HydroModels.get_var_names(model)), length(ts))

    # Test multi-node model run
    input_ntp_vec = repeat([input_ntp], 10)
    node_names = [Symbol(:node_, i) for i in 1:10]
    node_params = NamedTuple{Tuple(node_names)}(repeat([params], 10))
    node_initstates = NamedTuple{Tuple(node_names)}(repeat([initstates], 10))
    node_pas = ComponentVector(params=node_params, initstates=node_initstates)

    # Test output as vector of NamedTuples
    result_ntp_vec = model(input_ntp_vec, node_pas, timeidx=ts, output_type=:namedtuple)
    @test result_ntp_vec isa Vector
    @test result_ntp_vec[1] isa NamedTuple

    # Test output as 3D array
    result_mat_vec = model(input_ntp_vec, node_pas, timeidx=ts, output_type=:array)
    @test size(result_mat_vec) == (length(HydroModels.get_var_names(model)), length(input_ntp_vec), length(ts))
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
    df = DataFrame(CSV.File("data/m50/01013500.csv"))
    ts = collect(1:10000)
    prcp_vec = df[ts, "Prcp"]
    temp_vec = df[ts, "Temp"]
    dayl_vec = df[ts, "Lday"]
    snowpack_vec = df[ts, "SnowWater"]
    soilwater_vec = df[ts, "SoilWater"]
    qobs_vec = df[ts, "Flow"]

    inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
    means, stds = mean.(inputs), std.(inputs)
    (prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec) = [
        @.((tmp_vec - mean) / std) for (tmp_vec, mean, std) in zip(inputs, means, stds)
    ]

    #! define the snow pack reservoir
    snow_funcs = [
        SimpleFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        SimpleFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        SimpleFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    snow_dfuncs = [StateFlux([snowfall] => [melt], snowpack)]
    snow_ele = HydroBucket(:exphydro_snow, funcs=snow_funcs, dfuncs=snow_dfuncs)

    #! define the ET NN and Q NN
    et_nn = Lux.Chain(Lux.Dense(3 => 16, Lux.tanh), Lux.Dense(16 => 16, Lux.leakyrelu), Lux.Dense(16 => 1, Lux.leakyrelu), name=:etnn)
    et_nn_p = Vector(ComponentVector(LuxCore.initialparameters(StableRNG(42), et_nn)))
    q_nn = Lux.Chain(Lux.Dense(2 => 16, Lux.tanh), Lux.Dense(16 => 16, Lux.leakyrelu), Lux.Dense(16 => 1, Lux.leakyrelu), name=:qnn)
    q_nn_p = Vector(ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn)))

    #! get init parameters for each NN
    et_nn_flux = NeuralFlux([norm_snw, norm_slw, norm_temp] => [log_evap_div_lday], et_nn)
    q_nn_flux = NeuralFlux([norm_slw, norm_prcp] => [log_flow], q_nn)

    #! define the soil water reservoir
    soil_funcs = [
        #* normalize
        StdMeanNormFlux(
            [snowpack, soilwater, prcp, temp] => [norm_snw, norm_slw, norm_prcp, norm_temp],
            [[snowpack_mean, snowpack_std], [soilwater_mean, soilwater_std], [prcp_mean, prcp_std], [temp_mean, temp_std]]
        ),
        et_nn_flux,
        q_nn_flux,
    ]

    state_expr = rainfall + melt - step_func(soilwater) * lday * log_evap_div_lday - step_func(soilwater) * exp(log_flow)
    soil_dfuncs = [StateFlux([soilwater, rainfall, melt, lday, log_evap_div_lday, log_flow], soilwater, Num[], expr=state_expr)]
    soil_ele = HydroBucket(:m50_soil, funcs=soil_funcs, dfuncs=soil_dfuncs)

    #! define the Exp-Hydro model
    model = HydroModel(:m50, components=[snow_ele, soil_ele])

    @test Set(HydroModels.get_input_names(model)) == Set([:prcp, :temp, :lday])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :snowpack_std, :snowpack_mean, :soilwater_std, :soilwater_mean, :prcp_std, :prcp_mean, :temp_std, :temp_mean])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater])
    @test Set(HydroModels.get_nn_names(model)) == Set([:etnn, :qnn])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :rainfall, :snowfall, :melt, :log_evap_div_lday, :log_flow, :norm_snw, :norm_slw, :norm_temp, :norm_prcp])
    @test Set(HydroModels.get_var_names(model)) == Set([:prcp, :temp, :lday, :pet, :rainfall, :snowfall, :snowpack, :soilwater, :melt,
        :log_evap_div_lday, :log_flow, :norm_snw, :norm_slw, :norm_temp, :norm_prcp])

    base_params = (Df=2.674, Tmax=0.17, Tmin=-2.09)
    var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
    var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
    nn_params = (etnn=et_nn_p, qnn=q_nn_p)
    params = reduce(merge, [base_params, var_means, var_stds])
    initstates = (snowpack=0.0, soilwater=1303.00)
    pas = ComponentVector(initstates=initstates, params=params, nn=nn_params)
    input_ntp = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
    # Run the model for a single set of inputs and parameters
    result_ntp = model(input_ntp, pas, timeidx=ts)
    @test result_ntp isa NamedTuple

    # Run the model and get results as a matrix
    result_mat = model(input_ntp, pas, timeidx=ts, output_type=:matrix)
    @test size(result_mat) == (length(HydroModels.get_var_names(model)), length(ts))

    # Prepare inputs and parameters for multiple nodes
    input_ntp_vec = repeat([input_ntp], 10)  # 10 identical input sets
    node_names = [Symbol(:node_, i) for i in 1:10]
    node_params = NamedTuple{Tuple(node_names)}(repeat([params], 10))
    node_initstates = NamedTuple{Tuple(node_names)}(repeat([initstates], 10))
    node_pas = ComponentVector(params=node_params, initstates=node_initstates, nn=nn_params)

    # Run the model for multiple nodes and get results as a vector of NamedTuples
    result_ntp_vec = model(input_ntp_vec, node_pas, timeidx=ts, output_type=:namedtuple)
    @test result_ntp_vec isa Vector
    @test result_ntp_vec[1] isa NamedTuple

    # Run the model for multiple nodes and get results as a 3D array
    result_mat_vec = model(input_ntp_vec, node_pas, timeidx=ts, output_type=:array)
    @test size(result_mat_vec) == (length(HydroModels.get_var_names(model)), length(input_ntp_vec), length(ts))
end