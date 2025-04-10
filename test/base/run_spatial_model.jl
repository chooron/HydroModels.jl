step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

@parameters Tmin Tmax Df Smax f Qmax area_coef lag
@variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow q q_gen q_routed s_river

snow_ele = @hydrobucket :surface begin
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

soil_ele = @hydrobucket :soil begin
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

convertflux = @hydroflux q ~ flow * area_coef

@testset "test grid route hydro model (multiple hydrology nodes based on exp-hydro)" begin

    #! load data
    ts = collect(1:100)
    df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input_mat = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00, s_river=0.0)
    params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09, lag=0.2, area_coef=1.0)
    pas = ComponentVector(params=params)

    #! define the routing method
    flwdir = [1 4 8; 1 4 4; 1 1 2]
    positions = [[1, 1], [1, 2], [1, 3], [2, 1], [2, 2], [2, 3], [3, 1], [3, 2], [3, 3]]
    route = @hydroroute begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag) + q
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(flwdir, positions)
    end
    #! define the Exp-Hydro model
    model = @hydromodel :exphydro begin
        snow_ele
        soil_ele
        convertflux
        route
    end

    @test Set(HydroModels.get_input_names(model)) == Set([:temp, :lday, :prcp])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax, :lag, :area_coef])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater, :s_river])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :q, :q_routed])
    @test Set(reduce(union, HydroModels.get_var_names(model))) == Set([:temp, :lday, :prcp, :pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :snowpack, :soilwater, :s_river, :q, :q_routed])

    input_arr = repeat(reshape(input_mat, size(input_mat)[1], 1, size(input_mat)[2]), 1, 9, 1)
    node_params = ComponentVector(
        f=fill(0.0167, 9), Smax=fill(1709.46, 9), Qmax=fill(18.47, 9), Df=fill(2.674, 9),
        Tmax=fill(0.17, 9), Tmin=fill(-2.09, 9), lag=fill(0.2, 9), area_coef=fill(1.0, 9)
    )
    node_initstates = ComponentVector(
        snowpack=fill(0.0, 9), soilwater=fill(0.0, 9), s_river=fill(0.0, 9)
    )
    node_pas = ComponentVector(params=node_params)

    config = (timeidx=ts, ptyidx=1:9, styidx=1:9)
    result_mat_vec = model(input_arr, node_pas; initstates=node_initstates, config=config)
    @test size(result_mat_vec) == (length(HydroModels.get_state_names(model)) + length(HydroModels.get_output_names(model)), 9, length(ts))
end

@testset "test vector route hydro model (spatial hydrology model based on vector route and exp-hydro)" begin
    @parameters Tmin Tmax Df Smax f Qmax
    @variables prcp temp lday pet snowpack soilwater rainfall snowfall evap melt baseflow surfaceflow flow flow_routed

    #! load data
    ts = collect(1:100)
    df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input_mat = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')

    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00, s_river=0.0)
    params = ComponentVector(f=0.0167, Smax=1709.46, Qmax=18.47, Df=2.674, Tmax=0.17, Tmin=-2.09, lag=0.2, k=0.5, x=0.2, area_coef=1.0)
    pas = ComponentVector(params=params)

    #! define the routing method
    network = DiGraph(9)
    add_edge!(network, 1, 2)
    add_edge!(network, 2, 5)
    add_edge!(network, 3, 5)
    add_edge!(network, 4, 5)
    add_edge!(network, 5, 8)
    add_edge!(network, 6, 9)
    add_edge!(network, 7, 8)
    add_edge!(network, 8, 9)

    @variables q q_routed s_river q_gen
    @parameters area_coef lag

    discharge_route = @hydroroute :exphydro_routed begin
        fluxes = begin
            @hydroflux q_routed ~ s_river / (1 + lag) + q
        end
        dfluxes = begin
            @stateflux s_river ~ q - q_routed
        end
        aggr_func = HydroModels.build_aggr_func(network)
    end

    #! define the Exp-Hydro model
    model = @hydromodel :exphydro begin
        snow_ele
        soil_ele
        convertflux
        discharge_route
    end

    @test Set(HydroModels.get_input_names(model)) == Set([:temp, :lday, :prcp])
    @test Set(HydroModels.get_param_names(model)) == Set([:Tmin, :Tmax, :Df, :Smax, :f, :Qmax, :lag, :area_coef])
    @test Set(HydroModels.get_state_names(model)) == Set([:snowpack, :soilwater, :s_river])
    @test Set(HydroModels.get_output_names(model)) == Set([:pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :flow, :q, :q_routed])
    @test Set(reduce(union, HydroModels.get_var_names(model))) == Set(
        [:temp, :lday, :prcp, :pet, :snowfall, :rainfall, :melt, :evap, :baseflow, :surfaceflow, :s_river, :flow, :snowpack, :soilwater, :q, :q_routed]
    )

    input_arr = repeat(reshape(input_mat, size(input_mat)[1], 1, size(input_mat)[2]), 1, 9, 1)
    node_params = ComponentVector(
        f=fill(0.0167, 9), Smax=fill(1709.46, 9), Qmax=fill(18.47, 9), Df=fill(2.674, 9),
        Tmax=fill(0.17, 9), Tmin=fill(-2.09, 9), lag=fill(0.2, 9), area_coef=fill(1.0, 9)
    )
    node_initstates = ComponentVector(
        snowpack=fill(0.0, 9), soilwater=fill(0.0, 9), s_river=fill(0.0, 9)
    )
    node_pas = ComponentVector(params=node_params)

    config = (timeidx=ts,)
    result_mat_vec = model(input_arr, node_pas; initstates=node_initstates, config=config)
    @test size(result_mat_vec) == (length(HydroModels.get_state_names(model)) + length(HydroModels.get_output_names(model)), 9, length(ts))
end
