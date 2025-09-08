Df_v, Tmax_v, Tmin_v = 2.674548848, 0.175739196, -2.092959084
params = ComponentVector(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v)
init_states = ComponentVector(snowpack=0.0)

ts = collect(1:10)
df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

@variables temp lday prcp pet snowfall rainfall melt snowpack
@parameters Tmin Tmax Df


@testset "test multi hydro element (with state)" begin
    snow_single_ele = @hydrobucket :snow begin
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

    snow_multi_ele1 = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
        hru_types = collect(1:10)
    end


    snow_multi_ele2 = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
            @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
        end
        dfluxes = begin
            @stateflux snowpack ~ snowfall - melt
        end
        hru_types = [1, 2, 2, 2, 1, 3, 3, 2, 3, 2]
    end

    @testset "test run with multiple nodes input (independent parameters)" begin
        node_num = 10
        node_names = [Symbol(:node, i) for i in 1:node_num]
        node_params = ComponentVector(Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num))
        node_states = ComponentVector(snowpack=fill(0.0, node_num))
        input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_single_ele)]))
        node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
        node_input = permutedims(node_input, (2, 3, 1))
        node_output = snow_multi_ele1(node_input, ComponentVector(params=node_params); initstates=node_states)
        single_output = snow_single_ele(input, ComponentVector(params=(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v)), initstates=init_states, timeidx=ts, solver=ManualSolver(mutable=true))
        target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
        @test node_output == target_output
    end

    @testset "test run with multiple nodes input (share parameters)" begin
        # share parameters
        node_num = 10
        node_names = [Symbol(:node, i) for i in 1:node_num]
        node_params = ComponentVector(Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num))
        node_states = ComponentVector(snowpack=fill(0.0, node_num))

        node_pas = ComponentVector(params=node_params)
        input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_single_ele)]))
        node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], 10))
        node_input = permutedims(node_input, (2, 3, 1))
        node_output = snow_multi_ele2(node_input, node_pas; initstates=node_states)
        single_output = snow_single_ele(input, ComponentVector(params=(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v)), initstates=init_states, timeidx=ts, solver=ManualSolver(mutable=true))
        target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
        @test node_output == target_output
    end
end

@testset "test multi hydro element (without state)" begin
    @variables temp lday prcp pet snowfall rainfall melt snowpack
    @parameters Tmin Tmax Df

    snow_single_ele = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
    end

    snow_multi_ele1 = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
        hru_types = collect(1:10)
    end
    snow_multi_ele2 = @hydrobucket :snow begin
        fluxes = begin
            @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
            @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
            @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        end
        hru_types = [1, 2, 3, 3, 2, 1, 2, 1, 3, 2]
    end

    pas = ComponentVector(params=(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v))
    @testset "test run with multiple nodes input (independent parameters)" begin
        node_num = 10
        node_names = [Symbol(:node, i) for i in 1:node_num]
        node_params = ComponentVector(
            Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num)
        )
        node_states = ComponentVector(snowpack=fill(0.0, node_num))

        node_pas = ComponentVector(params=node_params[HydroModels.get_param_names(snow_single_ele)], initstates=node_states[HydroModels.get_state_names(snow_single_ele)])
        input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_single_ele)]))
        node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
        node_input = permutedims(node_input, (2, 3, 1))
        node_output = snow_multi_ele1(node_input, node_pas; initstates=node_states)
        single_output = snow_single_ele(input, pas, initstates=init_states, timeidx=ts, solver=ManualSolver(mutable=true))
        target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
        @test node_output == target_output
    end

    @testset "test run with multiple nodes input (share parameters)" begin
        # share parameters
        node_num = 3
        node_names = [Symbol(:node, i) for i in 1:node_num]
        node_params = ComponentVector(Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num))
        node_states = ComponentVector(snowpack=fill(0.0, node_num))

        node_pas = ComponentVector(params=node_params[HydroModels.get_param_names(snow_single_ele)])
        input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_single_ele)]))
        node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], 10))
        node_input = permutedims(node_input, (2, 3, 1))
        node_output = snow_multi_ele2(node_input, node_pas; initstates=node_states)
        single_output = snow_single_ele(input, pas, initstates=init_states, timeidx=ts, solver=ManualSolver(mutable=true))
        target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
        @test node_output == target_output
    end
end