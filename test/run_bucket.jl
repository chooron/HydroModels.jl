@testset "test hydrobucket" begin
    Df_v, Tmax_v, Tmin_v = 2.674548848, 0.175739196, -2.092959084
    params = ComponentVector(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v)
    init_states = ComponentVector(snowpack=0.0)

    ts = collect(1:10)
    df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

    input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')
    dtype = eltype(input[1])
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

    @testset "test hydro element (basic element, Snowpack in Exp-Hydro)" begin
        @variables temp lday prcp pet snowfall rainfall melt snowpack
        @parameters Tmin Tmax Df

        snow_fluxes = [
            HydroModels.HydroFlux([temp, lday] => [pet],
                exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
            HydroModels.HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin],
                exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
            HydroModels.HydroFlux([snowpack, temp] => [melt], [Tmax, Df],
                exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
        ]
        snow_dfluxes = [HydroModels.StateFlux([snowfall] => [melt], snowpack)]
        snow_ele = HydroModels.HydroBucket(fluxes=snow_fluxes, dfluxes=snow_dfluxes)
        @testset "test hydro element info" begin
            @test Set(HydroModels.get_input_names(snow_ele)) == Set((:temp, :lday, :prcp))
            @test Set(HydroModels.get_param_names(snow_ele)) == Set((:Tmin, :Tmax, :Df))
            @test Set(HydroModels.get_output_names(snow_ele)) == Set((:pet, :snowfall, :rainfall, :melt))
            @test Set(HydroModels.get_state_names(snow_ele)) == Set((:snowpack,))
        end
        pas = ComponentVector(params=params[HydroModels.get_param_names(snow_ele)])
        result = snow_ele(input, pas; initstates=init_states, timeidx=ts, solver=ManualSolver{true}())
        ele_state_and_output_names = vcat(HydroModels.get_state_names(snow_ele), HydroModels.get_output_names(snow_ele))
        result = NamedTuple{Tuple(ele_state_and_output_names)}(eachslice(result, dims=1))

        # @testset "test first output for hydro element" begin
        #     snowpack0 = init_states[:snowpack]
        #     pet0 = snow_fluxes[1]([input_ntp.temp[1], input_ntp.lday[1]], ComponentVector(params=ComponentVector()))[1]
        #     snowfall0, rainfall0 = snow_fluxes[2]([input_ntp.prcp[1], input_ntp.temp[1]], ComponentVector(params=(Tmin=params.Tmin,)))
        #     melt0 = snow_fluxes[3]([snowpack0, input_ntp.temp[1]], ComponentVector(params=(Tmax=params.Tmax, Df=params.Df)))[1]
        #     @test snowpack0 == result.snowpack[1]
        #     @test snowfall0 == result.snowfall[1]
        #     @test rainfall0 == result.rainfall[1]
        #     @test melt0 == result.melt[1]
        # end

        # @testset "test all of the output" begin
        #     param_func, nn_param_func = HydroModels._get_parameter_extractors(snow_ele, pas)
        #     itpfunc_list = map((var) -> LinearInterpolation(var, ts, extrapolate=true), eachrow(input))
        #     ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]
        #     du_func = HydroModels._get_du_func(snow_ele, ode_input_func, param_func, nn_param_func)
        #     solver = ManualSolver{true}()
        #     initstates_mat = collect(pas[:initstates][HydroModels.get_state_names(snow_ele)])
        #     #* solve the problem by call the solver
        #     snowpack_vec = solver(du_func, pas, initstates_mat, ts)[1, :]
        #     pet_vec = snow_fluxes[1](Matrix(reduce(hcat, [input_ntp.temp, input_ntp.lday])'), ComponentVector(params=ComponentVector()))[1, :]
        #     snow_fluxes_2_output = snow_fluxes[2](Matrix(reduce(hcat, [input_ntp.prcp, input_ntp.temp])'), ComponentVector(params=(Tmin=params.Tmin,)))
        #     snowfall_vec, rainfall_vec = snow_fluxes_2_output[1, :], snow_fluxes_2_output[2, :]
        #     melt_vec = snow_fluxes[3](Matrix(reduce(hcat, [snowpack_vec, input_ntp.temp])'), ComponentVector(params=(Tmax=params.Tmax, Df=params.Df)))[1, :]
        #     @test reduce(vcat, pet_vec) == collect(result.pet)
        #     @test reduce(vcat, snowfall_vec) == collect(result.snowfall)
        #     @test reduce(vcat, rainfall_vec) == collect(result.rainfall)
        #     @test reduce(vcat, melt_vec) == collect(result.melt)
        # end

        @testset "test run with multiple nodes input (independent parameters)" begin
            node_num = 10
            node_names = [Symbol(:node, i) for i in 1:node_num]
            node_params = ComponentVector(
                Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num)
            )
            node_states = ComponentVector(snowpack=fill(0.0, node_num))

            node_pas = ComponentVector(params=node_params[HydroModels.get_param_names(snow_ele)], initstates=node_states[HydroModels.get_state_names(snow_ele)])
            input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_ele)]))
            node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], length(node_names)))
            node_input = permutedims(node_input, (2, 3, 1))
            config = (ptyidx=1:10, styidx=1:10, timeidx=ts)
            node_output = snow_ele(node_input, node_pas; initstates=node_states, config...)
            single_output = snow_ele(input, pas, initstates=init_states, timeidx=ts, solver=ManualSolver{true}())
            target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
            @test node_output == target_output
        end

        @testset "test run with multiple nodes input (share parameters)" begin
            # share parameters
            node_num = 3
            node_names = [Symbol(:node, i) for i in 1:node_num]
            node_params = ComponentVector(Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num))
            node_states = ComponentVector(snowpack=fill(0.0, node_num))

            node_pas = ComponentVector(params=node_params[HydroModels.get_param_names(snow_ele)])
            input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_ele)]))
            node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], 10))
            node_input = permutedims(node_input, (2, 3, 1))
            config = (ptyidx=[1, 2, 2, 2, 1, 3, 3, 2, 3, 2], styidx=[1, 2, 2, 2, 1, 3, 3, 2, 3, 2], timeidx=ts)
            node_output = snow_ele(node_input, node_pas; initstates=node_states, config...)
            single_output = snow_ele(input, pas, initstates=init_states, timeidx=ts, solver=ManualSolver{true}())
            target_output = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([single_output], 10)), (1, 3, 2))
            @test node_output == target_output
        end
    end
end