@testset "test compact with the DifferentialEquations.jl solvers (in Exphydro snow bucket)" begin
    @variables temp lday prcp pet snowfall rainfall melt snowpack
    @parameters Tmin Tmax Df

    Df_v, Tmax_v, Tmin_v = 2.674548848, 0.175739196, -2.092959084
    params = ComponentVector(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v)
    init_states = ComponentVector(snowpack=0.0)

    ts = collect(1:100)
    df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
    input_ntp = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])

    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

    snow_ele = @hydrobucket :snow begin
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

    single_input = Matrix(reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_ele)]))')
    single_pas = ComponentVector(params=(Df=Df_v, Tmax=Tmax_v, Tmin=Tmin_v))

    node_num = 10
    node_params = ComponentVector(Df=fill(Df_v, node_num), Tmax=fill(Tmax_v, node_num), Tmin=fill(Tmin_v, node_num))
    node_pas = ComponentVector(params=node_params)
    node_states = ComponentVector(snowpack=fill(0.0, node_num))
    input_arr = reduce(hcat, collect(input_ntp[HydroModels.get_input_names(snow_ele)]))
    node_input = reduce((m1, m2) -> cat(m1, m2, dims=3), repeat([input_arr], node_num))
    node_input = permutedims(node_input, (2, 3, 1))

    @testset "test the ODEProblem and DiscreteProblem Solvers in single node" begin
        ode_result_1 = snow_ele(single_input, single_pas; initstates=init_states, timeidx=ts, solver=HydroModelSolvers.ODESolver())
        ode_result_2 = snow_ele(
            single_input, single_pas;
            initstates=init_states, timeidx=ts,
            solver=HydroModelSolvers.ODESolver(alg=BS3(), sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()), kwargs=Dict(:reltol => 1e-6, :abstol => 1e-6))
        )
        disc_result = snow_ele(single_input, single_pas; initstates=init_states, timeidx=ts, solver=HydroModelSolvers.DiscreteSolver())
        @test size(ode_result_1) == size(ode_result_2) == size(disc_result)
    end

    @testset "test the ODEProblem and DiscreteProblem Solvers in multiple node" begin
        ode_result_1 = snow_ele(
            node_input, node_pas; initstates=node_states,
            ptyidx=1:10, styidx=1:10, timeidx=ts, solver=HydroModelSolvers.ODESolver()
        )
        ode_result_2 = snow_ele(
            node_input, node_pas;
            initstates=node_states, timeidx=ts, ptyidx=1:10, styidx=1:10,
            solver=HydroModelSolvers.ODESolver(alg=BS3(), sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()), kwargs=Dict(:reltol => 1e-6, :abstol => 1e-6))
        )
        disc_result = snow_ele(node_input, node_pas;
            initstates=node_states, timeidx=ts, ptyidx=1:10, styidx=1:10, solver=HydroModelSolvers.DiscreteSolver()
        )
        @test size(ode_result_1) == size(ode_result_2) == size(disc_result)
    end
end

@testset "test compact with the DifferentialEquations.jl solvers (in Exphydro snow bucket)" begin
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

    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

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

    base_params = (Df=2.674, Tmax=0.17, Tmin=-2.09)
    var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
    var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
    init_states = ComponentVector(snowpack=0.0, soilwater=1303.00)
    single_pas = ComponentVector(params=reduce(merge, [base_params, var_means, var_stds]), nns=(etnn=et_nn_p, qnn=q_nn_p))
    input_ntp = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
    single_input = Matrix(reduce(hcat, collect(input_ntp[HydroModels.get_input_names(model)]))')

    # Prepare inputs and parameters for multiple nodes
    node_input = repeat(reshape(single_input, size(single_input)[1], 1, size(single_input)[2]), 1, 10, 1)
    node_params = ComponentVector(Df=fill(2.674, 10), Tmax=fill(0.17, 10), Tmin=fill(-2.09, 10),
        snowpack_std=fill(var_stds[:snowpack_std], 10), snowpack_mean=fill(var_means[:snowpack_mean], 10),
        soilwater_std=fill(var_stds[:soilwater_std], 10), soilwater_mean=fill(var_means[:soilwater_mean], 10),
        prcp_std=fill(var_stds[:prcp_std], 10), prcp_mean=fill(var_means[:prcp_mean], 10),
        temp_std=fill(var_stds[:temp_std], 10), temp_mean=fill(var_means[:temp_mean], 10))
    node_initstates = ComponentVector(snowpack=fill(0.0, 10), soilwater=fill(1303.00, 10))
    node_pas = ComponentVector(params=node_params,  nns=(etnn=et_nn_p, qnn=q_nn_p))

    @testset "test the ODEProblem and DiscreteProblem Solvers in single node" begin
        ode_result_1 = model(single_input, single_pas; initstates=init_states, config=(timeidx=ts, solver=HydroModelSolvers.ODESolver()))
        ode_result_2 = model(
            single_input, single_pas;
            initstates=init_states,
            config=(timeidx=ts, solver=HydroModelSolvers.ODESolver(alg=BS3(),
                sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()), kwargs=Dict(:reltol => 1e-6, :abstol => 1e-6)))
        )
        disc_result = model(single_input, single_pas; initstates=init_states, timeidx=ts, solver=HydroModelSolvers.DiscreteSolver())
        @test size(ode_result_1) == size(ode_result_2) == size(disc_result)
    end

    @testset "test the ODEProblem and DiscreteProblem Solvers in multiple node" begin
        ode_result_1 = model(
            node_input, node_pas; initstates=node_initstates,
            config=(ptyidx=1:10, styidx=1:10, timeidx=ts, solver=HydroModelSolvers.ODESolver())
        )
        ode_result_2 = model(
            node_input, node_pas;
            initstates=node_initstates,
            config=(timeidx=ts, ptyidx=1:10, styidx=1:10,
                solver=HydroModelSolvers.ODESolver(alg=BS3(), sensealg=BacksolveAdjoint(autojacvec=ZygoteVJP()), kwargs=Dict(:reltol => 1e-6, :abstol => 1e-6))
            )
        )
        disc_result = model(node_input, node_pas;
            initstates=node_initstates, config=(timeidx=ts, ptyidx=1:10, styidx=1:10, solver=HydroModelSolvers.DiscreteSolver())
        )
        @test size(ode_result_1) == size(ode_result_2) == size(disc_result)
    end

end