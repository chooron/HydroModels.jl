@testset "test ode solved results" begin
    file_path = "../data/exphydro/01013500.csv"
    data = CSV.File(file_path)
    df = DataFrame(data)
    ts = collect(1:10000)
    lday_vec = df[ts, "dayl(day)"]
    prcp_vec = df[ts, "prcp(mm/day)"]
    temp_vec = df[ts, "tmean(C)"]
    flow_vec = df[ts, "flow(mm)"]
    input_ntp = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
    init_states = ComponentVector(snowpack=0.0)
    params = ComponentVector(Df=2.674548848, Tmax=0.175739196, Tmin=-2.092959084)
    pas = ComponentVector(params=params, initstates=init_states)
    prcp_itp = LinearInterpolation(prcp_vec, ts)
    temp_itp = LinearInterpolation(temp_vec, ts)
    input = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

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

    function snowpack_bucket!(du, u, p, t)
        snowpack_ = u[1]
        Df, Tmax, Tmin = p.Df, p.Tmax, p.Tmin
        prcp_, temp_ = prcp_itp(t), temp_itp(t)
        snowfall_ = step_func(Tmin - temp_) * prcp_
        melt_ = step_func(temp_ - Tmax) * step_func(snowpack_) * min(snowpack_, Df * (temp_ - Tmax))
        du[1] = snowfall_ - melt_
    end
    prob = ODEProblem(snowpack_bucket!, [init_states.snowpack], (ts[1], ts[end]), params)
    sol = solve(prob, Tsit5(), saveat=ts, reltol=1e-3, abstol=1e-3)
    num_u = length(prob.u0)
    manual_result = [sol[i, :] for i in 1:num_u]
    ele_params_idx = [getaxes(pas[:params])[1][nm].idx for nm in HydroModels.get_param_names(snow_ele)]
    paramfunc = (p) -> [p[:params][idx] for idx in ele_params_idx]
    param_func, nn_param_func = HydroModels._get_parameter_extractors(snow_ele, pas)
    itpfunc_list = map((var) -> LinearInterpolation(var, ts, extrapolate=true), eachrow(input))
    ode_input_func = (t) -> [itpfunc(t) for itpfunc in itpfunc_list]
    du_func = HydroModels._get_du_func(snow_ele, ode_input_func, param_func, nn_param_func)
    solver = ODESolver(alg=Tsit5(), reltol=1e-3, abstol=1e-3)
    initstates_mat = collect(pas[:initstates][HydroModels.get_state_names(snow_ele)])
    #* solve the problem by call the solver
    solved_states = solver(du_func, pas, initstates_mat, ts)
    @test manual_result[1] == solved_states[1, :]
end