@testset "optimize exphydro" begin
    step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5
    # define variables and parameters
    @variables temp lday pet prcp snowfall rainfall snowpack melt
    @parameters Tmin Tmax Df Smax Qmax f

    @variables soilwater pet evap baseflow surfaceflow flow rainfall

    # define model components
    fluxes_1 = [
        HydroFlux([temp, lday] => [pet], exprs=[29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)]),
        HydroFlux([prcp, temp] => [snowfall, rainfall], [Tmin], exprs=[step_func(Tmin - temp) * prcp, step_func(temp - Tmin) * prcp]),
        HydroFlux([snowpack, temp] => [melt], [Tmax, Df], exprs=[step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))]),
    ]
    dfluxes_1 = [StateFlux([snowfall] => [melt], snowpack),]
    bucket_1 = HydroBucket(name=:surface, fluxes=fluxes_1, dfluxes=dfluxes_1)

    fluxes_2 = [
        HydroFlux([soilwater, pet] => [evap], [Smax], exprs=[step_func(soilwater) * pet * min(1.0, soilwater / Smax)]),
        HydroFlux([soilwater] => [baseflow], [Smax, Qmax, f], exprs=[step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))]),
        HydroFlux([soilwater] => [surfaceflow], [Smax], exprs=[max(0.0, soilwater - Smax)]),
        HydroFlux([baseflow, surfaceflow] => [flow], exprs=[baseflow + surfaceflow]),
    ]
    dfluxes_2 = [StateFlux([rainfall, melt] => [evap, flow], soilwater)]
    bucket_2 = HydroBucket(name=:soil, fluxes=fluxes_2, dfluxes=dfluxes_2)
    model = HydroModel(name=:exphydro, components=[bucket_1, bucket_2])
    ioadapter = NamedTupleIOAdapter(model)

    # predefine the parameters
    f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084

    # load data
    file_path = "../data/exphydro/01013500.csv"
    data = CSV.File(file_path)
    df = DataFrame(data)
    ts = collect(1:10000)
    lday_vec = df[ts, "dayl(day)"]
    prcp_vec = df[ts, "prcp(mm/day)"]
    temp_vec = df[ts, "tmean(C)"]
    flow_vec = df[ts, "flow(mm)"]

    tunable_pas = ComponentVector(params=ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin))
    const_pas = ComponentVector(initstates=ComponentVector(snowpack=0.0, soilwater=1300.0))

    # parameters optimization
    input = (prcp=prcp_vec, lday=lday_vec, temp=temp_vec)
    input_matrix = Matrix(reduce(hcat, collect(input))')
    output = (flow=flow_vec,)
    config = (solver=ODESolver(), interp=LinearInterpolation)

    @testset "HydroOptimizer" begin
        # build optimizer
        hydro_opt = HydroOptimizer(component=ioadapter, maxiters=100)
        lb_list = [0.0, 100.0, 10.0, 0.0, 0.0, -3.0]
        ub_list = [0.1, 2000.0, 50.0, 5.0, 3.0, 0.0]
        config = (solver=ODESolver(), interp=LinearInterpolation)
        opt_params, loss_df = hydro_opt([input], [output], tunable_pas=tunable_pas, const_pas=const_pas, config=[config], lb=lb_list, ub=ub_list, return_loss_df=true)
    end

    @testset "GradOptimizer" begin
        config = (solver=ODESolver(sensealg=InterpolatingAdjoint(autodiff=true)), interp=LinearInterpolation)
        grad_opt = GradOptimizer(component=ioadapter, maxiters=100)
        opt_params, loss_df = grad_opt([input], [output], tunable_pas=tunable_pas, const_pas=const_pas, config=[config], return_loss_df=true)
    end

    # @testset "BatchOptimizer" begin
    #     config = (solver=ODESolver(sensealg=InterpolatingAdjoint(autodiff=ZygoteVJP())), interp=LinearInterpolation)
    #     batch_opt = BatchOptimizer(component=ioadapter, maxiters=100, solve_alg=Adam(1e-2), adtype=Optimization.AutoZygote())
    #     opt_pas, loss_df = batch_opt(
    #         fill(input, 5), fill(output, 5),
    #         tunable_pas=tunable_pas,
    #         const_pas=const_pas, config=fill(config, 5),
    #         return_loss_df=true
    #     )
    # end
end