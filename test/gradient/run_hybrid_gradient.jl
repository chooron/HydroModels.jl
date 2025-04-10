@testset "test m50 gradient computation based on Zygote" begin
    include("../models/m50.jl")

    #! load data
    df = DataFrame(CSV.File("../data/m50/01013500.csv"))
    ts = collect(1:100)
    prcp_vec = df[ts, "Prcp"]
    temp_vec = df[ts, "Temp"]
    dayl_vec = df[ts, "Lday"]
    snowpack_vec = df[ts, "SnowWater"]
    soilwater_vec = df[ts, "SoilWater"]
    qobs_vec = df[ts, "Flow"]

    inputs = [prcp_vec, temp_vec, snowpack_vec, soilwater_vec]
    means, stds = mean.(inputs), std.(inputs)
    (prcp_norm_vec, temp_norm_vec, snowpack_norm_vec, soilwater_norm_vec) = [@.((inp - mean) / std) for (inp, mean, std) in zip(inputs, means, stds)]

    et_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), ep_nn))
    q_nn_p = ComponentVector(LuxCore.initialparameters(StableRNG(42), q_nn))

    base_params = (Df=2.674, Tmax=0.17, Tmin=-2.09)
    var_stds = NamedTuple{Tuple([Symbol(nm, :_std) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(stds)
    var_means = NamedTuple{Tuple([Symbol(nm, :_mean) for nm in [:prcp, :temp, :snowpack, :soilwater]])}(means)
    nn_params = (epnn=et_nn_p, qnn=q_nn_p)
    params = reduce(merge, [base_params, var_means, var_stds])
    initstates = ComponentVector(snowpack=0.0, soilwater=1303.00)
    pas = ComponentVector(params=params, nns=nn_params)
    input_ntp = (prcp=prcp_vec, lday=dayl_vec, temp=temp_vec)
    input_mat = Matrix(reduce(hcat, collect(input_ntp[HydroModels.get_input_names(m50_model)]))')

    # Run the model to get output
    output = m50_model(input_mat, pas, initstates=initstates,
        config=(timeidx=ts, solver=HydroModelSolvers.ODESolver(sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())))
    )

    # Calculate gradient using Zygote
    gradient_result = Zygote.gradient(pas) do p
        output = m50_model(input_mat, p, initstates=initstates,
            config=(timeidx=ts, solver=HydroModelSolvers.ODESolver(sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())))
        )
        output[end, :] |> sum
    end

    # Test that gradient calculation runs successfully
    @test !isnothing(gradient_result)
    @test !isnothing(gradient_result[1])
end