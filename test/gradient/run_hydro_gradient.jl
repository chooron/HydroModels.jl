@testset "test gradient computation for hydrological model" begin
    include("../models/exphydro.jl")

    f, Smax, Qmax, Df, Tmax, Tmin = 0.01674478, 1709.461015, 18.46996175, 2.674548848, 0.175739196, -2.092959084
    # 0.0167, 1709.46, 18.47, 2.6745, 0.1757, -2.093
    params = ComponentVector(f=f, Smax=Smax, Qmax=Qmax, Df=Df, Tmax=Tmax, Tmin=Tmin)
    init_states = ComponentVector(snowpack=10.0, soilwater=1303.004248)
    pas = ComponentVector(params=params)
    file_path = "../data/exphydro/01013500.csv"
    df = CSV.File(file_path) |> DataFrame
    ts = collect(1:100)

    # single node input
    input = (lday=df[ts, "dayl(day)"], temp=df[ts, "tmean(C)"], prcp=df[ts, "prcp(mm/day)"])
    input_arr = Matrix(reduce(hcat, collect(input[HydroModels.get_input_names(bucket_1)]))')
    config = (solver=HydroModels.ManualSolver(mutable=false), timeidx=ts)

    @testset "test hydro gradient computation under manual solver" begin
        # Test Zygote gradient calculation
        zygote_grad = Zygote.gradient(pas) do p
            bucket_1(input_arr, p, inistates=init_states, solver=HydroModels.ManualSolver(mutable=false)) |> sum
        end

        # Test that gradient is not nothing and contains values
        @test !isnothing(zygote_grad)

        # Test ForwardDiff gradient calculation
        forwarddiff_grad = ForwardDiff.gradient(pas) do p
            bucket_1(input_arr, p, inistates=init_states, solver=HydroModels.ManualSolver(mutable=false)) |> sum
        end

        # Test that gradient is not nothing and contains values
        @test !isnothing(forwarddiff_grad)
    end

    @testset "test hydro gradient computation under ode solver" begin
        # Test Zygote gradient calculation
        zygote_grad = Zygote.gradient(pas) do p
            bucket_1(input_arr, p, inistates=init_states, interp=LinearInterpolation, solver=HydroModelSolvers.ODESolver()) |> sum
        end

        # Test that gradient is not nothing and contains values
        @test !isnothing(zygote_grad)

        # Test ForwardDiff gradient calculation
        forwarddiff_grad = ForwardDiff.gradient(pas) do p
            bucket_1(input_arr, p, inistates=init_states, interp=LinearInterpolation, solver=HydroModelSolvers.ODESolver()) |> sum
        end

        # Test that gradient is not nothing and contains values
        @test !isnothing(forwarddiff_grad)
    end

    @testset "test hydro gradient computation under discrete solver" begin
        # Test Zygote gradient calculation
        zygote_grad = Zygote.gradient(pas) do p
            bucket_1(input_arr, p, inistates=init_states, interp=LinearInterpolation, solver=HydroModelSolvers.DiscreteSolver()) |> sum
        end

        # Test that gradient is not nothing and contains values
        @test !isnothing(zygote_grad)

        # Test ForwardDiff gradient calculation
        forwarddiff_grad = ForwardDiff.gradient(pas) do p
            bucket_1(input_arr, p, inistates=init_states, interp=LinearInterpolation, solver=HydroModelSolvers.DiscreteSolver()) |> sum
        end

        # Test that gradient is not nothing and contains values
        @test !isnothing(forwarddiff_grad)
    end
end