using Test
using HydroModels
using SciMLBase
using OrdinaryDiffEq
using OrdinaryDiffEqFunctionMap
using DiffEqCallbacks
using SciMLSensitivity

Base.retry_load_extensions()

@testset "SciML extension boundaries" begin
    ode_ext = Base.get_extension(HydroModels, :HydroModelsOrdinaryDiffEqExt)
    fmap_ext = Base.get_extension(HydroModels, :HydroModelsOrdinaryDiffEqFunctionMapExt)
    callbacks_ext = Base.get_extension(HydroModels, :HydroModelsDiffEqCallbacksExt)

    @test ode_ext !== nothing
    @test fmap_ext !== nothing
    @test callbacks_ext !== nothing
    @test !isdefined(ode_ext, :SciMLSensitivity)

    du_func(u, p, t) = p .* u
    params = [0.5]
    initstates = [1.0]
    timeidx = collect(0.0:1.0:3.0)

    ode_output = hydrosolve(
        HydroModels.ODESolver,
        du_func,
        params,
        initstates,
        timeidx,
        (interpolator=HydroModels.ConstantInterpolation,),
    )
    @test size(ode_output) == (1, 4)
    @test eltype(ode_output) == Float64
    @test ode_output ≈ [1.0 1.6487212698173226 2.71828175675071 4.481688763572211]

    ode_sense_output = hydrosolve(
        HydroModels.ODESolver,
        du_func,
        params,
        initstates,
        timeidx,
        (interpolator=HydroModels.ConstantInterpolation, sense_alg=ForwardDiffSensitivity()),
    )
    @test ode_sense_output ≈ ode_output

    for T in (Float32, Float64)
        local input = reshape(T.([1, 2, 3, 4]), 1, :)
        local bucket = HydroModels.HydroBucket(
            (x, u, p) -> [u[1] + x[1]],
            (x, u, p) -> [p[1]],
            name=:callback_test,
            inputs=[:forcing],
            outputs=[:flux],
            states=[:storage],
            params=[:rate],
        )
        local prob, saved_values = SciMLBase.ODEProblem(
            bucket,
            input;
            params=T[0.5],
            timeidx=collect(1:4),
        )
        local sol = solve(prob, Tsit5(); saveat=collect(1:4))

        @test saved_values isa SavedValues{T,Vector{T}}
        @test eltype(saved_values.t) == T
        @test all(value -> value isa Vector{T}, saved_values.saveval)
        @test length(saved_values.saveval) > 0
        @test saved_values.saveval[1] !== saved_values.saveval[end]
        @test size(Array(sol)) == (1, 4)
    end

    discrete_output = hydrosolve(
        HydroModels.DiscreteSolver,
        (u, p, t) -> u .+ p,
        [2.0],
        [0.0],
        collect(0.0:1.0:3.0),
        NamedTuple(),
    )
    @test discrete_output == [2.0 6.0 14.0 30.0]

    @test fmap_ext.FunctionMap === FunctionMap
    @test fmap_ext.FunctionMap{true}() isa FunctionMap{true}
end
