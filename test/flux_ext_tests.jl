using ComponentArrays
using Flux
using OrdinaryDiffEq
using OrdinaryDiffEqFunctionMap
using SciMLBase
using HydroModels
using StableRNGs
using Test

@testset "Flux extension" begin
    @test !isnothing(Base.get_extension(HydroModels, :HydroModelsFluxExt))

    @variables flux_x flux_y flux_out
    model = Flux.Chain(Flux.Dense(2 => 3, tanh), Flux.Dense(3 => 1))
    flux = HydroModels.NeuralFlux(
        [flux_x, flux_y], [flux_out], model; chain_name=:flux_net, norm=x -> x
    )
    @test HydroModels.get_nn_names(flux) == [:flux_net]

    params = HydroModels.get_initial_params(flux; rng=StableRNG(11))
    flat_params = params.nns.flux_net.params
    @test length(flat_params) == length(first(Flux.destructure(model)))
    @test all(isfinite, flat_params)

    x = Float32[0.2 0.4 -0.1; 1.0 -0.5 0.7]
    expected = flux.chain.rebuild(flat_params)(x)
    @test flux(x, params) ≈ expected

    saved = ComponentVector(nns=(flux_net=ComponentVector(params=copy(flat_params)),))
    @test flux(x, saved) ≈ expected

    perturbed = ComponentVector(nns=(flux_net=ComponentVector(params=flat_params .+ 0.01),))
    @test !isapprox(flux(x, perturbed), expected)

    model64 = Flux.f64(Flux.Chain(Flux.Dense(2 => 2, tanh), Flux.Dense(2 => 1)))
    flux64 = HydroModels.NeuralFlux(
        [flux_x, flux_y], [flux_out], model64; chain_name=:flux_net64
    )
    params64 = HydroModels.get_initial_params(flux64; rng=StableRNG(12))
    x64 = Float64.(x)
    @test eltype(flux64(x64, params64)) == Float64
    @test all(isfinite, flux64(x64, params64))

    @test_throws ArgumentError HydroModels.NeuralFlux(
        [flux_x], [flux_out], Flux.Dense(1 => 1)
    )
    @test_throws ArgumentError HydroModels.NeuralFlux(
        [flux_x], [flux_out], Flux.Chain(Flux.Dropout(0.1)); chain_name=:dropout_net
    )

    bucket = HydroModels.create_simple_neural_bucket(
        Val(:flux);
        name=:flux_bucket,
        n_inputs=1,
        n_states=1,
        n_outputs=1,
        inputs=[:forcing],
        states=[:storage],
        outputs=[:runoff],
    )
    bucket_params = HydroModels.get_initial_params(bucket; rng=StableRNG(13))
    @test HydroModels.get_nn_names(bucket) ==
        [:flux_bucket_flux, :flux_bucket_state, :flux_bucket_output]
    input = Float32[0.2 0.3 0.4]
    bucket_result = bucket(input, bucket_params)
    @test size(bucket_result) == (2, 3)
    @test all(isfinite, bucket_result)

    ode_flux = flux
    ode_params = params
    ode_rhs = (u, p, t) -> vec(ode_flux(reshape(Float32[u[1], t], 2, 1), p))
    ode_config = HydroModels.HydroConfig(solver=HydroModels.ODESolver, solve_alg=Tsit5())
    ode_solution = HydroModels.hydrosolve(
        HydroModels.ODESolver,
        ode_rhs,
        ode_params,
        Float32[0.1],
        Float32[0, 1, 2],
        ode_config,
    )
    @test size(ode_solution) == (1, 3)
    @test all(isfinite, ode_solution)

    discrete_solution = HydroModels.hydrosolve(
        HydroModels.DiscreteSolver,
        ode_rhs,
        ode_params,
        Float32[0.1],
        Float32[0, 1, 2],
        HydroModels.HydroConfig(
            solver=HydroModels.DiscreteSolver, solve_alg=FunctionMap{true}()
        ),
    )
    @test size(discrete_solution) == (1, 3)
    @test all(isfinite, discrete_solution)

    bucket_rhs = function (u, p, t)
        nn_params = HydroModels._get_neural_bucket_params(bucket, p)
        _, next_state, _ = HydroModels._neural_bucket_step(
            bucket,
            Float32[t],
            u,
            nn_params,
            (flux=nothing, state=nothing, output=nothing),
        )
        next_state .- u
    end
    bucket_solution = HydroModels.hydrosolve(
        HydroModels.DiscreteSolver,
        bucket_rhs,
        bucket_params,
        Float32[0.1],
        Float32[0, 1, 2],
        HydroModels.HydroConfig(
            solver=HydroModels.DiscreteSolver, solve_alg=FunctionMap{true}()
        ),
    )
    @test size(bucket_solution) == (1, 3)
    @test all(isfinite, bucket_solution)
end
