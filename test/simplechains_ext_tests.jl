using ComponentArrays
using OrdinaryDiffEq
using OrdinaryDiffEqFunctionMap
using SciMLBase
using SimpleChains
using HydroModels
using StableRNGs
using Test

@testset "SimpleChains extension" begin
    @test !isnothing(Base.get_extension(HydroModels, :HydroModelsSimpleChainsExt))

    @variables sc_x sc_y sc_out
    chain = SimpleChain(
        SimpleChains.static(2), TurboDense(tanh, 3), TurboDense(identity, 1)
    )
    flux = HydroModels.NeuralFlux(
        [sc_x, sc_y], [sc_out], chain; chain_name=:simple_net, norm=x -> x
    )
    @test HydroModels.get_nn_names(flux) == [:simple_net]

    params = HydroModels.get_initial_params(flux; rng=StableRNG(21))
    flat_params = params.nns.simple_net.params
    expected_params = SimpleChains.init_params(chain, Float32; rng=StableRNG(21))
    @test length(flat_params) == length(expected_params)
    @test all(isfinite, flat_params)

    x = Float32[0.2 0.4 -0.1; 1.0 -0.5 0.7]
    expected = chain(x, flat_params)
    @test Array(flux(x, params)) ≈ Array(expected)

    x64 = Float64.(x)
    @test all(isfinite, flux(x64, params))
    @test length(
        HydroModels.get_initial_params(flux; rng=StableRNG(21)).nns.simple_net.params
    ) == length(flat_params)

    saved = ComponentVector(nns=(simple_net=ComponentVector(params=copy(flat_params)),))
    @test Array(flux(x, saved)) ≈ Array(expected)
    perturbed = ComponentVector(
        nns=(simple_net=ComponentVector(params=flat_params .+ 0.01),)
    )
    @test !isapprox(Array(flux(x, perturbed)), Array(expected))

    @test_throws ArgumentError HydroModels.NeuralFlux(
        [sc_x], [sc_out], SimpleChain(SimpleChains.static(1), TurboDense(identity, 1))
    )

    bucket = HydroModels.create_neural_bucket(
        Val(:simplechains);
        name=:simple_bucket,
        n_inputs=1,
        n_states=1,
        n_outputs=1,
        hidden_size=3,
        inputs=[:forcing],
        states=[:storage],
        outputs=[:runoff],
    )
    bucket_params = HydroModels.get_initial_params(bucket; rng=StableRNG(22))
    @test HydroModels.get_nn_names(bucket) ==
        [:simple_bucket_flux, :simple_bucket_state, :simple_bucket_output]
    @test size(bucket(Float32[0.2 0.3 0.4], bucket_params)) == (2, 3)

    rhs = (u, p, t) -> vec(flux(reshape(Float32[u[1], t], 2, 1), p))
    solution = HydroModels.hydrosolve(
        HydroModels.ODESolver,
        rhs,
        params,
        Float32[0.1],
        Float32[0, 1, 2],
        HydroModels.HydroConfig(solver=HydroModels.ODESolver, solve_alg=Tsit5()),
    )
    @test size(solution) == (1, 3)
    @test all(isfinite, solution)

    discrete_solution = HydroModels.hydrosolve(
        HydroModels.DiscreteSolver,
        rhs,
        params,
        Float32[0.1],
        Float32[0, 1, 2],
        HydroModels.HydroConfig(
            solver=HydroModels.DiscreteSolver, solve_alg=FunctionMap{true}()
        ),
    )
    @test size(discrete_solution) == (1, 3)

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
end
