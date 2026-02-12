using Test
using HydroModels
using ComponentArrays
using Lux
using Random

@testset "get_initial_params tests" begin

    @testset "HydroFlux with parameters" begin
        # Create a simple flux with parameters
        @variables prcp q
        @parameters k [bounds = (0.1, 10.0)]
        @parameters b [bounds = (0.0, 5.0)]

        flux = @hydroflux q ~ k * prcp + b

        # Get initial parameters
        params = get_initial_params(flux; rng=Random.MersenneTwister(42))

        # Check structure
        @test haskey(getaxes(params)[1], :params)
        @test haskey(params.params, :k)
        @test haskey(params.params, :b)

        # Check bounds (should be within specified ranges)
        @test 0.1 <= params.params.k <= 10.0
        @test 0.0 <= params.params.b <= 5.0

        println("HydroFlux params: ", params)
    end

    @testset "NeuralFlux with neural network" begin
        @variables prcp temp et

        # Create neural network
        nn = Chain(
            Dense(2 => 10, tanh),
            Dense(10 => 1),
            name=:et_net
        )

        neural_flux = @neuralflux et ~ nn([prcp, temp])

        # Get initial parameters
        params = get_initial_params(neural_flux; rng=Random.MersenneTwister(42))

        # Check structure
        @test haskey(getaxes(params)[1], :nns)
        @test haskey(params.nns, :et_net)
        @test params.nns.et_net isa ComponentVector

        # Check neural network has parameters
        @test length(params.nns.et_net) > 0

        println("NeuralFlux params structure: ", keys(params))
        println("Neural network param count: ", length(params.nns.et_net))
    end

    @testset "Hybrid component with params and nns" begin
        @variables prcp temp snowpack soilwater et q
        @parameters k [bounds = (0.1, 5.0)]
        @parameters alpha [bounds = (0.0, 1.0)]

        # Create neural flux
        nn = Chain(Dense(2 => 8, relu), Dense(8 => 1), name=:et_nn)
        et_flux = @neuralflux et ~ nn([temp, snowpack])

        # Create traditional flux
        q_flux = @hydroflux q ~ k * soilwater + alpha * prcp

        # Create bucket (simplified)
        bucket = @hydrobucket :test_bucket begin
            fluxes = [et_flux, q_flux]
        end

        # Get initial parameters
        params = get_initial_params(bucket; rng=Random.MersenneTwister(42))

        # Check structure
        @test haskey(getaxes(params)[1], :params)
        @test haskey(getaxes(params)[1], :nns)

        # Check traditional parameters
        @test haskey(params.params, :k)
        @test haskey(params.params, :alpha)
        @test 0.1 <= params.params.k <= 5.0
        @test 0.0 <= params.params.alpha <= 1.0

        # Check neural network parameters
        @test haskey(params.nns, :et_nn)
        @test params.nns.et_nn isa ComponentVector

        println("Hybrid bucket params: ", keys(params))
        println("Traditional params: ", keys(params.params))
        println("Neural network params: ", keys(params.nns))
    end

    @testset "get_nn_initial_params" begin
        @variables prcp temp et

        nn = Chain(Dense(2 => 10, tanh), Dense(10 => 1), name=:my_net)
        neural_flux = @neuralflux et ~ nn([prcp, temp])

        # Get only neural network parameters
        nn_params = get_nn_initial_params(neural_flux; rng=Random.MersenneTwister(42))

        @test !isnothing(nn_params)
        @test haskey(nn_params, :my_net)
        @test nn_params.my_net isa ComponentVector

        println("NN params only: ", keys(nn_params))
    end

    @testset "Component without parameters" begin
        @variables prcp q

        # Flux without parameters
        flux = @hydroflux q ~ prcp * 2.0

        params = get_initial_params(flux)

        # Should return empty ComponentVector
        @test params isa ComponentVector
        @test length(params) == 0

        println("Empty params: ", params)
    end

    @testset "Reproducibility with RNG" begin
        @variables prcp q
        @parameters k [bounds = (0.1, 10.0)]

        flux = @hydroflux q ~ k * prcp

        # Get parameters with same seed
        rng1 = Random.MersenneTwister(123)
        rng2 = Random.MersenneTwister(123)

        params1 = get_initial_params(flux; rng=rng1)
        params2 = get_initial_params(flux; rng=rng2)

        @test params1.params.k == params2.params.k

        println("Reproducible params: ", params1.params.k, " == ", params2.params.k)
    end

end
