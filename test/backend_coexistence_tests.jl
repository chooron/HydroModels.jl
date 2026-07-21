using Flux
using Lux
using LuxCore
using SimpleChains
using HydroModels
using Test

@testset "Neural backend coexistence" begin
    @test !isnothing(Base.get_extension(HydroModels, :HydroModelsLuxExt))
    @test !isnothing(Base.get_extension(HydroModels, :HydroModelsFluxExt))
    @test !isnothing(Base.get_extension(HydroModels, :HydroModelsSimpleChainsExt))

    @variables x1 x2 y
    lux_chain = Lux.Chain(layer=Lux.Dense(2 => 1), name=:lux_backend)
    flux_chain = Flux.Chain(Flux.Dense(2 => 1))
    simple_chain = SimpleChain(SimpleChains.static(2), TurboDense(identity, 1))

    lux_flux = HydroModels.NeuralFlux([x1, x2], [y], lux_chain; chain_name=:lux_backend)
    flux_flux = HydroModels.NeuralFlux([x1, x2], [y], flux_chain; chain_name=:flux_backend)
    simple_flux = HydroModels.NeuralFlux(
        [x1, x2], [y], simple_chain; chain_name=:simple_backend
    )

    @test HydroModels.get_nn_names(lux_flux) == [:lux_backend]
    @test HydroModels.get_nn_names(flux_flux) == [:flux_backend]
    @test HydroModels.get_nn_names(simple_flux) == [:simple_backend]
    @test length(HydroModels.get_initial_params(lux_flux).nns.lux_backend) > 0
    @test length(HydroModels.get_initial_params(flux_flux).nns.flux_backend.params) > 0
    @test length(HydroModels.get_initial_params(simple_flux).nns.simple_backend.params) > 0
end
