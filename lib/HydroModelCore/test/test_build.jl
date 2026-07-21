using Test
using HydroModelCore
using Symbolics

@testset "Build configuration" begin
    @test SAFE_CONFIG.mode == Safe
    @test FAST_CONFIG.mode == Fast
    @test AUTODIFF_CONFIG.mode == AutoDiff
    @test HydroModelCore.DEFAULT_BUILD_CONFIG === SAFE_CONFIG
    @test HydroModelCore.is_ad_safe(SAFE_CONFIG)
    @test HydroModelCore.is_ad_safe(AUTODIFF_CONFIG)
    @test !HydroModelCore.is_ad_safe(FAST_CONFIG)
end

@testset "Dimension and assignment expressions" begin
    @test HydroModelCore.make_index_expr(:inputs, 1, Dim0()) == :(inputs[1])
    @test HydroModelCore.make_index_expr(:inputs, 1, Dim1()) == :(inputs[1, :])
    @test HydroModelCore.make_index_expr(:inputs, 1, Dim2()) == :(inputs[1, :, :])

    safe = HydroModelCore.generate_var_assignments([:x], :inputs, Dim1(), SAFE_CONFIG)
    fast = HydroModelCore.generate_var_assignments([:x], :inputs, Dim1(), FAST_CONFIG)
    @test safe == [:(x = inputs[1, :])]
    @test fast[1].head == :macrocall
    @test fast[1].args[1] == Symbol("@inbounds")
end

@testset "Generated scalar function" begin
    @variables x k
    infos = HydroModelCore.HydroInfos(
        inputs=[:x], outputs=[:y], params=[:k], states=Symbol[], nns=Symbol[]
    )
    flux_func = build_flux_func([k * x], infos)
    @test flux_func([3.0], (params=(k=2.0,),)) == [6.0]
end

@testset "Build analysis" begin
    config = HydroModelCore.FunctionBuildConfig(
        :((inputs, pas)),
        HydroModelCore.generate_var_assignments([:x], :inputs, Dim0()),
        [:(y = inputs[1])],
        :(return [y])
    )
    stats = HydroModelCore.analyze_generated_code(config)
    @test stats[:num_assignments] == 1
    @test stats[:num_computations] == 1
    @test stats[:ad_safe] == true
end
