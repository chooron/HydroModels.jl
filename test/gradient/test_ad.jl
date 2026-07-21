using Test
using ComponentArrays
using DataInterpolations
using DifferentiationInterface
using ForwardDiff
using Mooncake
using OrdinaryDiffEq
using SciMLBase
using SciMLSensitivity
using HydroModels

@variables forcing storage runoff
@parameters k bias

const gradient_bucket = @hydrobucket begin
    fluxes = begin
        @hydroflux runoff ~ k * forcing
    end
    dfluxes = begin
        @stateflux storage ~ runoff + bias
    end
end

const gradient_input = reshape(collect(1.0:6.0), 1, :)
const gradient_times = collect(1:6)
const gradient_initstates = [0.25]

gradient_params(v) = ComponentVector(params=(k=v[1], bias=v[2]))

struct ScaledInputFlux
    scale::Float64
end

(flux::ScaledInputFlux)(inputs, params) = (flux.scale .* params.k .* inputs[1],)

function algebraic_loss(v)
    output = gradient_bucket(
        gradient_input,
        gradient_params(v),
        (solver=ImmutableSolver, timeidx=gradient_times);
        initstates=gradient_initstates,
    )
    sum(output)
end

function component_merge_loss(θ)
    base = ComponentVector(a=θ.a, b=θ.b)
    override = ComponentVector(b=θ.b_override, c=θ.c)
    merged = merge_componentvectors(base, override)
    return merged.a^2 + merged.b^2 + merged.c^2
end

@testset "Callable interpolation and functional components" begin
    prebuilt = HydroModels.LinearInterpolation(gradient_input, gradient_times)
    @test hydrointerp(HydroModels.LinearInterpolation, gradient_input, gradient_times)(2.5) ≈ prebuilt(2.5)
    @test hydrointerp(prebuilt, zeros(1, 1), [1]) === prebuilt

    output = gradient_bucket(
        gradient_input,
        gradient_params([0.7, 0.1]),
        (solver=ImmutableSolver, interpolator=prebuilt, timeidx=gradient_times),
        initstates=gradient_initstates,
    )
    @test all(isfinite, output)

    flux = HydroFlux(ScaledInputFlux(2.0); inputs=[:forcing], outputs=[:runoff], params=[:k])
    @test flux(gradient_input, ComponentVector(params=(k=0.5,))) == gradient_input
end

@testset "ComponentVector overlays route Mooncake gradients" begin
    θ = ComponentVector(a=2.0, b=3.0, b_override=5.0, c=7.0)
    value, gradient = value_and_gradient(component_merge_loss, AutoMooncake(), θ)
    gradient_cv = ComponentVector(gradient.fields.data, getaxes(θ))

    @test value == 78.0
    @test gradient_cv.a == 4.0
    @test gradient_cv.b == 0.0
    @test gradient_cv.b_override == 10.0
    @test gradient_cv.c == 14.0

    base = ComponentVector(params=(a=1.0, b=2.0), nns=(net=ComponentVector(weight=3.0),))
    override = ComponentVector(params=(b=5.0,), nns=(net=ComponentVector(weight=7.0),))
    merged = merge_componentvectors(base, override; strict=true)
    @test NamedTuple(merged) == (params=(a=1.0, b=5.0), nns=(net=(weight=7.0,),))
    @test_throws ArgumentError merge_componentvectors(base, ComponentVector(params=(unknown=1.0,),); strict=true)
end

function ode_loss(v, interpolator, sense_alg)
    output = gradient_bucket(
        gradient_input,
        gradient_params(v),
        (
            solver=ODESolver,
            interpolator=interpolator,
            timeidx=gradient_times,
            sense_alg=sense_alg,
        );
        initstates=gradient_initstates,
    )
    sum(output)
end

function central_difference(f, x; h=1e-5)
    result = similar(x)
    for i in eachindex(x)
        xp = copy(x)
        xm = copy(x)
        xp[i] += h
        xm[i] -= h
        result[i] = (f(xp) - f(xm)) / (2h)
    end
    result
end

@testset "ForwardDiff and Mooncake gradients" begin
    params = [0.7, 0.1]
    fd = ForwardDiff.gradient(algebraic_loss, params)
    _, mooncake = value_and_gradient(algebraic_loss, AutoMooncake(), params)
    finite = central_difference(algebraic_loss, params)

    @test isapprox(fd, finite; rtol=1e-7, atol=1e-9)
    @test isapprox(mooncake, finite; rtol=1e-7, atol=1e-9)
    @test Base.get_extension(HydroModels, :HydroModelsMooncakeExt) !== nothing
end

@testset "ODE gradients with DataInterpolations" begin
    params = [0.7, 0.1]
    interpolator = DataInterpolations.LinearInterpolation
    fd_sense = ForwardDiffSensitivity()
    mooncake_sense = GaussAdjoint(autojacvec=SciMLSensitivity.MooncakeVJP())

    fd_loss = v -> ode_loss(v, interpolator, fd_sense)
    mooncake_loss = v -> ode_loss(v, interpolator, mooncake_sense)
    finite = central_difference(fd_loss, params; h=1e-4)
    fd = ForwardDiff.gradient(fd_loss, params)
    _, mooncake = value_and_gradient(mooncake_loss, AutoMooncake(), params)

    @test isapprox(fd, finite; rtol=1e-4, atol=1e-6)
    @test isapprox(mooncake, finite; rtol=1e-4, atol=1e-6)

    constant_interpolator = HydroModels.ConstantInterpolation
    constant_output = gradient_bucket(
        gradient_input,
        gradient_params(params),
        (
            solver=ODESolver,
            interpolator=constant_interpolator,
            timeidx=gradient_times,
            sense_alg=mooncake_sense,
        );
        initstates=gradient_initstates,
    )
    @test all(isfinite, constant_output)
end
