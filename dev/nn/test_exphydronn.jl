using Lux
using Random
using StableRNGs
using CUDA, cuDNN
using BenchmarkTools
using ComponentArrays

# smooting step function
step_fct(x) = (tanh.(5.0f0 .* x) .+ 1.0f0) .* 0.5f0
# snow precipitation
Ps(P, T, Tmin) = step_fct(Tmin .- T) .* P
# rain precipitation
Pr(P, T, Tmin) = step_fct(T .- Tmin) .* P
# snow melt
M(S0, T, Df, Tmax) = step_fct(T .- Tmax) .* Df .* (T .- Tmax)
# evapotranspiration
function PET(T, Lday)
    @. 29.8f0 * Lday * 0.611f0 * exp((17.3f0 * T) / (T + 237.3f0)) / (T + 273.2f0)
end

function ET(S1, T, Lday, Smax)
    return step_fct(S1 .- Smax) .* PET(T, Lday) .+
           step_fct(Smax .- S1) .* PET(T, Lday) .* (S1 ./ Smax)
end

# base flow
function Qb(S1, f, Smax, Qmax)
    return step_fct(S1 .- Smax) .* Qmax .+
           step_fct(Smax .- S1) .* Qmax .* exp.(-f .* (Smax .- S1))
end

# peak flow
Qs(S1, Smax) = step_fct(S1 .- Smax) .* (S1 .- Smax)

@kwdef struct ExpHydroLuxLayer <: Lux.AbstractRecurrentCell
    hidden_dims::Int
end

function initialodestate(l, input)
    s0 = similar(input, Float32, 1, l.hidden_dims)
    s1 = similar(input, Float32, 1, l.hidden_dims)
    fill!(s0, 0)
    fill!(s1, 0)
    return [s0, s1]
end

function (l::ExpHydroLuxLayer)(input::AbstractArray, ps, st::NamedTuple)
    carry = initialodestate(l, input)
    return LuxCore.apply(l, (input, carry), ps, st)
end

function (l::ExpHydroLuxLayer)(input::Tuple, ps, st::NamedTuple)
    forcing, states = input
    p_vec, t_vec, lday_vec = forcing[1:1, :], forcing[2:2, :], forcing[3:3, :]
    s0_vec, s1_vec = states[1], states[2]
    qout_vec = Qb(s1_vec, ps.f, ps.Smax, ps.Qmax) .+ Qs(s1_vec, ps.Smax)
    melt_vec = M(s0_vec, t_vec, ps.Df, ps.Tmax)
    ds0_vec = Ps(p_vec, t_vec, ps.Tmin) .- melt_vec
    ds1_vec = Pr(p_vec, t_vec, ps.Tmin) .+ melt_vec .- ET(s1_vec, t_vec, lday_vec, ps.Tmax) .- qout_vec
    return (qout_vec, [s0_vec .+ ds0_vec, s1_vec .+ ds1_vec]), st
end

function LuxCore.initialparameters(rng::AbstractRNG, l::ExpHydroLuxLayer)
    return (
        f=rand(rng, Float32, 1, l.hidden_dims) .* 0.1f0,
        Smax=rand(rng, Float32, 1, l.hidden_dims) .* (2000.0f0 .- 100.0f0) .+ 100.0f0,
        Qmax=rand(rng, Float32, 1, l.hidden_dims) .* (50.0f0 .- 10.0f0) .+ 10.0f0,
        Df=rand(rng, Float32, 1, l.hidden_dims) .* 5.0f0,
        Tmax=rand(rng, Float32, 1, l.hidden_dims) .* 3.0f0,
        Tmin=rand(rng, Float32, 1, l.hidden_dims) .* -3.0f0,
    )
end

function ExpHydroLuxModel(hidden_dims::Int)
    return Chain(
        Recurrence(
            ExpHydroLuxLayer(hidden_dims); return_sequence=true, ordering=TimeLastIndex()
        ),
        x -> stack(x; dims=3),
    )
end

dev = gpu_device()
model = ExpHydroLuxModel(10)
ps, st = Lux.setup(Random.default_rng(), model) |> dev
input = rand(3, 10, 10000) |> dev
output = model(input, ps, st)
