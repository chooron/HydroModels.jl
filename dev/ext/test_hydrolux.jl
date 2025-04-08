#* Multi-Nodes在原有的架构下执行起来相对困难,我认为应该在Lux.jl的架构下进行重新设计
using Lux
using Lux: foldl_init
using Random
using StableRNGs
using Optimisers
using CUDA
using cuDNN
using Zygote
using MLUtils
using BenchmarkTools
using CSV
using DataFrames

# smooting step function
step_fct(x) = (tanh.(5.0 .* x) .+ 1.0) .* 0.5
# snow precipitation
Ps(P, T, Tmin) = step_fct(Tmin .- T) .* P
# rain precipitation
Pr(P, T, Tmin) = step_fct(T .- Tmin) .* P
# snow melt
M(S0, T, Df, Tmax) = step_fct(T .- Tmax) .* Df .* (T .- Tmax)
# evapotranspiration
PET(T, Lday) = @. 29.8 * Lday * 0.611 * exp((17.3 * T) / (T + 237.3)) / (T + 273.2)
ET(S1, T, Lday, Smax) = step_fct(S1 .- Smax) .* PET(T, Lday) .+ step_fct(Smax .- S1) .* PET(T, Lday) .* (S1 ./ Smax)
# base flow
Qb(S1, f, Smax, Qmax) = step_fct(S1 .- Smax) .* Qmax .+ step_fct(Smax .- S1) .* Qmax .* exp.(-f .* (Smax .- S1))
# peak flow
Qs(S1, Smax) = step_fct(S1 .- Smax) .* (S1 .- Smax)

@kwdef struct ExpHydroLuxLayer <: Lux.AbstractLuxLayer
    hidden_dims::Int
end

function (l::ExpHydroLuxLayer)(input, ps, st::NamedTuple)
    forcing, states = input
    p_vec, t_vec, lday_vec = forcing[1, :], forcing[2, :], forcing[3, :]
    s0_vec, s1_vec = states[1], states[2]
    qout_vec = Qb(s1_vec, ps.f, ps.Smax, ps.Qmax) .+ Qs(s1_vec, ps.Smax)
    melt_vec = M(s0_vec, t_vec, ps.Df, ps.Tmax)
    ds0_vec = Ps(p_vec, t_vec, ps.Tmin) .- melt_vec
    ds1_vec = Pr(p_vec, t_vec, ps.Tmin) .+ melt_vec .- ET(s1_vec, t_vec, lday_vec, ps.Tmax) .- qout_vec
    return (qout_vec, [s0_vec .+ ds0_vec, s1_vec .+ ds1_vec]), st
end

function LuxCore.initialparameters(rng::AbstractRNG, l::ExpHydroLuxLayer)
    return (f=rand(rng, Float64, l.hidden_dims) .* 0.1,
        Smax=rand(rng, Float64, l.hidden_dims) .* (2000.0 .- 100.0) .+ 100.0,
        Qmax=rand(rng, Float64, l.hidden_dims) .* (50.0 .- 10.0) .+ 10.0,
        Df=rand(rng, Float64, l.hidden_dims) .* 5.0,
        Tmax=rand(rng, Float64, l.hidden_dims) .* 3.0,
        Tmin=rand(rng, Float64, l.hidden_dims) .* -3.0)
end

LuxCore.initialstates(::AbstractRNG, ::ExpHydroLuxLayer) = NamedTuple()

@kwdef struct ExpHydroLuxModel{L1,L2} <: Lux.AbstractLuxContainerLayer{(:exphydro, :fc)}
    exphydro::L1
    fc::L2
    device
end

function ExpHydroLuxModel(hidden_dims::Int, device)
    return ExpHydroLuxModel(ExpHydroLuxLayer(hidden_dims), Lux.Dense(hidden_dims, 1, leakyrelu), device)
end

function LuxCore.initialparameters(rng::AbstractRNG, l::ExpHydroLuxModel)
    return (
        exphydro=LuxCore.initialparameters(rng, l.exphydro),
        fc=LuxCore.initialparameters(rng, l.fc)
    ) |> l.device
end

function LuxCore.initialstates(rng::AbstractRNG, l::ExpHydroLuxModel)
    return (
        exphydro=LuxCore.initialstates(rng, l.exphydro),
        fc=LuxCore.initialstates(rng, l.fc)
    ) |> l.device
end

function initialodestate(l::ExpHydroLuxLayer, dev)
    return [zeros(Float64, l.hidden_dims) |> dev, zeros(Float64, l.hidden_dims) |> dev]
end

function (r::ExpHydroLuxModel)(input, ps, st::NamedTuple)
    function recur_op(::Nothing, input)
        (out, carry), state = LuxCore.apply(r.exphydro, (input, initialodestate(r.exphydro, r.device)), ps.exphydro, st)
        return [out], carry, state
    end
    function recur_op((outputs, carry, state), input)
        (out, carry), state = LuxCore.apply(r.exphydro, (input, carry), ps.exphydro, state)
        return vcat(outputs, [out]), carry, state
    end
    results = foldl_init(recur_op, eachslice(input, dims=3))
    output = stack(first(results), dims=2)
    fc_output, st = LuxCore.apply(r.fc, output, ps.fc, last(results))
    return fc_output, st
end

hidden_dims = 1
dev = cpu_device() # cpu_device gpu_device
model = ExpHydroLuxModel(hidden_dims, dev)
ps, st = Lux.setup(StableRNG(42), model) |> dev

test_data = rand(Float64, hidden_dims) |> dev
tmp_input = rand(Float64, 3, hidden_dims, 10000) |> dev
@btime output, state = model(tmp_input, ps, st)

# CUDA.@profile model(tmp_input, ps, st)
# @btime Zygote.gradient((ps) -> model(tmp_input, ps, st)[1] |> sum, ps)
# Lux.sigmoid

# CUDA.@profile sin.(test_data .+ 1.0) + sin.(test_data .* 2.0)

# layer = ExpHydroLuxLayer(hidden_dims)
# ps, st = Lux.setup(StableRNG(42), layer) |> dev
# input = rand(Float64, 3, hidden_dims) |> dev
# states = [zeros(Float64, hidden_dims) |> dev, zeros(Float64, hidden_dims) |> dev]
# output, st = layer((input, states), ps, st)
# CUDA.@profile layer((input, states), ps, st)

# @btime layer((input, states), ps, st)
# CUDA.@profile Zygote.gradient((ps) -> layer((input, states), ps, st)[1] |> sum, ps)
