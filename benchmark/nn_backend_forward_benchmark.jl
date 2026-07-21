using BenchmarkTools
using ComponentArrays
using Flux
using HydroModels
using Lux
using LuxCore
using OrdinaryDiffEq
using Random
using SimpleChains

const N = 16
const X = rand(Float32, N)

lux_chain = Lux.Chain(;
    hidden=Lux.Dense(N => 32, tanh), output=Lux.Dense(32 => 1), name=:benchmark_lux
)
lux_ps = ComponentVector(LuxCore.initialparameters(MersenneTwister(1), lux_chain))
lux_st = LuxCore.initialstates(MersenneTwister(1), lux_chain)

flux_model = Flux.Chain(Flux.Dense(N => 32, tanh), Flux.Dense(32 => 1))
flux_ps, flux_rebuild = Flux.destructure(flux_model)

simple_chain = SimpleChain(
    SimpleChains.static(N), TurboDense(tanh, 32), TurboDense(identity, 1)
)
simple_ps = SimpleChains.init_params(simple_chain, Float32; rng=MersenneTwister(1))

println("parameter counts")
println("  Lux: ", length(lux_ps))
println("  Flux: ", length(flux_ps))
println("  SimpleChains: ", length(simple_ps))

println("\nwarm forward benchmarks")
lux_bench = @benchmark LuxCore.apply($lux_chain, $X, $lux_ps, $lux_st)
flux_bench = @benchmark $flux_rebuild($flux_ps)($X)
simple_bench = @benchmark $simple_chain($X, $simple_ps)
println("Lux direct forward")
display(lux_bench)
println("Flux rebuild-each-call forward")
display(flux_bench)
println("SimpleChains direct forward")
display(simple_bench)

flux_runtime = flux_rebuild(flux_ps)
flux_runtime_bench = @benchmark $flux_runtime($X)
println("\nFlux rebuilt runtime-only benchmark")
display(flux_runtime_bench)

function _ode_input(u, t)
    x = copy(X)
    x[1] = u[1]
    x[2] = Float32(t)
    return x
end

function solve_lux(u, p, t)
    y, _ = LuxCore.apply(lux_chain, _ode_input(u, t), lux_ps, lux_st)
    return y
end

function solve_flux(u, p, t)
    return flux_rebuild(flux_ps)(_ode_input(u, t))
end

function solve_simplechains(u, p, t)
    return simple_chain(_ode_input(u, t), simple_ps)
end

println("\nODE RHS repeated-call benchmarks")
for (label, rhs) in
    (("Lux", solve_lux), ("Flux", solve_flux), ("SimpleChains", solve_simplechains))
    rhs_calls = Ref(0)
    counted_rhs = (u, p, t) -> begin
        rhs_calls[] += 1
        rhs(u, p, t)
    end
    prob = ODEProblem(
        (du, u, p, t) -> (du .= counted_rhs(u, p, t)), Float32[0.1], (0.0f0, 1.0f0), nothing
    )
    solve(prob, Tsit5(); save_everystep=false)
    one_solve_rhs_calls = rhs_calls[]
    bench = @benchmark solve($prob, Tsit5(); save_everystep=false)
    println(label, " (RHS calls in one warm solve: ", one_solve_rhs_calls, ")")
    display(bench)
end

println("\nInterpretation: compare Flux's rebuild-each-call trial with the
runtime-only trial to estimate rebuild overhead. BenchmarkTools output is
machine-dependent and is not a unit-test threshold.")
