using Zygote
using OrdinaryDiffEq
using SciMLSensitivity
using MLUtils
using ComponentArrays
using DataInterpolations
# using RuntimeGeneratedFunctions
# RuntimeGeneratedFunctions.init(@__MODULE__)

# func_expr = :((inputs, states, pas) -> begin
#     #= e:\JlCode\HydroModels\src\utils\build.jl:157 =#
#     #= e:\JlCode\HydroModels\src\utils\build.jl:158 =#
#     Base.@_inline_meta #= e:\JlCode\HydroModels\src\utils\build.jl:127 =#
#     #= e:\JlCode\HydroModels\src\utils\build.jl:159 =#
#     prcp = inputs[1]
#     temp = inputs[2]
#     lday = inputs[3]
#     snowpack = states[1]
#     Tmin = pas.params.Tmin
#     Df = pas.params.Df
#     Tmax = pas.params.Tmax
#     snowfall = broadcast(*, broadcast(*, 0.5, prcp), broadcast(+, 1.0, broadcast(tanh, broadcast(*, 5.0, broadcast(+, Tmin, broadcast(*, -1, temp))))))
#     melt = broadcast(*, broadcast(*, broadcast(*, 0.25, broadcast(+, 1.0, broadcast(tanh, broadcast(*, 5.0, snowpack)))), broadcast(+, 1.0, broadcast(tanh, broadcast(*, 5.0, broadcast(+, broadcast(*, -1, Tmax), temp))))), broadcast(min, snowpack, broadcast(*, Df, broadcast(+, broadcast(*, -1, Tmax), temp))))
#     return stack([broadcast(+, broadcast(*, -1, melt), snowfall)], dims=1)
# end)

# tmp_func = @RuntimeGeneratedFunction(func_expr)

# step_func = (x) -> (tanh(5.0 * x) + 1.0) * 0.5

tmp_func2 = (inputs, states, pas) -> begin
    prcp = inputs[1]
    temp = inputs[2]
    lday = inputs[3]
    snowpack = states[1]
    Tmin = pas[1]
    Df = pas[2]
    Tmax = pas[3]
    snowfall = @. step_func(Tmin - temp) * prcp
    melt = @. step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))
    return stack([snowfall .- melt], dims=1)
end

Df, Tmax, Tmin = 2.674548848, 0.175739196, -2.092959084
ps = ComponentVector(params=(Df=Df, Tmax=Tmax, Tmin=Tmin))
ps_axes = getaxes(ps)
u0 = rand(1, 10)

node_input = rand(3, 10, 20)
itpfuncs = LinearInterpolation.(eachslice(node_input, dims=1), Ref(1:20))
itpfunc = LinearInterpolation(reshape(node_input, :, 20), 1:20)
tmp_input_func = (t) -> ntuple(i -> itpfuncs[i](t), length(itpfuncs))

function multi_ode_func!(du, u, p, t)
    inputs = reshape(itpfunc(t), 3, 10)
    prcp = inputs[1]
    temp = inputs[2]
    lday = inputs[3]
    snowpack = u[1,:]
    Tmin = p[1]
    Df = p[2]
    Tmax = p[3]
    snowfall = @. step_func(Tmin - temp) * prcp
    melt = @. step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))
    du .= stack([snowfall .- melt], dims=1)
    # du .= tmp_func2(eachslice(tmp_input, dims=1), eachslice(u, dims=1), p)
end

tspan = (1.0, 20.0)
prob = ODEProblem(multi_ode_func!, u0, tspan, Vector(ps))
sol = solve(prob, Tsit5(), sensealg=GaussAdjoint(autojacvec=EnzymeVJP())) |> Array

Zygote.gradient(Vector(ps)) do p
    prob = ODEProblem(multi_ode_func!, u0, tspan, p)
    sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array |> sum
end

# function test_func(t, p)
#     tmp_input = reshape(itpfunc(t), 3, 10)
#     tmp_func2(eachslice(tmp_input, dims=1), eachslice(u0, dims=1), p) |> sum
# end

# ps_d = zeros(eltype(ps), length(ps))
# grad = Enzyme.autodiff(Reverse, test_func, Active, Active(3.5), Duplicated(Vector(ps), ps_d))[1]