using Zygote
using OrdinaryDiffEq
using SciMLSensitivity
using MLUtils
using ComponentArrays
using DataInterpolations

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

tmp_func = (inputs, states, pas) -> begin
    prcp = inputs[1]
    temp = inputs[2]
    lday = inputs[3]
    snowpack = states[1]
    Tmin = pas.params.Tmin
    Df = pas.params.Df
    Tmax = pas.params.Tmax
    snowfall = step_func(Tmin - temp) * prcp
    melt = step_func(temp - Tmax) * min(snowpack, Df * (temp - Tmax))
    return [snowfall - melt]
end


Df, Tmax, Tmin = 2.67, 0.17, -2.09
ps = ComponentVector(params=(Df=Df, Tmax=Tmax, Tmin=Tmin))
ps_axes = getaxes(ps)
u0 = zeros(eltype(ps), 1)

node_input = rand(3, 20)
timeidx = collect(1:20)
itpfunc = LinearInterpolation(node_input, timeidx)

function ode_func!(du, u, p, t)
    du .= tmp_func(itpfunc(t), u, ComponentVector(p, ps_axes))
end

tspan = (1.0, 20.0)
prob = ODEProblem(ode_func!, u0, tspan, Vector(ps))
sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP()))

Zygote.gradient(Vector(ps)) do p
    prob = ODEProblem(ode_func!, u0, tspan, p)
    sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array
    return sum(sol)
end
