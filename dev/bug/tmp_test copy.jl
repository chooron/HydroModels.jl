using Zygote
using OrdinaryDiffEq
using SciMLSensitivity
using MLUtils

# node_input = rand(3, 10)
# itpfunc = LinearInterpolation(node_input, 1:10)

# function lorenz(u, p, t)
#     tmp_input = itpfunc(t)
#     dx = 10.0 * (u[2] - u[1]) * p[1] * tmp_input[1]
#     dy = u[1] * (28.0 - u[3]) - u[2] * p[2] * tmp_input[2]
#     dz = u[1] * u[2] - (8 / 3) * u[3] * p[3] * tmp_input[3]
#     [dx, dy, dz]
# end

# tspan = (1.0, 10.0)
# ps = rand(3)
# u0 = rand(3)
# prob = ODEProblem(lorenz, u0, tspan, ps)
# sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array |> sum

# Zygote.gradient(ps) do p
#     prob = ODEProblem(lorenz, u0, tspan, p)
#     sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array |> sum
# end


node_input = rand(3, 5, 10)
itpfuncs = LinearInterpolation.(eachslice(node_input, dims=1), Ref(1:10))
itpfuncs1 = LinearInterpolation(node_input[1, :, :], 1:10)
itpfuncs2 = LinearInterpolation(node_input[2, :, :], 1:10)
itpfuncs3 = LinearInterpolation(node_input[3, :, :], 1:10)

function multi_snow_func(u, p, t)
    # tmp_input = [itpfunc(t) for itpfunc in itpfuncs]
    tmp_input = (itpfuncs1(t), itpfuncs2(t), itpfuncs3(t))
    dx = @. 10.0 * (u[2, :] - u[1, :]) * p[1, :] * tmp_input[1]
    dy = @. u[1, :] * (28.0 - u[3, :]) - u[2, :] * p[2, :] * tmp_input[2]
    dz = @. u[1, :] * u[2, :] - (8 / 3) * u[3, :] * p[3, :] * tmp_input[3]
    stack([dx, dy, dz], dims=1)
end

tspan = (1.0, 10.0)
ps = rand(3, 5)
u0 = rand(3, 5)
prob = ODEProblem(lorenz, u0, tspan, ps)
sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array |> sum

Zygote.gradient(ps) do p
    prob = ODEProblem(lorenz, u0, tspan, p)
    sol = solve(prob, Tsit5(), sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array |> sum
end