using DifferentialEquations
using SciMLSensitivity
using Zygote
using Enzyme
using DataInterpolations
using MLUtils

test_input = rand(3, 10, 101)
itp_funcs = LinearInterpolation.(eachslice(test_input, dims=1), Ref(0:100))

function test_func(u, p, t)
    tmp_i = [itpfunc(t * 100) for itpfunc in itp_funcs]
    u0 = u[1, :]
    return reduce(hcat, [p[1] .* u0 .+ tmp_i[1] .- tmp_i[2] .* tmp_i[3]])
end

u0 = ones(1, 10)
tspan = (0.0, 1.0)

prob = ODEProblem(test_func, u0, tspan, [1.03])
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP()))
sol.u[end] |> sum

Zygote.gradient(p -> begin
        test_func(u0, p, 0.5) |> sum
    end, [1.03])
stack(fill(ones(10, 100),3), dims=2)
# Zygote.gradient(p -> begin
#     prob = ODEProblem(test_func, u0, tspan, p)
#     sol = solve(prob, Tsit5(), reltol = 1e-6, abstol = 1e-6, sensealg = BacksolveAdjoint(autojacvec = EnzymeVJP()))
#     sol.u[end] |> sum
# end, [1.03])
