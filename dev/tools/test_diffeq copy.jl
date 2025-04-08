using SciMLSensitivity, DifferentialEquations, Zygote
using ComponentArrays
using Enzyme
using Optimization
using OptimizationOptimisers

p = ComponentVector(p1=1.5, p2=1.0, p3=3.0, p4=1.0)
p_axes = getaxes(p)

function fiip(du, u, p, t)
    p_ = ComponentVector(p, p_axes)
    du[1] = dx = p_.p1 * u[1] - p_.p2 * u[1] * u[2]
    du[2] = dy = -p_.p3 * u[2] + p_.p4 * u[1] * u[2]
end

u0 = [1.0; 1.0]
prob = ODEProblem(fiip, [1.0; 1.0], (0.0, 10.0), [1.5, 1.0, 3.0, 1.0])
base_sol = solve(prob, Tsit5(), saveat=0.1,) |> Array

objective(u, p) = begin
    u_ = Vector(u)
    sol = solve(prob, Tsit5(), u0=u0, p=u_, saveat=0.1, sensealg=BacksolveAdjoint(autojacvec=EnzymeVJP())) |> Array
    loss = sum((base_sol .- sol) .^ 2)
    println(loss)
    return loss
end

optf = Optimization.OptimizationFunction(objective, AutoZygote())
optprob = OptimizationProblem(optf, ComponentVector(p1=3.5, p2=1.4, p3=3.2, p4=1.5)) # [3.5, 1.4, 3.2, 1.5] ComponentVector(p1=3.5, p2=1.4, p3=3.2, p4=1.5)
optsolve = Optimization.solve(optprob, Adam(), maxiters=1000)
