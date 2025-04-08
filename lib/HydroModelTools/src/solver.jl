"""
A custom ODEProblem solver
"""
@kwdef struct ODESolver <: AbstractHydroSolver
    alg = Tsit5()
    sensealg = BacksolveAdjoint(autojacvec=EnzymeVJP())
    dev = identity
    kwargs = Dict(:reltol=>1e-3, :abstol=>1e-3)
end

function (solver::ODESolver)(
    du_func::Function,
    params::AbstractArray,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    ode_func! = (du, u, p, t) -> (du[:] = du_func(u, p, t))

    #* build problem
    prob = ODEProblem{true}(ode_func!, initstates, (timeidx[1], timeidx[end]), params)
    #* solve problem
    sol = solve(
        prob, solver.alg, saveat=timeidx,
        sensealg=solver.sensealg;
        solver.kwargs...
    )
    if convert_to_array
        if SciMLBase.successful_retcode(sol)
            sol_arr = Array(sol)
        else
            @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
            sol_arr = zeros(size(initstates)..., length(timeidx))
        end
        return sol_arr |> solver.dev
    else
        return sol
    end
end

"""
A custom ODEProblem solver
"""
@kwdef struct DiscreteSolver <: AbstractHydroSolver
    alg = FunctionMap{true}()
    dev = identity
end

function (solver::DiscreteSolver)(
    du_func::Function,
    params::AbstractArray,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    ode_func! = (du, u, p, t) -> (du[:] = du_func(u, p, t))
    #* build problem
    prob = DiscreteProblem(ode_func!, initstates, (timeidx[1], timeidx[end]), params)
    #* solve problem
    sol = solve(prob, solver.alg, saveat=timeidx) # , sensealg=solver.sensealg
    if convert_to_array
        if SciMLBase.successful_retcode(sol)
            sol_arr = Array(sol)
        else
            @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
            sol_arr = zeros(size(initstates)..., length(timeidx))
        end
        return sol_arr |> solver.dev
    else
        return sol
    end
end