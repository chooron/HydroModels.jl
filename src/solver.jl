"""
A custom ODEProblem solver
"""
@kwdef struct ODESolver <: AbstractSolver
    alg = Tsit5()
    sensealg = InterpolatingAdjoint()
    reltol = 1e-2
    abstol = 1e-2
    saveat = 1.0
end

function (solver::ODESolver)(
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    #* build problem
    # 虽然pas本身就包含了initstates但是initstates的构建方式因输入会有所不同
    prob = ODEProblem(
        ode_func!,
        initstates,
        (timeidx[1], timeidx[end]),
        pas
    )
    #* solve problem
    sol = solve(
        prob,
        solver.alg,
        saveat=timeidx,
        reltol=solver.reltol,
        abstol=solver.abstol,
        sensealg=solver.sensealg
    )
    if convert_to_array
        if SciMLBase.successful_retcode(sol)
            sol_arr = Array(sol)
        else
            @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
            sol_arr = zeros(size(initstates)..., length(timeidx))
        end
        return sol_arr
    else
        return sol
    end
end

"""
A custom ODEProblem solver
"""
@kwdef struct DiscreteSolver <: AbstractSolver
    alg = FunctionMap{true}()
    sensealg = InterpolatingAdjoint()
end

function (solver::DiscreteSolver)(
    ode_func!::Function,
    params::ComponentVector,
    initstates::AbstractArray,
    timeidx::AbstractVector;
    convert_to_array::Bool=true
)
    #* build problem
    # 虽然pas本身就包含了initstates但是initstates的构建方式因输入会有所不同
    prob = DiscreteProblem(
        ode_func!,
        initstates,
        (timeidx[1], timeidx[end]),
        params
    )
    #* solve problem
    sol = solve(
        prob,
        solver.alg,
        saveat=timeidx
    )
    if convert_to_array
        if SciMLBase.successful_retcode(sol)
            sol_arr = Array(sol)
        else
            @warn "ODE solver failed, please check the parameters and initial states, or the solver settings"
            sol_arr = zeros(size(initstates)..., length(timeidx))
        end
        return sol_arr
    else
        return sol
    end
end

@kwdef struct ManualSolver <: AbstractSolver
    #* 计算效率过差不予考虑
end

function (solver::ManualSolver)(
    ode_func!::Function,
    pas::ComponentVector,
    initstates::AbstractVector,
    timeidx::AbstractVector
)
    T = promote_type(eltype(pas), eltype(initstates))
    init_du = zeros(T, size(initstates))
    states_results = ()
    for t in timeidx
        ode_func!(init_du, initstates, pas, t)
        initstates = initstates .+ init_du
        (states_results..., initstates)
    end
    reduce(hcat, states_results)
end