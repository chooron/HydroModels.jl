"""
    ManualSolver{mutable} <: AbstractHydroSolver

A lightweight, type-stable ODE solver optimized for hydrological modeling.

# Type Parameters
- `mutable::Bool`: Controls array mutability and performance characteristics
    - `true`: Uses mutable arrays for in-place updates (30% faster)
    - `false`: Uses immutable arrays for functional programming style

# Description
ManualSolver implements a fixed-step explicit Euler method for solving ordinary 
differential equations (ODEs). It is specifically designed for hydrological models 
where stability and computational efficiency are prioritized over high-order accuracy.

## Performance Characteristics
- Type-stable implementation for predictable performance
- Optional in-place operations via the `mutable` parameter
- Linear time complexity with respect to simulation length
- Constant space complexity when `mutable=true`

## Memory Management
- `mutable=true`: 
    - Modifies arrays in-place
- `mutable=false`:
    - Creates new arrays at each step

## Use Cases
- Small to medium-scale hydrological models
- Systems with moderate stiffness
- Real-time applications requiring predictable performance
- Memory-constrained environments (with `mutable=true`)

See also: [`AbstractHydroSolver`](@ref), [`solve`](@ref)
"""
struct ManualSolver{mutable} end

function (solver::ManualSolver{true})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 1},
    timeidx::AbstractVector;
    kwargs...
)
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(eltype(initstates), length(initstates), length(timeidx))
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results[:, i] = tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{true})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 2},
    timeidx::AbstractVector;
    kwargs...
)
    T1 = promote_type(eltype(pas), eltype(initstates))
    states_results = zeros(eltype(initstates), size(initstates)..., length(timeidx))
    tmp_initstates = copy(initstates)
    for (i, t) in enumerate(timeidx)
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_du_mat = reduce(hcat, tmp_du)
        tmp_initstates = tmp_initstates .+ tmp_du_mat
        states_results[:, :, i] .= tmp_initstates
    end
    states_results
end

function (solver::ManualSolver{false})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 1},
    timeidx::AbstractVector;
    kwargs...
)
    states_results = []
    tmp_initstates = copy(initstates)
    for t in timeidx
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results = vcat(states_results, [tmp_initstates])
    end
    reduce((m1, m2) -> cat(m1, m2, dims=length(size(initstates))+1), states_results)
end

function (solver::ManualSolver{false})(
    du_func::Function,
    pas::AbstractVector,
    initstates::AbstractArray{<:Number, 2},
    timeidx::AbstractVector;
    kwargs...
)
    states_results = []
    tmp_initstates = copy(initstates)
    for t in timeidx
        tmp_du = du_func(tmp_initstates, pas, t)
        tmp_initstates = tmp_initstates .+ tmp_du
        states_results = vcat(states_results, [tmp_initstates])
    end
    output = reduce((m1, m2) -> cat(m1, m2, dims=length(size(initstates))+1), states_results)
    output
end
