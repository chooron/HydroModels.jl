"""
    UnitHydrograph{N, ST} <: AbstractHydrograph

Represents a Unit Hydrograph routing component.

# Fields
- `uh_func::Function`: Function defining the unit hydrograph shape based on time and parameters.
- `max_lag_func::Function`: Function calculating the maximum lag time based on parameters.
- `infos::NamedTuple`: Metadata (inputs, outputs, params).

# Constructor
```julia
UnitHydrograph(inputs, params; uh_pairs, max_lag, outputs=[], configs=(solvetype=:DISCRETE, suffix=:_lag), name=nothing)
```
- `inputs`, `params`, `outputs`: Vectors of `Num` variables.
- `uh_pairs`: Vector of `Pair` defining piecewise UH segments (`time_bound => expression`).
- `max_lag`: Upper bound for the first UH segment.
- `configs`: NamedTuple containing `solvetype` (`:DISCRETE` or `:SPARSE`) and `suffix` for default output names.
- `name`: Optional symbol identifier.

# Notes
- The type parameter `ST` reflects the `solvetype` (`:DISCRETE` or `:SPARSE`).
- Generates `uh_func` and `max_lag_func` from `uh_pairs` and `params`.
- Used to convolve input time series with the defined unit hydrograph.
"""
struct UnitHydrograph{N,ST} <: AbstractHydrograph
    "The unit hydrograph function"
    uh_func::Function
    "calculate max lag"
    max_lag_func::Function
    "A named tuple containing information about inputs, outputs, parameters, and states"
    infos::NamedTuple

    function UnitHydrograph(
        inputs::AbstractVector{T}, params::AbstractVector{T};
        uh_func::Function, max_lag_func::Function,
        configs::NamedTuple=(solvetype=:DISCRETE, suffix=:_lag, outputs=Num[]),
        name::Union{Symbol,Nothing}=nothing,
    ) where {T<:Num}
        #* Setup the name information of the hydroroutement
        input_names = tosymbol.(inputs)
        output_names = if length(get(configs, :outputs, Num[])) == 0
            Symbol.(input_names, get(configs, :suffix, :_lag))
        else
            tosymbol.(get(configs, :outputs, Num[]))
        end
        outputs = map(name -> only(@variables $name), output_names)
        infos = (; inputs=inputs, outputs=outputs, params=params)
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name
        solvetype = get(configs, :solvetype, :DISCRETE)
        @assert solvetype in [:DISCRETE, :SPARSE] "solvetype must be one of [:DISCRETE, :SPARSE]"
        return new{uh_name,solvetype}(uh_func, max_lag_func, infos)
    end

    function UnitHydrograph(
        inputs::AbstractVector{T}, params::AbstractVector{T},
        uh_pairs::AbstractVector{<:Pair}, max_lag=uh_pairs[1][1],
        configs::NamedTuple=(solvetype=:DISCRETE, suffix=:_lag, outputs=Num[]),
        name::Union{Symbol,Nothing}=nothing,
    ) where {T<:Num}
        uh_func, max_lag_func = build_uh_func(uh_pairs, params, max_lag)
        return UnitHydrograph(inputs, params, uh_func=uh_func, max_lag_func=max_lag_func, configs=configs, name=name)
    end
end

"""
    @unithydro name begin ... end

Macro to simplify the construction of a `UnitHydrograph` object.

# Usage
```julia
@unithydro :my_uh begin
    uh_func = begin
        # Define piecewise UH function: max_time => expression
        lag => 0.5 * (t / lag)^2.5
        2lag => 1.0 - 0.5 * abs(2 - t / lag)^2.5
    end
    uh_vars = [runoff] # Input variable(s)
    configs = (solvetype=:SPARSE, suffix=:_routed)
end
```
Defines a `UnitHydrograph` with the specified name, piecewise UH function definition (`uh_func`), input variables (`uh_vars`), and configuration (`configs`) within the `begin...end` block.

# Arguments
- `name`: (Optional) Symbol for the unit hydrograph name.
- `block`: A `begin...end` block containing assignments for `uh_func`, `uh_vars`, and optionally `configs`.

# Returns
- A `UnitHydrograph` instance.
"""
macro unithydro(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after unit hydrograph name"

    uh_func_expr, uh_vars_expr, configs_expr = nothing, nothing, nothing
    for arg in expr.args
        if arg isa LineNumberNode
            continue
        elseif Meta.isexpr(arg, :(=)) && arg.args[1] == :uh_func
            uh_func_expr = arg.args[2]
        elseif Meta.isexpr(arg, :(=)) && arg.args[1] == :uh_vars
            uh_vars_expr = arg.args[2]
        elseif Meta.isexpr(arg, :(=)) && arg.args[1] == :configs
            configs_expr = arg.args[2]
        end
    end

    @assert uh_func_expr !== nothing "Missing uh_func in unit hydrograph definition"
    @assert uh_vars_expr !== nothing "Missing uh_vars in unit hydrograph definition"
    @assert Meta.isexpr(uh_func_expr, :block) "Expected a begin...end block for uh_func"

    uh_pairs, cond_values = Pair[], []
    for expr in uh_func_expr.args
        if expr isa LineNumberNode
            continue
        elseif expr.args[1] == :(=>)
            push!(uh_pairs, expr.args[2] => expr.args[3])
            push!(cond_values, expr.args[3])
        end
    end

    configs = configs_expr !== nothing ? configs_expr : default_configs
    params_expr = Expr(:call, :reduce, :union, Expr(:call, :map,
        Expr(:->, :val, quote
            HydroModels.Num.(filter(x -> HydroModels.isparameter(x), HydroModels.get_variables(val)))
        end),
        Expr(:vect, cond_values...)
    ))

    return esc(quote
        UnitHydrograph(
            $uh_vars_expr, $params_expr,
            $uh_pairs, $(uh_pairs[1][1]),
            $configs, $(name)
        )
    end)
end

"""
    (flux::UnitHydrograph)(input::AbstractArray, pas::ComponentVector; kwargs...)

Applies the Unit Hydrograph convolution to the input time series.

# Arguments
- `input::AbstractArray{T,D}`: Input data. Shape can be `(variables, timesteps)` (D=2) or `(variables, nodes, timesteps)` (D=3).
- `pas::ComponentVector`: Parameters for the unit hydrograph function (`uh_func`).

# Keyword Arguments
- (for `:DISCRETE` solver): `solver`, `timeidx`, `interp`.
- (for 3D input): `ptyidx` mapping parameters to nodes.

# Returns
- `AbstractArray`: Routed output, with shape matching the input's time (and node, if D=3) dimensions: `(output_variables, timesteps)` or `(output_variables, nodes, timesteps)`.

# Example
```julia
# Assuming 'uh' is a UnitHydrograph instance
routed_flow = uh(input_runoff, parameters)
```

# Notes
- Behavior depends on `solvetype` (`:DISCRETE` or `:SPARSE`) set during construction.
- `:DISCRETE` uses an ODE solver approach (requires `solver` kwarg).
- `:SPARSE` uses direct convolution via sparse matrices.
- 3D input requires `ptyidx` and internally calls the 2D methods per node.
- Single time point input (Vector) is not supported.
"""
(::UnitHydrograph)(::AbstractVector, ::ComponentVector; kwargs...) = @error "UnitHydrograph is not support for single timepoint"

function (flux::UnitHydrograph{N,:DISCRETE})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T,N}
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    interp = get(kwargs, :interp, DirectInterpolation)
    interp_func = interp(input, timeidx)
    #* prepare the initial states
    uh_weight = map(t -> flux.uh_func(t, pas), 1:flux.max_lag_func(pas))[1:end-1]
    if length(uh_weight) == 0
        @warn "The unit hydrograph weight is empty, please check the unit hydrograph function"
        return input
    else
        update_func(i, u, p) = i .* p .+ [diff(u, dims=1); -u[end]]
        #* solve the problem
        sol = solver((u, p, t) -> stack(update_func.(interp_func(t), eachslice(u, dims=1), Ref(p)), dims=1),
            uh_weight ./ sum(uh_weight),
            zeros(size(input, 1), length(uh_weight)), timeidx
        )
        return sol[:, 1, :]
    end
end

function (flux::UnitHydrograph{N,:SPARSE})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T,N}
    uh_weight = map(t -> flux.uh_func(t, pas), 1:flux.max_lag_func(pas))[1:end-1]

    if length(uh_weight) == 0
        @warn "The unit hydrograph weight is empty, please check the unit hydrograph function"
        return input
    else
        function sparse_compute(input_vec)
            #* the weight of the unit hydrograph is normalized by the sum of the weights
            uh_result = [-(i - 1) => uh_wi .* input_vec ./ sum(uh_weight) for (i, uh_wi) in enumerate(uh_weight)]
            #* sum the matrix
            sum(spdiagm(uh_result...), dims=2)[1:end-length(uh_weight)+1]
        end
        return stack(sparse_compute.(eachslice(input, dims=1)), dims=1)
    end
end

function (uh::UnitHydrograph)(input::AbstractArray{T,3}, pas::AbstractArray; kwargs...) where {T}
    #* Extract the initial state of the parameters and routement in the pas variable
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    uh_param_names = get_param_names(uh)
    extract_params = stack(pas[:params][uh_param_names], dims=1)[ptyidx, :]
    extract_params_cv = map(eachindex(ptyidx)) do idx
        ComponentVector(params=NamedTuple{Tuple(uh_param_names)}(extract_params[idx, :]))
    end
    node_sols = uh.(eachslice(input, dims=2), extract_params_cv)
    sol_mat = reduce((m1, m2) -> cat(m1, m2, dims=1), node_sols)
    return reshape(sol_mat, 1, size(input)[2], size(input)[3])
end