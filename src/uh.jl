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
    "calculate weight of unit hydrograph"
    uh_func::Function
    max_lag::Function
    "A named tuple containing information about inputs, outputs, parameters, and states"
    infos::NamedTuple

    function UnitHydrograph(
        inputs::AbstractVector{T}, outputs::AbstractVector{T}, params::AbstractVector{T},
        uh_func::Function, max_lag::Function; name::Optional{Symbol}=nothing,
        kwargs...
    ) where {T<:Num}
        solvetype = get(kwargs, :solvetype, :DISCRETE)
        infos = (; inputs=inputs, outputs=outputs, params=params)
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name
        return new{uh_name, solvetype}(uh_func, max_lag, infos)
    end 

    function UnitHydrograph(
        inputs::AbstractVector{T}, outputs::AbstractVector{T}, params::AbstractVector{T};
        uh_pairs::AbstractVector{<:Pair}, name::Optional{Symbol}=nothing,
        kwargs...
    ) where {T<:Num}
        solvetype = get(kwargs, :solvetype, :DISCRETE)
        max_lag = get(kwargs, :max_lag, uh_pairs[1][1])
        uh_func, max_lag_func = build_uh_func(uh_pairs, params, max_lag)
        @assert length(inputs) == length(outputs) == 1 "only one input and one output is supported"
        @assert solvetype in [:DISCRETE, :SPARSE, :DSP] "solvetype must be one of [:DISCRETE, :SPARSE, :DSP]"
        solvetype == :DSP && @warn "The DSP solver is not supported for Zygote, please use :DISCRETE or :SPARSE instead, and this method required DSP.jl"
        wfunc(pas) = begin
            # weights = map(t -> uh_func(t, pas), 1:max_lag_func(pas))[1:end-1]
            # filter(x -> x > maximum(weights) * min_weight_prop, weights)
            lag_weights = [uh_func(t, pas) for t in 1:max_lag_func(pas)]
            lag_weights = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])
        end
        infos = (; inputs=inputs, outputs=outputs, params=params)
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name
        return new{uh_name,solvetype}(uh_func, max_lag_func, infos)
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
    @assert Meta.isexpr(expr, :block) "Expected a begin...end block"
    uh_func_expr, uh_vars_expr, kwargs_vec = nothing, nothing, Expr[]

    for arg in expr.args
        arg isa LineNumberNode && continue
        if Meta.isexpr(arg, :(=))
            key, value = arg.args[1], arg.args[2]
            if key == :uh_func
                uh_func_expr = value
            elseif key == :uh_vars
                uh_vars_expr = value
            else
                push!(kwargs_vec, Expr(:kw, key, value))
            end
        end
    end

    @assert uh_func_expr !== nothing "Missing uh_func in unit hydrograph definition"
    @assert uh_vars_expr !== nothing "Missing uh_vars in unit hydrograph definition"

    for var_name in extract_variables(uh_func_expr)
        if !@isdefined(var_name)
            expr_str = string(uh_func_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    for var_name in extract_variables(uh_vars_expr)
        if !@isdefined(var_name)
            expr_str = string(uh_vars_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end

    @assert Meta.isexpr(uh_vars_expr, :call) && uh_vars_expr.args[1] == :(=>) "uh_vars must be a single Pair, e.g., P => Q"
    uh_inputs_expr, uh_outputs_expr = Expr(:vect, uh_vars_expr.args[2]), Expr(:vect, uh_vars_expr.args[3])
    uh_pairs, cond_values = Pair[], []
    for uh_expr in uh_func_expr.args
        uh_expr isa LineNumberNode && continue
        if Meta.isexpr(uh_expr, :call) && uh_expr.args[1] == :(=>)
            push!(uh_pairs, uh_expr.args[2] => uh_expr.args[3])
            push!(cond_values, uh_expr.args[3])
        end
    end

    return esc(quote
        HNum = HydroModels.Num
        params_val = reduce(union, map(
            val -> HNum.(filter(x -> HydroModels.isparameter(x), HydroModels.get_variables(val))),
            [$(cond_values...)]
        ))
        UnitHydrograph(
            $uh_inputs_expr, $uh_outputs_expr, params_val;
            uh_pairs=$uh_pairs, max_lag=$(uh_pairs[1][1]), name=$(name),
            $(kwargs_vec...)
        )
    end)
end


"""
    build_uh_func(uh_pairs, params, max_lag)

Builds unit hydrograph functions.
"""
function build_uh_func(uh_pairs::AbstractVector{<:Pair}, params::AbstractVector, max_lag::Number)
    conditions_rev = vcat([0], reverse(first.(uh_pairs)))
    values_rev = reverse(last.(uh_pairs))
    params_assign_calls = generate_param_assignments(params=params)

    values_exprs = map(eachindex(values_rev)) do i
        :(
            if $(toexpr(conditions_rev[i])) < t < $(toexpr(conditions_rev[i+1]))
                return $(toexpr(values_rev[i]))
            end
        )
    end

    uh_func_expr = :(function (t, pas)
        $(params_assign_calls...)
        $(values_exprs...)
        return 0.0
    end)

    max_lag_expr = :(function (pas)
        $(params_assign_calls...)
        return ceil($(toexpr(max_lag)))
    end)

    return @RuntimeGeneratedFunction(uh_func_expr), @RuntimeGeneratedFunction(max_lag_expr)
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

function (uh::UnitHydrograph{N,:DISCRETE})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T,N}
    solver = get(kwargs, :solver, DiscreteSolver())
    interp = get(kwargs, :interp, DirectInterpolation)
    timeidx = collect(1:size(input, 2))
    interp_func = interp(input, timeidx)
    lag_weights = [uh.uh_func(t, pas) for t in 1:uh.max_lag(pas)]
    uh_weight = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])
    if length(uh_weight) == 0
        return input
    else
        update_func(i, u, p) = i .* p .+ [diff(u, dims=1); -u[end]]
        sol = solver((u, p, t) -> stack(update_func.(interp_func(t), eachslice(u, dims=1), Ref(p)), dims=1),
            uh_weight ./ sum(uh_weight),
            zeros(size(input, 1), length(uh_weight)), timeidx
        )
        return sol[:, 1, :]
    end
end

function (uh::UnitHydrograph{N,:SPARSE})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T,N}
    lag_weights = [uh.uh_func(t, pas) for t in 1:uh.max_lag(pas)]
    uh_weight = vcat([lag_weights[1]], (circshift(lag_weights, -1).-lag_weights)[1:end-1])
    if length(uh_weight) == 0
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

function (uh::UnitHydrograph{N,:DSP})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T,N}
    uh_weight = uh.wfunc(pas)
    if length(uh_weight) == 0
        @warn "The unit hydrograph weight is empty, please check the unit hydrograph function"
        return input
    else
        stack(conv.(eachslice(input, dims=1), Ref(uh_weight)), dims=1)
    end
end

function (uh::UnitHydrograph)(input::AbstractArray{T,3}, pas::AbstractArray; kwargs...) where {T}
    #* Extract the initial state of the parameters and routement in the pas variable
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    uh_param_names = get_param_names(uh)
    extract_params = eachrow(reshape(reduce(vcat, pas[:params][uh_param_names]), :, length(uh_param_names))[ptyidx, :])
    extract_params_cv = map(eachindex(ptyidx)) do idx
        ComponentVector(params=NamedTuple{Tuple(uh_param_names)}(extract_params[idx]))
    end
    node_sols = uh.(eachslice(input, dims=2), extract_params_cv)
    sol_mat = reduce((m1, m2) -> cat(m1, m2, dims=1), node_sols)
    return reshape(sol_mat, 1, size(input)[2], size(input)[3])
end