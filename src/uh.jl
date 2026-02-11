"""
Unit Hydrograph module - defines unit hydrograph routing components using convolution.
"""

"""
    UnitHydrograph{UF, MF, HT, NT} <: AbstractHydrograph

Represents a Unit Hydrograph routing component using convolution.

When `htypes` is `nothing`, this is a single-node component that accepts 2D input.
When `htypes` is a `Vector{Int}`, this is a multi-node component that accepts 3D input.

$(FIELDS)

# Type Parameters
- `UF`: UH function type
- `MF`: Maximum lag function type
- `HT`: HRU type (`Nothing` for 2D, `Vector{Int}` for 3D)
- `NT`: Metadata type
"""
struct UnitHydrograph{UF,MF,HT,NT} <: AbstractHydrograph
    "unit hydrograph name"
    name::Symbol
    "calculate weight of unit hydrograph"
    uh_func::UF
    "calculate the max lag of the unit hydrograph"
    max_lag::MF
    "HRU type indices (Nothing = single-node 2D, Vector{Int} = multi-node 3D)"
    htypes::HT
    "A named tuple containing information about inputs, outputs, parameters, and states"
    infos::NT

    # Direct function constructor
    function UnitHydrograph(
        inputs::AbstractVector,
        outputs::AbstractVector,
        params::AbstractVector,
        uh_func::Function,
        max_lag_func::Function;
        name::Optional{Symbol}=nothing,
        htypes::Optional{Vector{Int}}=nothing,
        kwargs...
    )
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name

        return new{typeof(uh_func),typeof(max_lag_func),typeof(htypes),typeof(infos)}(
            uh_name, uh_func, max_lag_func, htypes, infos
        )
    end

    # Symbolic constructor
    function UnitHydrograph(
        inputs::AbstractVector,
        outputs::AbstractVector,
        params::AbstractVector;
        uh_conds::AbstractVector{<:Pair},
        name::Optional{Symbol}=nothing,
        kwargs...
    )
        max_lag = get(kwargs, :max_lag, uh_conds[1][1])
        htypes = get(kwargs, :htypes, nothing)
        # Support old-style Vector{Int} htypes (convert empty to nothing)
        if htypes isa Vector{Int} && isempty(htypes)
            htypes = nothing
        end
        param_names = !isempty(params) ? tosymbol.(Num.(params)) : Symbol[]

        uh_func, max_lag_func = build_uh_func(uh_conds, param_names, max_lag)

        @assert length(inputs) == length(outputs) == 1 "Only one input and one output is supported"

        infos = HydroInfos(
            params=param_names,
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name

        return new{typeof(uh_func),typeof(max_lag_func),typeof(htypes),typeof(infos)}(
            uh_name, uh_func, max_lag_func, htypes, infos
        )
    end
end

# ============================================================================
# @unithydro macro
# ============================================================================

"""
    @unithydro [name] begin ... end

Macro to simplify the construction of a UnitHydrograph.

# Usage
```julia
@variables P, Q, lag
@unithydro :my_uh begin
    uh_func = begin
        lag => 0.5 * (t / lag)^2.5
        2lag => 1.0 - 0.5 * abs(2 - t / lag)^2.5
    end
    uh_vars = P => Q
end
```
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

    # Check variable definitions
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

    uh_inputs_expr = Expr(:vect, uh_vars_expr.args[2])
    uh_outputs_expr = Expr(:vect, uh_vars_expr.args[3])

    uh_conds_pair, cond_values = Pair[], []
    for uh_expr in uh_func_expr.args
        uh_expr isa LineNumberNode && continue
        if Meta.isexpr(uh_expr, :call) && uh_expr.args[1] == :(=>)
            push!(uh_conds_pair, uh_expr.args[2] => uh_expr.args[3])
            push!(cond_values, uh_expr.args[2])
            push!(cond_values, uh_expr.args[3])
        end
    end

    return esc(quote
        params_val = reduce(union, map(
            val -> filter(x -> HydroModels.isparameter(x), HydroModels.get_variables(val)),
            [$(cond_values...)]
        )) |> collect
        UnitHydrograph(
            $uh_inputs_expr, $uh_outputs_expr, params_val;
            uh_conds=$uh_conds_pair, max_lag=$(uh_conds_pair[1][1]), name=$(name),
            $(kwargs_vec...)
        )
    end)
end

# ============================================================================
# Convolution computation
# ============================================================================

"""
    hydro_conv(weights, input)

Efficient unit hydrograph convolution using NNlib.conv.
"""
function hydro_conv(weights::Vector{T}, input::Vector{T}) where T
    weights = weights / sum(weights)
    K = length(weights)
    if K == 0
        return T[]
    end
    N = length(input)
    if N == 0
        return zeros(T, N + 1)
    end
    input_reshaped = reshape(input, N, 1, 1)
    weights_reshaped = reshape(weights, K, 1, 1)
    output_reshaped = conv(input_reshaped, weights_reshaped; pad=(K - 1, 0), flipped=false)
    output = vec(output_reshaped)
    return output
end

# ============================================================================
# Functor methods
# ============================================================================

# 2D computation (single-node, htypes = Nothing)
function (uh::UnitHydrograph{UF,MF,Nothing,NT})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {UF,MF,NT,T}
    params = _as_componentvector(params)
    # Compute UH weights (Zygote-compatible)
    max_lag_val = uh.max_lag(params)
    lag_weights = [uh.uh_func(t, params) for t in 1:max_lag_val]

    uh_weight = if isempty(lag_weights)
        T[]
    else
        n = length(lag_weights)
        weights = similar(lag_weights)
        weights[1] = lag_weights[1]
        for i in 2:n
            weights[i] = lag_weights[i] - lag_weights[i-1]
        end
        weights
    end

    if length(uh_weight) == 0
        return input
    end
    uh_out = hydro_conv(uh_weight, input[1, :])
    return reshape(uh_out, 1, :)
end

# 3D computation (multi-node, htypes = Vector{Int})
function (uh::UnitHydrograph{UF,MF,Vector{Int},NT})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {UF,MF,NT,T}
    params = _as_componentvector(params)
    ptyidx = uh.htypes
    uh_param_names = get_param_names(uh)

    # Extract parameters for each node
    params_matrix = reshape(
        reduce(vcat, params[:params][uh_param_names]),
        :, length(uh_param_names)
    )[ptyidx, :]

    extract_params_cv = map(1:length(ptyidx)) do idx
        ComponentVector(params=NamedTuple{Tuple(uh_param_names)}(params_matrix[idx, :]))
    end

    # Apply UH to each node (using 2D method)
    node_sols = map(1:size(input, 2)) do i
        uh_2d = UnitHydrograph(
            collect(get_input_names(uh)),
            collect(get_output_names(uh)),
            collect(get_param_names(uh)),
            uh.uh_func, uh.max_lag;
            name=uh.name
        )
        uh_2d(input[:, i, :], extract_params_cv[i], config)
    end

    sol_mat = cat(node_sols..., dims=1)
    reshape(sol_mat, 1, size(input, 2), size(input, 3))
end

# Error: single-node UH receiving 3D input
function (uh::UnitHydrograph{UF,MF,Nothing,NT})(
    input::AbstractArray{T,3},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {UF,MF,NT,T}
    error("UnitHydrograph without htypes only accepts 2D input (variables × time).\n" *
          "For multi-node computation, provide htypes.\n" *
          "Got input shape: $(size(input))")
end

# Error: multi-node UH receiving 2D input
function (uh::UnitHydrograph{UF,MF,Vector{Int},NT})(
    input::AbstractArray{T,2},
    params::AbstractVector,
    config::ConfigType=default_config();
    kwargs...
) where {UF,MF,NT,T}
    error("UnitHydrograph with htypes only accepts 3D input (variables × nodes × time).\n" *
          "For single-node computation, omit htypes.\n" *
          "Got input shape: $(size(input))")
end

# Export interfaces
export UnitHydrograph, @unithydro
