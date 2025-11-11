"""
Unit Hydrograph module - defines unit hydrograph routing components, supporting both symbolic and functional construction approaches.
"""

"""
    UnitHydrograph{MS, UF, MF, HT, NT} <: AbstractHydrograph

Represents a Unit Hydrograph routing component for convolution.

$(FIELDS)

# Type Parameters
- `MS`: Whether in multi-node mode
- `UF`: UH function type
- `MF`: Maximum lag function type
- `HT`: HRU type vector type
- `NT`: Metadata type
"""
struct UnitHydrograph{MS,UF,MF,HT,NT} <: AbstractHydrograph
    "unit hydrograph name"
    name::Symbol
    "calculate weight of unit hydrograph"
    uh_func::UF
    "calculate the max lag of the unit hydrograph"
    max_lag::MF
    "nodes type"
    hru_types::HT
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
        hru_types::Vector{Int}=Int[],
        kwargs...
    )
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name

        return new{length(hru_types) > 1,typeof(uh_func),typeof(max_lag_func),typeof(hru_types),typeof(infos)}(
            uh_name, uh_func, max_lag_func, hru_types, infos
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
        hru_types = get(kwargs, :hru_types, Int[])
        param_names = !isempty(params) ? tosymbol.(Num.(params)) : Symbol[]

        uh_func, max_lag_func = build_uh_func(uh_conds, param_names, max_lag)

        @assert length(inputs) == length(outputs) == 1 "Only one input and one output is supported"

        infos = HydroInfos(
            params=param_names,
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name

        return new{length(hru_types) > 1,typeof(uh_func),typeof(max_lag_func),typeof(hru_types),typeof(infos)}(
            uh_name, uh_func, max_lag_func, hru_types, infos
        )
    end
end

"""
    @unithydro [name] begin ... end

Macro to simplify the construction of a UnitHydrograph.

# Usage
The macro takes an optional name and a begin...end block that defines the hydrograph.
In the `uh_func` block, `t` is a special variable representing time, and other variables are treated as parameters.

```jldoctest
julia> @variables P, Q, lag  # Define symbolic variables
julia> @unithydro :my_uh begin
           uh_func = begin
               # Piecewise UH definition: max_time => expression
               lag => 0.5 * (t / lag)^2.5
               2lag => 1.0 - 0.5 * abs(2 - t / lag)^2.5
           end
           uh_vars = P => Q  # Defines input => output variables
       end
```

# Arguments
- `name`: (Optional) A Symbol for the name of the UnitHydrograph.
- The begin...end block must contain:
    - `uh_func`: A block of max_time => expression pairs defining the piecewise unit hydrograph. The variable `t` can be used to represent time.
    - `uh_vars`: A Pair of input_variable => output_variable.
- Other assignments in the block are passed as keyword arguments to the UnitHydrograph constructor.

# Returns
- An instance of UnitHydrograph.
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

"""
    (uh::UnitHydrograph)(input, params, config; kwargs...)

Apply the unit hydrograph convolution to an input time series. This is the functor implementation for UnitHydrograph.

This method is dispatched based on the dimensionality of the input array and the MS type parameter 
of the struct, which indicates if it's a multi-node setup.

- **2D Input:** Performs convolution for a single location.
- **3D Input:** Performs convolution for multiple locations (nodes). It iterates over the nodes and applies the 2D method for each.

A method for 1D AbstractVector input is defined to throw an error, as single time points are not supported.
"""
function (uh::UnitHydrograph)(
    input::AbstractArray{T,2},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,2} where {T}
    # Compute UH weights (avoiding in-place modification)
    max_lag_val = uh.max_lag(params)
    lag_weights = [uh.uh_func(t, params) for t in 1:max_lag_val]

    # Compute difference weights
    uh_weight = if isempty(lag_weights)
        T[]
    else
        shifted = circshift(lag_weights, -1)
        vcat([lag_weights[1]], (shifted.-lag_weights)[1:end-1])
    end

    # If no weights, return input directly
    if length(uh_weight) == 0
        return input
    end
    uh_out = hydro_conv(uh_weight, input[1, :])
    return reshape(uh_out, 1, :)
end

# Multi-node UH
function (uh::UnitHydrograph{true})(
    input::AbstractArray{T,3},
    params::ComponentVector,
    config::ConfigType=default_config();
    kwargs...
)::AbstractArray{T,3} where {T}
    ptyidx = uh.hru_types
    uh_param_names = get_param_names(uh)

    # Extract parameters for each node
    params_matrix = reshape(
        reduce(vcat, params[:params][uh_param_names]),
        :, length(uh_param_names)
    )[ptyidx, :]

    extract_params_cv = map(1:length(ptyidx)) do idx
        ComponentVector(params=NamedTuple{Tuple(uh_param_names)}(params_matrix[idx, :]))
    end

    # Apply UH to each node
    node_sols = map(1:size(input, 2)) do i
        uh(input[:, i, :], extract_params_cv[i], config)
    end

    # Stack results efficiently (O(n) instead of O(n²))
    sol_mat = cat(node_sols..., dims=1)  # Splatting is more efficient than reduce
    reshape(sol_mat, 1, size(input)[2], size(input)[3])
end

# Export interfaces
export UnitHydrograph, @unithydro
