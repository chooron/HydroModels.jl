
"""
    UnitHydrograph{MS, UF, MF, HT, NT} <: AbstractHydrograph

Represents a Unit Hydrograph routing component for convolution.

$(FIELDS)
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

    function UnitHydrograph(
        inputs::AbstractVector, outputs::AbstractVector, params::AbstractVector,
        uh_func::Function, max_lag::Function;
        name::Optional{Symbol}=nothing, hru_types::Vector{Int}=Int[],
        kwargs...
    )
        infos = HydroInfos(
            inputs=!isempty(inputs) ? tosymbol.(inputs) : Symbol[],
            params=!isempty(params) ? tosymbol.(params) : Symbol[],
            outputs=!isempty(outputs) ? tosymbol.(outputs) : Symbol[]
        )
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name
        return new{length(hru_types) > 1,typeof(uh_func),typeof(max_lag_func),typeof(hru_types),typeof(infos)}(
            uh_name, uh_func, max_lag, hru_types, infos
        )
    end

    function UnitHydrograph(
        inputs::AbstractVector, outputs::AbstractVector, params::AbstractVector;
        uh_conds::AbstractVector{<:Pair}, name::Optional{Symbol}=nothing,
        kwargs...
    )
        max_lag = get(kwargs, :max_lag, uh_conds[1][1])
        hru_types = get(kwargs, :hru_types, Int[])
        param_names = !isempty(params) ? tosymbol.(Num.(params)) : Symbol[]
        uh_func, max_lag_func = build_uh_func(uh_conds, param_names, max_lag)
        @assert length(inputs) == length(outputs) == 1 "only one input and one output is supported"
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

A macro to simplify the construction of a `UnitHydrograph`.

# Usage
The macro takes an optional name and a `begin...end` block that defines the hydrograph.
In the `uh_func` block, `t` is a special variable representing time, and other variables are treated as parameters.

```julia
@variables P, Q, lag # Define symbolic variables
@unithydro :my_uh begin
    uh_func = begin
        # Piecewise UH definition: max_time => expression
        lag => 0.5 * (t / lag)^2.5
        2lag => 1.0 - 0.5 * abs(2 - t / lag)^2.5
    end
    uh_vars = P => Q # Defines input => output variables
end
```

# Arguments
- `name`: (Optional) A `Symbol` for the name of the `UnitHydrograph`.
- The `begin...end` block must contain:
    - `uh_func`: A block of `max_time => expression` pairs defining the piecewise unit hydrograph. The variable `t` can be used to represent time.
    - `uh_vars`: A `Pair` of `input_variable => output_variable`.
- Other assignments in the block are passed as keyword arguments to the `UnitHydrograph` constructor.

# Returns
- An instance of `UnitHydrograph`.
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
    uh_conds, cond_values = Pair[], []
    for uh_expr in uh_func_expr.args
        uh_expr isa LineNumberNode && continue
        if Meta.isexpr(uh_expr, :call) && uh_expr.args[1] == :(=>)
            push!(uh_conds, uh_expr.args[2] => uh_expr.args[3])
            push!(cond_values, uh_expr.args[3])
        end
    end

    return esc(quote
        params_val = reduce(union, map(
            val -> filter(x -> HydroModels.isparameter(x), HydroModels.get_variables(val)),
            [$(cond_values...)]
        ))
        UnitHydrograph(
            $uh_inputs_expr, $uh_outputs_expr, params_val;
            uh_conds=$uh_conds, max_lag=$(uh_conds[1][1]), name=$(name),
            $(kwargs_vec...)
        )
    end)
end


"""
$(TYPEDSIGNATURES)

Builds and returns a `uh_func` and a `max_lag_func` from piecewise conditions.

The generated `uh_func` computes weights based on time `t` and parameters `pas`.
The generated `max_lag_func` computes the maximum lag time from parameters `pas`.
Both functions are created using `@RuntimeGeneratedFunction` for performance.

Returns a `Tuple{Function, Function}` containing the generated `uh_func` and `max_lag_func`.
"""
function build_uh_func(uh_conds::AbstractVector{<:Pair}, params::AbstractVector{Symbol}, max_lag::Number)
    conditions_rev = vcat([0], reverse(first.(uh_conds)))
    values_rev = reverse(last.(uh_conds))
    params_assign_calls = generate_param_assignments(params=params)

    values_exprs = map(eachindex(values_rev)) do i
        :(
            if $(toexpr(conditions_rev[i])) ≤ t ≤ $(toexpr(conditions_rev[i+1]))
                return $(toexpr(values_rev[i]))
            end
        )
    end

    uh_func_expr = :(function (t, pas)
        $(params_assign_calls...)
        $(values_exprs...)
        return 1.0
    end)

    max_lag_expr = :(function (pas)
        $(params_assign_calls...)
        return ceil($(toexpr(max_lag)))
    end)

    return @RuntimeGeneratedFunction(uh_func_expr), @RuntimeGeneratedFunction(max_lag_expr)
end


"""
    (uh::UnitHydrograph)(input, pas; kwargs...)

Apply the unit hydrograph convolution to an input time series. This is the functor implementation for `UnitHydrograph`.

This method is dispatched based on the dimensionality of the `input` array and the `MS` type parameter of the struct, which indicates if it's a multi-node setup.

- **2D Input:** `(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...)`
  Performs convolution for a single location.
- **3D Input:** `(uh::UnitHydrograph{true})(input::AbstractArray{T,3}, pas::ComponentVector; kwargs...)`
  Performs convolution for multiple locations (nodes). It iterates over the nodes and applies the 2D method for each.

A method for 1D `AbstractVector` input is defined to throw an error, as single time points are not supported.
"""
function (uh::UnitHydrograph)(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T}
    solver = ManualSolver(mutable=true)
    timeidx = collect(1:size(input, 2))
    interp_func = DirectInterpolation(input, timeidx)
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

function (uh::UnitHydrograph{true})(input::AbstractArray{T,3}, pas::ComponentVector; kwargs...) where {T}
    ptyidx = uh.hru_types
    uh_param_names = get_param_names(uh)
    extract_params = eachrow(reshape(reduce(vcat, pas[:params][uh_param_names]), :, length(uh_param_names))[ptyidx, :])
    extract_params_cv = map(eachindex(ptyidx)) do idx
        ComponentVector(params=NamedTuple{Tuple(uh_param_names)}(extract_params[idx]))
    end
    node_sols = uh.(eachslice(input, dims=2), extract_params_cv)
    sol_mat = reduce((m1, m2) -> cat(m1, m2, dims=1), node_sols)
    return reshape(sol_mat, 1, size(input)[2], size(input)[3])
end
