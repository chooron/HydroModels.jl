"""
    UHFunction{uhtype} <: Function

Represents a unit hydrograph function for routing water through a hydrological system.

# Fields
- `uhtype::Symbol`: A symbol indicating the type of unit hydrograph function. Supported types are:
  - `:UH_1_HALF`: A triangular unit hydrograph with a single peak at lag time
  - `:UH_2_FULL`: A triangular unit hydrograph with a peak at 2*lag time
  - `:CUSTOM`: A custom unit hydrograph function defined by the user
- `func::Union{Nothing, Function}`: For `:CUSTOM` type, the function that computes the unit hydrograph value
- `max_lag::Union{Nothing, Any}`: For `:CUSTOM` type, the maximum lag time multiplier

# Constructors
- `UHFunction(uhtype::Symbol)`: Creates a predefined unit hydrograph function of the specified type
- `UHFunction(func::Function, max_lag)`: Creates a custom unit hydrograph function

# Methods
- `(uh::UHFunction{uhtype})(t, lag)`: Computes the unit hydrograph value at time `t` given the lag time `lag`
- `get_uh_tmax(uh::UHFunction{uhtype}, lag)`: Returns the maximum time required for computing the unit hydrograph with the given lag time `lag`

# Examples
```julia
# Create a predefined unit hydrograph function
uh_func = UHFunction(:UH_2_FULL)

# Create a custom unit hydrograph function
custom_uh = UHFunction((t, lag) -> t < lag ? (t/lag)^2 : 1.0, 1.0)
```

"""
struct UHFunction{uhtype}
    func::Union{Nothing,Function}
    max_lag::Union{Nothing,Any}

    function UHFunction(uhtype::Symbol)
        return new{uhtype}(nothing, nothing)
    end

    function UHFunction(func::Function, max_lag)
        return new{:CUSTOM}(func, max_lag)
    end
end

function (uh::UHFunction{:UH_1_HALF})(t, lag)
    if t - lag > 0
        typeof(lag)(1)
    else
        (t / lag)^2.5
    end
end

get_uh_tmax(::UHFunction{:UH_1_HALF}, lag) = ceil(lag)

function (uh::UHFunction{:UH_2_FULL})(t, lag)
    if t - lag * 2 > 0
        typeof(lag)(1)
    elseif t - lag > 0
        (1 - 0.5 * abs(2 - t / lag)^2.5)
    else
        (0.5 * abs(t / lag)^2.5)
    end
end


"""
    UnitHydrograph{solvetype} <: AbstractRouteFlux

Represents a unit hydrograph model for routing water through a hydrological system.

# Fields
- `inputs::Vector{Num}`: A vector of input variables (Num).
- `outputs::Vector{Num}`: A vector of output variables (Num).
- `params::Vector{Num}`: A vector of parameter variables (Num).
- `uhfunc::Function`: The unit hydrograph function.
- `meta::HydroMeta`: A named tuple containing information about inputs, outputs, parameters, and states.

# Constructor
    UnitHydrograph(input::Num, param::Num, uhfunc::Function; solvetype::Symbol=:unithydro1)

# Arguments
- `input::Num`: The input variable.
- `param::Num`: The parameter variable.
- `uhfunc::Function`: The unit hydrograph function.
- `solvetype::Symbol`: The solver type (default is `:unithydro1`).

# Description
UnitHydrograph represents a unit hydrograph flux model for routing water through a hydrological system.
It uses a unit hydrograph function to transform input flows into routed output flows.

The structure supports different solving methods, specified by the `solvetype` parameter.
Currently, it implements two solver types:
- `:unithydro1`: Uses a discrete problem approach to calculate the routed flow.
- `:unithydro2`: Uses a sparse matrix approach for more efficient computation, especially for longer time series.
The choice of solver type can affect both the performance and memory usage of the model.

This flux model is particularly useful in hydrological modeling for representing the
temporal distribution of runoff as it moves through a watershed. It can account for the
lag and attenuation of flow as water travels through the system.

The `uhfunc` field holds the unit hydrograph function, which typically takes a parameter
(like time) and returns weights that describe how an input is distributed over time in the output.

When called, the UnitHydrograph object applies the unit hydrograph to the input flow series,
effectively convolving the input with the unit hydrograph to produce the routed output flow.

This structure is designed to be flexible and can be integrated into larger hydrological models
to represent various routing processes in different parts of a water system.

"""
struct UnitHydrograph{N,ST} <: AbstractHydrograph
    "The unit hydrograph function"
    uh_func::Function
    "calculate max lag"
    max_lag_func::Function
    "A named tuple containing information about inputs, outputs, parameters, and states"
    infos::NamedTuple

    function UnitHydrograph(
        inputs::AbstractVector{T},
        params::AbstractVector{T};
        uh_pairs::AbstractVector{<:Pair},
        max_lag::Number=uh_func[1][1],
        outputs::AbstractVector{T}=T[],
        configs::NamedTuple=(solvetype=:DISCRETE, suffix=:_lag),
        name::Union{Symbol,Nothing}=nothing,
    ) where {T<:Num}
        #* Setup the name information of the hydroroutement
        input_names = tosymbol.(inputs)
        output_names = length(outputs) == 0 ? Symbol.(input_names, configs.suffix) : tosymbol.(outputs)
        outputs = map(name -> only(@variables $name), output_names)
        infos = (; inputs=inputs, outputs=outputs, params=params)

        uh_func, max_lag_func = build_uh_func(uh_pairs, params, max_lag)
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name

        @assert configs.solvetype in [:DISCRETE, :SPARSE] "solvetype must be one of [:DISCRETE, :SPARSE]"
        return new{uh_name,configs.solvetype}(uh_func, max_lag_func, infos)
    end
end

"""
    @unithydro name begin
        uh_func = begin
            2lag => (1 - 0.5 * abs(2 - t / lag)^2.5)
            lag => (0.5 * abs(t / lag)^2.5)
        end
        uh_vars = [q1, q2]
        configs = (solvetype=:DISCRETE, suffix=:_lag)
    end

Creates a UnitHydrograph with the specified name, unit hydrograph function, variables, and configuration.

# Arguments
- `name`: Symbol for the unit hydrograph name (optional)
- `uh_func`: Block defining the unit hydrograph function as key-value pairs
  - Keys represent upper time bounds (in terms of lag parameter)
  - Values are expressions for calculating weights at time t
  - The lower bound for each interval is either 0 (for the smallest key) or the next smaller key
- `uh_vars`: Array of variables to route through the unit hydrograph
- `configs`: Named tuple of configuration options
  - `solvetype`: Solution method (`:DISCRETE` or `:SPARSE`)
  - `suffix`: Suffix to append to variable names for output variables

# Example
```julia
uh = @unithydro :maxbas_uh begin
    uh_func = begin
        2lag => (1 - 0.5 * (2 - t / lag)^2.5)
        lag => (0.5 * (t / lag)^2.5)
    end
    uh_vars = [q]
    configs = (solvetype=:DISCRETE, suffix=:_lag)
end
```

This creates a unit hydrograph named `:maxbas_uh` that routes the variable `q` and produces `q_lag`.
The unit hydrograph function has two segments:
- For t between 0 and lag: weight = (0.5 * (t / lag)^2.5)
- For t between lag and 2lag: weight = (1 - 0.5 * (2 - t / lag)^2.5)
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
            Num.(filter(x -> isparameter(x), get_variables(val)))
        end),
        Expr(:vect, cond_values...)
    ))

    return esc(quote
        UnitHydrograph(
            $uh_vars_expr, $params_expr;
            uh_pairs=$(uh_pairs), max_lag=$(uh_pairs[1][1]),
            name=$(name), configs=$configs
        )
    end)
end

"""
    (flux::UnitHydrograph)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the unit hydrograph flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step (not supported, will throw an error).
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API.

# Returns
- For matrix input: A matrix where each column is the result of applying the unit hydrograph to the corresponding input column.
- For array input: A array of routed outputs, with dimensions (output_var_names, node_names, ts_len).

# Notes
- The behavior differs based on the `solvetype` specified during the `UnitHydrograph` construction:
  - `:unithydro1` uses a discrete problem solver approach.
  - `:unithydro2` uses a sparse matrix convolution approach.
- Vector input is not supported and will throw an error.
"""

(::UnitHydrograph)(::AbstractVector, ::ComponentVector; kwargs...) = @error "UnitHydrograph is not support for single timepoint"

function (flux::UnitHydrograph{N,:DISCRETE})(input::AbstractArray{T,2}, pas::ComponentVector; kwargs...) where {T,N}
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    interp = get(kwargs, :interp, DirectInterpolation)
    interp_func = interp(input, timeidx)
    #* prepare the initial states
    uh_weight = map(t -> flux.uh_func(t, pas), 1:flux.max_lag_func(pas))[1:end-1]
    println(uh_weight)
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
    println(uh_weight)
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

function (uh::UnitHydrograph)(input::AbstractArray{T,3}, params::AbstractArray; kwargs...) where {T}
    #* Extract the initial state of the parameters and routement in the pas variable
    ptyidx = get(kwargs, :ptyidx, 1:size(input, 2))
    params = Vector(params)
    extract_params = reshape(view(params, ptyidx), 1, :)
    node_sols = uh.(eachslice(input, dims=2), eachslice(extract_params, dims=2))
    sol_mat = reduce((m1, m2) -> cat(m1, m2, dims=1), node_sols)
    return reshape(sol_mat, 1, size(input)[2], size(input)[3])
end