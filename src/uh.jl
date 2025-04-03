"""
    UHFunction{uhtype} <: Function

Represents a unit hydrograph function for routing water through a hydrological system.

# Fields
- `uhtype::Symbol`: A symbol indicating the type of unit hydrograph function. Currently, only `:UH_1_HALF` and `:UH_2_FULL` are supported.

# Methods
- `(uh::UHFunction{uhtype})(t, lag)`: Computes the unit hydrograph value at time `t` given the lag time `lag`.
- `get_uh_tmax(uh::UHFunction{uhtype}, lag)`: Returns the maximum time required for computing the unit hydrograph with the given lag time `lag`.

"""
struct UHFunction{uhtype}
    function UHFunction(uhtype::Symbol)
        return new{uhtype}()
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

get_uh_tmax(::UHFunction{:UH_2_FULL}, lag) = 2 * ceil(lag)


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
    uhfunc::UHFunction
    "A named tuple containing information about inputs, outputs, parameters, and states"
    infos::NamedTuple

    function UnitHydrograph(
        input::T,
        output::T,
        params::Vector{T};
        uhfunc::UF,
        name::Union{Symbol,Nothing}=nothing,
        solvetype::Symbol=:SPARSE,
    ) where {T<:Num,UF<:UHFunction}
        @assert solvetype in [:DISCRETE, :SPARSE] "solvetype must be one of [:DISCRETE, :SPARSE]"
        #* Setup the name information of the hydroroutement
        infos = (;inputs = [input], outputs = [output], params = params)
        uh_name = isnothing(name) ? Symbol("##uh#", hash(infos)) : name
        return new{uh_name, solvetype}(uhfunc, infos)
    end

    function UnitHydrograph(
        fluxes::Pair{Vector{T},Vector{T}},
        params::Vector{T};
        uhfunc::UF,
        solvetype::Symbol=:SPARSE,
    ) where {T<:Num,UF<:UHFunction}
        input = fluxes[1][1]
        output = fluxes[2][1]
        return UnitHydrograph(input, output, params, uhfunc=uhfunc, solvetype=solvetype)
    end
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

function (flux::UnitHydrograph{N,:DISCRETE})(input::AbstractArray{T,2}, params::AbstractVector; kwargs...) where {T,N}
    solver = get(kwargs, :solver, ManualSolver(mutable=true))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    input_vec = input[1, :]
    #* convert the lagflux to a discrete problem
    lag_du_func(u, p, t) = input_vec[Int(t)] .* p[:weight] .+ [diff(u); -u[end]]
    #* prepare the initial states
    lag = Vector(params)[1]
    uh_weight = map(t -> flux.uhfunc(t, lag), 1:get_uh_tmax(flux.uhfunc, lag))[1:end-1]
    if length(uh_weight) == 0
        @warn "The unit hydrograph weight is empty, please check the unit hydrograph function"
        return input
    else
        #* solve the problem
        sol = solver(lag_du_func, ComponentVector(weight=uh_weight ./ sum(uh_weight)), zeros(length(uh_weight)), timeidx)
        reshape(sol[1, :], 1, length(input_vec))
    end
end

function (flux::UnitHydrograph{N,:SPARSE})(input::AbstractArray{T,2}, params::AbstractVector; kwargs...) where {T,N}
    input_vec = input[1, :]
    lag = Vector(params)[1]
    uh_weight = map(t -> flux.uhfunc(t, lag), 1:get_uh_tmax(flux.uhfunc, lag))[1:end-1]

    if length(uh_weight) == 0
        @warn "The unit hydrograph weight is empty, please check the unit hydrograph function"
        return input
    else
        #* the weight of the unit hydrograph is normalized by the sum of the weights
        uh_result = [-(i - 1) => uh_wi .* input_vec ./ sum(uh_weight) for (i, uh_wi) in enumerate(uh_weight)]
        #* construct the sparse matrix
        uh_sparse_matrix = spdiagm(uh_result...)
        #* sum the matrix
        sum_route = sum(uh_sparse_matrix, dims=2)[1:end-length(uh_weight)+1]
        reshape(sum_route, 1, length(input_vec))
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