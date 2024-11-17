"""
    HydroFlux

Represents a simple flux component in a hydrological model.

# Arguments
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `exprs::Vector{Num}`: A vector of expressions describing the formulas for output variables.
- `meta::HydroMeta`: Contains metadata about the flux, including input, output, and parameter names.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `params::Vector{Num}`: A vector of parameter variables.
- `exprs::Vector{Num}`: A vector of expressions describing the formulas for output variables.
- `func::Function`: A compiled function that calculates the flux.
- `meta::HydroMeta`: Contains metadata about the flux, including input, output, and parameter names.

# Constructors
    HydroFlux(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num}, meta::HydroMeta)
    HydroFlux(flux_names::Pair{Vector{Symbol},Vector{Symbol}}, param_names::Vector{Symbol}=Symbol[]; flux_funcs::Vector{<:Function}=Function[])

# Description
`HydroFlux` is a structure that encapsulates a simple flux calculation in a hydrological model. 
It can be constructed either by providing explicit inputs, outputs, parameters, and expressions, 
or by specifying names for fluxes and parameters along with optional flux functions.

The structure automatically compiles the provided expressions or functions into an efficient 
calculation function, which can be used to compute flux values given input and parameter values.

This structure is particularly useful for representing straightforward hydrological processes 
where the relationship between inputs and outputs can be expressed as simple mathematical formulas.
"""
struct HydroFlux <: AbstractHydroFlux
    "Vector of input variables"
    inputs::Vector{Num}
    "Vector of output variables"
    outputs::Vector{Num}
    "Vector of parameter variables"
    params::Vector{Num}
    "Vector of expressions describing the formulas for output variables"
    exprs::Vector{Num}
    "Compiled function that calculates the flux"
    func::Function
    "Metadata about the flux, including input, output, and parameter names"
    meta::HydroMeta

    function HydroFlux(
        inputs::Vector{Num},
        outputs::Vector{Num},
        params::Vector{Num};
        exprs::Vector{Num}=Num[],
    )
        #* name the flux
        output_names = Symbolics.tosymbol.(outputs, escape=false)
        flux_name = Symbol(Symbol(reduce((x, y) -> Symbol(x, y), output_names)), :_sflux)
        #* construct meta
        meta = HydroMeta(name=flux_name, inputs=inputs, outputs=outputs, params=params)
        #* if no expression provided, use the hydrology formula library to build the flux
        if length(exprs) == 0
            @info "No expression provided, using the hydrology formula library (`HydroModelLibrary.jl`) to build the flux"
            #* Get the corresponding calculation formula according to the input and output parameter names
            hydrofunc = HydroEquation(input_names, output_names, param_names)
            exprs = HydroModelLibrary.expr(hydrofunc)
        end
        #* build flux function
        flux_func = build_flux_func(inputs, outputs, params, exprs)

        return new(
            inputs,
            outputs,
            params,
            exprs,
            flux_func,
            meta
        )
    end

    #* construct hydro flux with input fluxes and output fluxes
    HydroFlux(fluxes::Pair{Vector{Num},Vector{Num}}, params::Vector{Num}=Num[]; exprs::Vector{Num}=Num[]) = HydroFlux(fluxes[1], fluxes[2], params, exprs=exprs)
end

"""
    (flux::AbstractHydroFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the simple flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API
    - `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).

# Returns
- For vector input: The result of applying the flux function to the input and parameters.
- For matrix input: A matrix where each column is the result of applying the flux function to the corresponding input column.
- For 3D array input: A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).
"""
function (flux::AbstractHydroFlux)(input::Vector, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, 1)
    params_vec = collect([pas[:params][nm] for nm in get_param_names(flux)])
    flux.func(input, params_vec, timeidx)
end

function (flux::AbstractHydroFlux)(input::Matrix, pas::ComponentVector; kwargs...)
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    # assert the input params must include all the parameters in the flux
    @assert length(timeidx) == size(input, 2) "Time index length does not match the number of time steps"
    @assert all(nm in keys(pas[:params]) for nm in get_param_names(flux)) "Input parameters do not match the flux parameters, the flux parameters should be: $(get_param_names(flux))"
    params_vec = collect([pas[:params][nm] for nm in get_param_names(flux)])
    reduce(hcat, flux.func.(eachslice(input, dims=2), Ref(params_vec), timeidx))
end

function (flux::AbstractHydroFlux)(input::Array, pas::ComponentVector; kwargs...)
    #* get kwargs
    ptypes = get(kwargs, :ptypes, collect(keys(pas[:params])))
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 3)))
    @assert length(timeidx) == size(input, 3) "Time index length does not match the number of time steps"

    #* extract params and nn params
    params_collect = [pas[:params][ptype] for ptype in ptypes]
    #* check params input is correct
    for (ptype, params_item) in zip(ptypes, params_collect)
        @assert all(param_name in keys(params_item) for param_name in get_param_names(flux)) "Missing required parameters. Expected all of $(get_param_names(flux)), but got $(keys(params_item)) at param type: $ptype."
    end
    params_vec = collect([collect([params_item[pname] for pname in get_param_names(flux)]) for params_item in params_collect])

    #* array dims: (var_names * node_names * ts_len)
    flux_output_vec = [reduce(hcat, flux.func.(eachslice(input[:, i, :], dims=2), Ref(params_vec[i]), timeidx)) for i in eachindex(ptypes)]
    flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    permutedims(flux_output_arr, (3, 1, 2))
end


"""
    NeuralFlux <: AbstractNeuralFlux

Represents a neural network-based flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `outputs::Vector{Num}`: A vector of output variables.
- `nnparam::Symbolics.Arr`: An array of neural network parameters.
- `expr::Symbolics.Arr{Num,1}`: Expressions describing the formulas for output variables.
- `func::Function`: A compiled function that calculates the flux using the neural network.
- `meta::HydroMeta`: Contains metadata about the flux, including input, output, and neural network parameter names.
- `nninfos::NamedTuple`: Contains information about the neural network's input and output structure.

# Constructors
    # 1. Construct a NeuralFlux with specified input/output fluxes and a neural network, the neural network should specify the name
    NeuralFlux(fluxes::Pair{Vector{Num},Vector{Num}}, chain::Lux.AbstractExplicitContainerLayer)

# Description
`NeuralFlux` is a structure that encapsulates a neural network-based flux calculation in a hydrological model. 
It combines symbolic mathematics with neural networks to create a flexible and powerful representation of complex hydrological processes.

The structure automatically handles the integration of the neural network into the symbolic framework, 
allowing for seamless incorporation of machine learning models into traditional hydrological equations.

This structure is particularly useful for representing complex, non-linear relationships in hydrological systems 
where traditional equations may be insufficient or unknown.
"""
struct NeuralFlux <: AbstractNeuralFlux
    "Vector of input variables"
    inputs::Vector{Num}
    "Vector of output variables"
    outputs::Vector{Num}
    "Array of neural network parameters"
    nnparam::Symbolics.Arr
    "Array of expressions describing the formulas for output variables"
    expr::Symbolics.Arr{Num,1}
    "Compiled function that calculates the flux using the neural network"
    func::Function
    "Metadata about the flux, including input, output, and neural network parameter names"
    meta::HydroMeta
    "Information about the neural network's input and output structure"
    nninfos::NamedTuple

    function NeuralFlux(
        inputs::Vector{Num},
        outputs::Vector{Num},
        chain, # ::LuxCore.AbstractExplicitContainerLayer
    )
        #* assert the chain has a name
        @assert chain.name isa Symbol "Neural network chain should have a name with Symbol type"
        #* Get the neural network name (neural flux param name) and object
        chain_name = chain.name
        #* Initialize model parameter type for model parameter dimension definition
        init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
        params_axes = getaxes(init_params)

        #* Define parameter variables according to initialization parameters: Define type as Vector{parameter length}
        chain_params = first(@parameters $chain_name[1:length(init_params)] = Vector(init_params))
        #* Use Symbolics.array_term to define the slow-building parameter variables:
        #* when the model is called, it is rebuilt into the ComponentVector type according to
        #* the axes of `init_params` and the Vector type of the parameter as the calculation parameter input
        lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, params_axes, size=size(chain_params))

        #* Constructing neural network input and output variables
        #* The input and output of the model can only be of type Symbolics.Arr{Num, 1},
        #* so it cannot be defined based on outputs and output_vars
        nn_input_name = Symbol(chain_name, :_input)
        nn_output_name = Symbol(chain_name, :_output)
        nn_input = first(@variables $(nn_input_name)[1:length(inputs)])
        nn_output = first(@variables $(nn_output_name)[1:length(outputs)])

        #* Constructing a calculation expression based on a neural network
        flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)
        nn_func = (x, p) -> LuxCore.stateless_apply(chain, x, ComponentVector(p, params_axes))

        #* neuralflux meta
        meta = HydroMeta(name=Symbol(chain_name, :_nflux), inputs=inputs, outputs=outputs, params=params, nns=[chain_name])
        nninfos = (inputs=nn_input, outputs=nn_output, paramlen=length(init_params))

        new(
            inputs, outputs, chain_params,
            flux_expr, nn_func,
            meta, nninfos,
        )
    end

    #* construct neural flux with input fluxes and output fluxes
    NeuralFlux(fluxes::Pair{Vector{Num},Vector{Num}}, chain) = NeuralFlux(fluxes[1], fluxes[2], chain)
end

"""
    (flux::AbstractFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the flux model (simple or neural) to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API

# Returns
- For vector input: The result of applying the flux function to the input and parameters.
- For matrix input: A matrix where each column is the result of applying the flux function to the corresponding input column.
- For 3D array input: A 3D array of flux outputs, with dimensions (output_var_names, node_names, ts_len).

# Note
For neural flux models, the parameters are accessed from `pas[:nn]` instead of `pas[:params]`.
"""
function (flux::AbstractNeuralFlux)(input::Vector, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][get_nn_names(flux)[1]]
    flux.func(input, nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::Matrix, pas::ComponentVector; kwargs...)
    nn_params_vec = pas[:nn][get_nn_names(flux)[1]]
    flux.func(input', nn_params_vec)
end

function (flux::AbstractNeuralFlux)(input::Array, pas::ComponentVector; kwargs...)
    nn_params = pas[:nn][get_nn_names(flux)[1]]
    #* array dims: (ts_len * node_names * var_names)
    flux_output_vec = [flux.func(input[:, i, :], nn_params) for i in 1:size(input)[2]]
    flux_output_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), flux_output_vec)
    permutedims(flux_output_arr, (1, 3, 2))
end

"""
    StateFlux <: AbstractStateFlux

Represents a state flux component in a hydrological model.

# Fields
- `inputs::Vector{Num}`: A vector of input variables.
- `state::Num`: The state variable.
- `params::Vector{Num}`: A vector of parameter variables.
- `expr::Num`: The expression describing the state variable's formula.
- `func::Function`: A function to calculate the state flux.
- `meta::HydroMeta`: Contains metadata about inputs, state, parameters, and neural networks.

# Constructors
    # 1. Detailed specification of inputs, state, parameters, and state expression
    StateFlux(fluxes::Vector{Num}, state::Num, params::Vector{Num}=Num[]; expr::Num)
    # 2. Automatic construction of state expression as the difference between sum of input fluxes and sum of output fluxes
    StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num)

# Description
StateFlux is a structure that represents a state flux in a hydrological model. It encapsulates 
the relationship between input fluxes, output fluxes, and a state variable. The structure 
provides methods to calculate state changes based on the provided expressions and parameters.

The first constructor allows for detailed specification of inputs, state, parameters, and the 
state expression. The second constructor automatically constructs a state expression as the 
difference between sum of input fluxes and sum of output fluxes. The third constructor is used 
for simple state transitions.

This structure is particularly useful in building complex hydrological models where state 
variables evolve over time based on various input and output fluxes.
"""
struct StateFlux <: AbstractStateFlux
    "A map of input names (Symbol) and its variables (Num)"
    inputs::Vector{Num}
    "A map of state names (Symbol) and its variables (Num)"
    state::Num
    "A map of parameters names (Symbol) and its variables (Num)"
    params::Vector{Num}
    "flux expressions to descripe the formula of the state variable"
    expr::Num
    "flux expressions to descripe the formula of the output variable"
    func::Function
    "bucket information: keys contains: input, output, param, state"
    meta::HydroMeta

    function StateFlux(
        fluxes::Vector{Num}, state::Num, params::Vector{Num}=Num[]; expr::Num
    )
        #* Convert to a symbol based on the variable
        state_name = Symbolics.tosymbol(state, escape=false)
        meta = HydroMeta(name=Symbol(state_name, :_stflux), inputs=fluxes, states=[state], params=params)
        state_func = build_flux_func(fluxes, [state], params, [expr])
        return new(
            fluxes,
            state,
            params,
            expr,
            state_func,
            meta
        )
    end
    #* construct state flux with input fluxes and output fluxes
    StateFlux(fluxes::Pair{Vector{Num},Vector{Num}}, state::Num) = StateFlux(vcat(fluxes[1], fluxes[2]), state, expr=sum(fluxes[1]) - sum(fluxes[2]))
    #* construct state flux with state variables
    StateFlux(states::Pair{Num,Num}) = StateFlux([states[2]], states[1], expr=states[2] - states[1])
end

"""
    (flux::AbstractStateFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

Apply the state flux model to input data of various dimensions.

# Arguments
- `input`: Input data, which can be:
  - `Vector`: A vector of input values for a single time step.
  - `Matrix`: A matrix of input values, where each column represents a different time step.
  - `Array`: A 3D array of input values, with dimensions (var_names, node_names, ts_len).
- `pas::ComponentVector`: A component vector containing parameter values.
- `ptypes::AbstractVector{Symbol}`: A vector of symbols representing parameter categories (only used for `Array` input).
- `kwargs...`: Additional keyword arguments (unused in this function), provided for compatibility with the component callable function API.

# Returns
This function does not actually return a value, as state flux models cannot be run directly.

# Notes
- State flux models cannot be run directly and will throw an error if attempted.
- To use state flux models, they should be incorporated into a HydroFlux or other composite flux model.
"""
(::AbstractStateFlux)(::Vector, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::Matrix, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"
(::AbstractStateFlux)(::Array, ::ComponentVector; kwargs...) = @error "State Flux cannot run directly, please using HydroFlux to run"

"""
    UnitHydroFlux{solvetype} <: AbstractRouteFlux

Represents a unit hydrograph flux model for routing water through a hydrological system.

# Fields
- `inputs::Vector{Num}`: A vector of input variables (Num).
- `outputs::Vector{Num}`: A vector of output variables (Num).
- `params::Vector{Num}`: A vector of parameter variables (Num).
- `uhfunc::Function`: The unit hydrograph function.
- `meta::HydroMeta`: A named tuple containing information about inputs, outputs, parameters, and states.

# Constructor
    UnitHydroFlux(input::Num, param::Num, uhfunc::Function; solvetype::Symbol=:unithydro1)

# Arguments
- `input::Num`: The input variable.
- `param::Num`: The parameter variable.
- `uhfunc::Function`: The unit hydrograph function.
- `solvetype::Symbol`: The solver type (default is `:unithydro1`).

# Description
UnitHydroFlux represents a unit hydrograph flux model for routing water through a hydrological system. 
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

When called, the UnitHydroFlux object applies the unit hydrograph to the input flow series, 
effectively convolving the input with the unit hydrograph to produce the routed output flow.

This structure is designed to be flexible and can be integrated into larger hydrological models 
to represent various routing processes in different parts of a water system.

"""
struct UnitHydroFlux{solvetype} <: AbstractUnitHydroFlux
    "A vector of input variables (Num)"
    inputs::Vector{Num}
    "A vector of output variables (Num)"
    outputs::Vector{Num}
    "A vector of parameter variables (Num)"
    params::Vector{Num}
    "The unit hydrograph function"
    uhfunc::UHFunction
    "A named tuple containing information about inputs, outputs, parameters, and states"
    meta::HydroMeta

    function UnitHydroFlux(
        input::Num,
        param::Num;
        uhtype::Symbol=:UH_1_HALF,
        output::Union{Num,Nothing}=nothing,
        solvetype::Symbol=:DISCRETE,
    )
        input_name = Symbolics.tosymbol(input, escape=false)
        param_name = Symbolics.tosymbol(param, escape=false)
        if isnothing(output)
            output_name = Symbol(input_name, :_routed)
            output = first(@variables $output_name)
        else
            output_name = Symbolics.tosymbol(output, escape=false)
        end

        uhfunc = UHFunction(uhtype)
        #* Setup the name information of the hydroroutement
        meta = HydroMeta(inputs=[input_name], outputs=[output_name], params=[param_name], name=Symbol(output_name, :_uh_flux))

        return new{solvetype}(
            [input],
            [output],
            [param],
            uhfunc,
            meta
        )
    end
end

"""
    (flux::UnitHydroFlux)(input::Union{Vector,Matrix,Array}, pas::ComponentVector; ptypes::AbstractVector{Symbol}=Symbol[], kwargs...)

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
- The behavior differs based on the `solvetype` specified during the `UnitHydroFlux` construction:
  - `:unithydro1` uses a discrete problem solver approach.
  - `:unithydro2` uses a sparse matrix convolution approach.
- Vector input is not supported and will throw an error.
"""

(::UnitHydroFlux)(::Vector, ::ComponentVector; kwargs...) = @error "UnitHydroFlux is not support for single timepoint"

function (flux::UnitHydroFlux{:DISCRETE})(input::Matrix, pas::ComponentVector; kwargs...)
    solver = get(kwargs, :solver, DiscreteSolver())
    timeidx = get(kwargs, :timeidx, collect(1:size(input, 2)))
    input_vec = input[1, :]
    #* convert the lagflux to a discrete problem
    lag_prob!(du, u, p, t) = begin
        du[:] = input_vec[Int(t)] .* p[:weight] .+ [diff(u); -u[end]]
        nothing
    end
    #* prepare the initial states
    lag = pas[:params][get_param_names(flux)[1]]
    uh_weight = map(t -> flux.uhfunc(t, lag), 1:get_uh_tmax(flux.uhfunc, lag))
    initstates = input_vec[1] .* uh_weight ./ sum(uh_weight)
    #* solve the problem
    sol = solver(lag_prob!, ComponentVector(weight=uh_weight ./ sum(uh_weight)), initstates, timeidx)
    reshape(sol[1, :], 1, length(input_vec))
end

function (flux::UnitHydroFlux{:SPARSE})(input::Matrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    lag = pas[:params][get_param_names(flux)[1]]
    uh_weight = map(t -> flux.uhfunc(t, lag), 1:get_uh_tmax(flux.uhfunc, lag))
    #* the weight of the unit hydrograph is normalized by the sum of the weights
    uh_result = [-(i - 1) => uh_wi .* input_vec ./ sum(uh_weight) for (i, uh_wi) in enumerate(uh_weight)]
    #* construct the sparse matrix
    uh_sparse_matrix = spdiagm(uh_result...)
    #* sum the matrix
    sum_route = sum(uh_sparse_matrix, dims=2)[1:end-length(uh_weight)+1]
    reshape(sum_route, 1, length(input_vec))
end

# todo: 卷积计算的结果与前两个计算结果不太一致
function (flux::UnitHydroFlux{:INTEGRAL})(input::Matrix, pas::ComponentVector; kwargs...)
    input_vec = input[1, :]
    itp_method = get(kwargs, :interp, LinearInterpolation)
    itp = itp_method(input_vec, collect(1:length(input_vec)), extrapolate=true)
    #* construct the unit hydrograph function based on the interpolation method and parameter
    lag = pas[:params][get_param_names(flux)[1]]
    tmax = get_uh_tmax(flux.uhfunc, lag)
    uh_sum = solve(IntegralProblem(flux.uhfunc, (0, tmax), lag), QuadGKJL()).u
    uh_itg_func = (x, p) -> flux.uhfunc(x, lag) * itp(p - x) / uh_sum
    #* solve the integral problem
    prob = IntegralProblem(uh_itg_func, (0, tmax), 1.0)
    routed_result = map(1:length(input_vec)) do t
        prob = remake(prob, p=t)
        sol = solve(prob, QuadGKJL())
        sol.u
    end
    reshape(routed_result, 1, length(input_vec))
end

function (uh::AbstractUnitHydroFlux)(input::Array, pas::ComponentVector; kwargs...)
    #* array dims: (variable dim, num of node, sequence length)
    #* Extract the initial state of the parameters and routement in the pas variable
    ptypes = get(kwargs, :ptypes, collect(keys(pas[:params])))
    pytype_params = [pas[:params][ptype] for ptype in ptypes]

    sols = map(eachindex(ptypes)) do (idx)
        tmp_pas = ComponentVector(params=pytype_params[idx])
        node_sols = reduce(hcat, uh(input[:, idx, :], tmp_pas))
        node_sols
    end
    sol_arr = reduce((m1, m2) -> cat(m1, m2, dims=3), sols)
    return permutedims(sol_arr, (1, 3, 2))
end