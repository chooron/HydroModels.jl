"""
Utility functions module - provides sorting, parameter expansion, variable extraction,
and parameter conversion utilities.
"""

# ============================================================================
# Parameter conversion utilities
# ============================================================================

"""
    _as_componentvector(p::ComponentVector)

Identity for ComponentVector - return as-is.
"""
@inline _as_componentvector(p::ComponentVector) = p

"""
    _as_componentvector(p::AbstractVector, axes)

Convert a plain vector to ComponentVector using stored axes.
"""
@inline _as_componentvector(p::AbstractVector, axes) = ComponentVector(p, axes)

"""
    sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})

Topologically sort flux components based on their variable dependencies.

This function builds a directed graph representing dependencies between flux components
(where edges connect input variables to output variables) and performs a topological sort
to determine the correct calculation order.

# Examples
```jldoctest
julia> fluxes = [flux1, flux2, flux3];
julia> sorted_fluxes = sort_fluxes(fluxes);
# sorted_fluxes now contains components in correct calculation order
```
"""
function sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})
    # Extract all input and output names
    input_names = reduce(union, get_input_names.(fluxes); init=Symbol[])
    output_names = reduce(union, get_output_names.(fluxes); init=Symbol[])
    
    # Filter: exclude outputs from inputs, exclude inputs from outputs
    input_names = setdiff(input_names, output_names)
    output_names = setdiff(output_names, input_names)
    
    # Build mapping from output names to flux instances (with size hint for performance)
    num_outputs = sum(length ∘ get_output_names, fluxes; init=0)
    fluxes_map = Dict{Symbol,AbstractComponent}()
    sizehint!(fluxes_map, num_outputs)
    for flux in fluxes
        for name in get_output_names(flux)
            fluxes_map[name] = flux
        end
    end
    
    # Build directed computation graph
    var_names = vcat(input_names, output_names)
    var_indices = Dict(nm => idx for (idx, nm) in enumerate(var_names))
    digraph = SimpleDiGraph(length(var_names))
    
    for flux in fluxes
        inp_names = get_input_names(flux)
        out_names = get_output_names(flux)
        for inp_nm in inp_names
            haskey(var_indices, inp_nm) || continue
            for out_nm in out_names
                haskey(var_indices, out_nm) || continue
                add_edge!(digraph, var_indices[inp_nm], var_indices[out_nm])
            end
        end
    end
    
    # Extract fluxes based on topological order
    sorted_fluxes = AbstractComponent[]
    seen_fluxes = Set{AbstractComponent}()
    
    for idx in topological_sort(digraph)
        var_name = var_names[idx]
        if var_name in output_names && haskey(fluxes_map, var_name)
            flux = fluxes_map[var_name]
            if flux ∉ seen_fluxes
                push!(sorted_fluxes, flux)
                push!(seen_fluxes, flux)
            end
        end
    end
    
    sorted_fluxes
end

"""
    sort_components(components::AbstractVector{<:AbstractComponent})

Topologically sort hydrological components based on their variable dependencies.

This function builds a dependency graph between components (inputs to outputs/states)
and uses topological sort to find the correct execution order.

# Examples
```jldoctest
julia> components = [bucket1, flux1, route1];
julia> sorted_components = sort_components(components);
# sorted_components now contains components in correct calculation order
```
"""
function sort_components(components::AbstractVector{<:AbstractComponent})
    input_names, output_names, state_names = get_var_names(components)
    
    # Build mapping from output/state names to components (with size hint)
    total_vars = length(output_names) + length(state_names)
    component_map = Dict{Symbol,AbstractComponent}()
    sizehint!(component_map, total_vars)
    for comp in components
        _, out_names, st_names = get_var_names(comp)
        for name in vcat(out_names, st_names)
            component_map[name] = comp
        end
    end
    
    # Build directed graph
    var_names = reduce(union, [input_names, output_names, state_names])
    var_indices = Dict(nm => idx for (idx, nm) in enumerate(var_names))
    digraph = SimpleDiGraph(length(var_names))
    
    for comp in components
        inp_names, out_names, st_names = get_var_names(comp)
        out_st_names = vcat(out_names, st_names)
        
        for inp_nm in inp_names
            haskey(var_indices, inp_nm) || continue
            for out_nm in out_st_names
                haskey(var_indices, out_nm) || continue
                add_edge!(digraph, var_indices[inp_nm], var_indices[out_nm])
            end
        end
    end
    
    # Extract sorted components
    sorted_components = AbstractComponent[]
    seen_components = Set{AbstractComponent}()
    
    for idx in topological_sort(digraph)
        var_name = var_names[idx]
        if haskey(component_map, var_name)
            comp = component_map[var_name]
            if comp ∉ seen_components
                push!(sorted_components, comp)
                push!(seen_components, comp)
            end
        end
    end
    
    sorted_components
end

"""
    expand_component_params(params::ComponentVector, param_names::AbstractVector, htypes::AbstractVector)

Expand specified parameters of a ComponentVector using a given index.

This is useful for spatial modeling, where a parameter vector needs to be expanded to match the grid size.

# Arguments
- `params`: ComponentVector containing parameters
- `param_names`: Parameter names to expand
- `htypes`: HRU type index vector

# Returns
- Expanded ComponentVector
"""
function expand_component_params(
    params::ComponentVector,
    param_names::AbstractVector{Symbol},
    htypes::AbstractVector{Int}
)
    # If no parameters need expansion, return directly
    isempty(param_names) && return params

    # Extract and expand parameters with bounds checking
    params_nt = NamedTuple(params[:params])
    expanded_params = map(param_names) do pname
        if !haskey(params_nt, pname)
            throw(ArgumentError("Parameter $pname does not exist in params"))
        end
        param_vec = params_nt[pname]
        # Validate indices
        if !isempty(htypes) && (minimum(htypes) < 1 || maximum(htypes) > length(param_vec))
            throw(BoundsError(param_vec, htypes))
        end
        param_vec[htypes]
    end

    # Build new parameter NamedTuple
    new_params = (; params=NamedTuple{Tuple(param_names)}(expanded_params))

    # Merge and return
    merge(params |> NamedTuple, new_params) |> ComponentVector
end

"""
    get_default_states(component::AbstractComponent, input::AbstractArray{T,2}) where {T}

Get a ComponentVector of default (zeroed) states for a given component.

The shape of the states depends on the input array dimensions:
- 2D input `(time, features)`: Returns a vector of zeros for each state.
- 3D input `(time, locations, features)`: Returns a NamedTuple of state vectors for each location.

# Arguments
- `component`: Hydrological component
- `input`: Input array

# Returns
- ComponentVector of default states
"""
function get_default_states(component::AbstractComponent, input::AbstractArray{T,2}) where {T}
    state_names = get_state_names(component)
    isempty(state_names) && return ComponentVector{T}()
    
    states_axes = Axis(NamedTuple{Tuple(state_names)}(1:length(state_names)))
    ComponentVector(zeros(T, length(state_names)), states_axes)
end

function get_default_states(component::AbstractComponent, input::AbstractArray{T,3}) where {T}
    state_names = get_state_names(component)
    isempty(state_names) && return ComponentVector{T}()
    
    num_nodes = size(input, 2)
    state_vectors = ntuple(_ -> zeros(T, num_nodes), length(state_names))
    ComponentVector(NamedTuple{Tuple(state_names)}(state_vectors))
end

"""
    extract_variables(expr)

Extract all variable symbols from a Julia expression.

# Arguments
- `expr`: Julia expression

# Returns
- Set containing all variable symbols
"""
function extract_variables(expr)
    vars = Set{Symbol}()
    _extract_vars!(expr, vars)
    return vars
end

# Helper function: recursively extract variables
_extract_vars!(::Any, ::Set{Symbol}) = nothing

function _extract_vars!(s::Symbol, vars::Set{Symbol})
    push!(vars, s)
end

function _extract_vars!(ex::Expr, vars::Set{Symbol})
    if ex.head == :call
        # Function call: extract arguments
        for arg in ex.args[2:end]
            _extract_vars!(arg, vars)
        end
    elseif ex.head == :.
        # Field access: only extract base variable
        _extract_vars!(ex.args[1], vars)
    elseif ex.head in (:(=), :function, :macro, :->, :kw)
        # Skip left side of assignments, function definitions, etc.
        return nothing
    else
        # Other expressions: recursively process all arguments
        for arg in ex.args
            _extract_vars!(arg, vars)
        end
    end
end

"""
    @hydroflux_for [name] for ... end

Macro: Create multiple HydroFlux equations through loop unrolling.

# Examples
```jldoctest
julia> @variables P[1:3], Q[1:3], k[1:3];
julia> @hydroflux_for :multi_flux for i in 1:3
           Q[i] ~ k[i] * P[i]
       end
```
"""
macro hydroflux_for(args...)
    name = length(args) == 1 ? nothing : args[1]
    loop_expr = length(args) == 1 ? args[1] : args[2]
    
    if !Meta.isexpr(loop_expr, :for)
        return :(error("@hydroflux_for macro must be used with a 'for' loop expression"))
    end
    
    # Parse loop variable and range
    loop_var = loop_expr.args[1].args[1]
    range_expr = loop_expr.args[1].args[2]
    
    # Evaluate range
    range_val = if Meta.isexpr(range_expr, :call) && range_expr.args[1] == :(:)
        if length(range_expr.args) == 3
            range_expr.args[2]:range_expr.args[3]
        elseif length(range_expr.args) == 4
            range_expr.args[2]:range_expr.args[3]:range_expr.args[4]
        else
            eval(range_expr)
        end
    else
        eval(range_expr)
    end
    
    # Extract loop body
    loop_body = loop_expr.args[2]
    if !Meta.isexpr(loop_body, :block)
        loop_body = Expr(:block, loop_body)
    end
    
    # Unroll loop
    equations = filter(x -> !(x isa LineNumberNode) && !Meta.isexpr(x, :line), loop_body.args)
    all_equations = []
    
    for i_val in range_val
        for eq in equations
            new_eq = deepcopy(eq)
            replace_loop_var!(new_eq, loop_var, i_val)
            push!(all_equations, new_eq)
        end
    end
    
    vect_eqs_expr = Expr(:tuple, all_equations...)
    
    # Check variable definitions
    for var_name in extract_variables(vect_eqs_expr)
        if !@isdefined(var_name)
            expr_str = string(vect_eqs_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end
    
    return esc(:(HydroFlux(exprs=$vect_eqs_expr, name=$name)))
end

# Helper function: recursively replace loop variable
function replace_loop_var!(expr, var_name, value)
    if expr isa Expr
        for i in 1:length(expr.args)
            if expr.args[i] == var_name
                expr.args[i] = value
            else
                replace_loop_var!(expr.args[i], var_name, value)
            end
        end
    end
end

"""
    get_initial_params(component::AbstractComponent; rng=Random.default_rng(), eltype=Float64)

Get initial parameters for a component based on its infos metadata.

Returns a ComponentVector with:
- `params`: NamedTuple of traditional parameters (initialized to 1.0)
- `nns`: NamedTuple of neural network parameters (randomly initialized)

Uses multiple dispatch for different component types.

# Examples
```julia
# For any component with params and/or nns
flux = @hydroflux q ~ k * p
params = get_initial_params(flux)
# Returns: ComponentVector(params=(k=1.0,))

# For component with neural networks
neural_flux = @neuralflux et ~ nn([prcp, temp])
params = get_initial_params(neural_flux)
# Returns: ComponentVector(params=(), nns=(nn_name=ComponentVector(...),))
```
"""
function get_initial_params(component::AbstractComponent; rng=Random.default_rng(), eltype=Float64)
    param_names = get_param_names(component)
    nn_names = get_nn_names(component)

    # Initialize traditional parameters to 1.0
    params_nt = if !isempty(param_names)
        NamedTuple{Tuple(param_names)}(ones(eltype, length(param_names)))
    else
        NamedTuple()
    end

    # Initialize neural network parameters
    nns_nt = if !isempty(nn_names)
        _get_nn_params(component, nn_names, rng)
    else
        NamedTuple()
    end

    # Build ComponentVector
    if !isempty(param_names) && !isempty(nn_names)
        return ComponentVector(params=params_nt, nns=nns_nt)
    elseif !isempty(param_names)
        return ComponentVector(params=params_nt)
    elseif !isempty(nn_names)
        return ComponentVector(nns=nns_nt)
    else
        return ComponentVector()
    end
end

# Helper function to extract neural network parameters
function _get_nn_params(component::AbstractNeuralFlux, nn_names, rng)
    if !isnothing(component.chain)
        ps = LuxCore.initialparameters(rng, component.chain)
        return NamedTuple{Tuple(nn_names)}((ComponentVector(ps),))
    end
    return NamedTuple()
end

function _get_nn_params(component::AbstractHydroBucket, nn_names, rng)
    if hasproperty(component, :neural_fluxes) && !isempty(component.neural_fluxes)
        all_nn_params = Dict{Symbol, ComponentVector}()
        for nf in component.neural_fluxes
            nf_nn_names = get_nn_names(nf)
            nf_params = _get_nn_params(nf, nf_nn_names, rng)
            merge!(all_nn_params, pairs(nf_params))
        end
        return NamedTuple(all_nn_params)
    end
    return NamedTuple()
end

function _get_nn_params(component::AbstractComponent, nn_names, rng)
    # For NeuralBucket or other components with network fields
    if isdefined(Main, :NeuralBucket) && component isa Main.NeuralBucket
        flux_ps = LuxCore.initialparameters(rng, component.flux_network)
        state_ps = LuxCore.initialparameters(rng, component.state_network)
        output_ps = LuxCore.initialparameters(rng, component.output_network)

        return NamedTuple{Tuple(nn_names)}((
            ComponentVector(flux_ps),
            ComponentVector(state_ps),
            ComponentVector(output_ps)
        ))
    end

    # For composite components (HydroModel, etc.)
    if hasproperty(component, :components)
        all_nn_params = Dict{Symbol, ComponentVector}()
        for comp in component.components
            comp_nn_names = get_nn_names(comp)
            if !isempty(comp_nn_names)
                comp_params = _get_nn_params(comp, comp_nn_names, rng)
                merge!(all_nn_params, pairs(comp_params))
            end
        end
        return NamedTuple(all_nn_params)
    end

    return NamedTuple()
end

"""
    get_nn_initial_params(component::AbstractComponent; rng=Random.default_rng())

Extract initial neural network parameters from a component.

Supports:
- `NeuralFlux`: Returns parameters for the neural network
- `NeuralBucket`: Returns parameters for flux, state, and output networks
- `HydroBucket`: Returns parameters from embedded neural flux components
- `HydroModel`: Recursively collects parameters from all sub-components

Returns a NamedTuple with neural network names as keys and ComponentVector parameters as values.
Returns `nothing` if no neural networks are found.

# Examples
```julia
# For a model with NeuralFlux
model = HydroModel([neural_flux, bucket, ...])
nn_params = get_nn_initial_params(model)
# Returns: (nn_name = ComponentVector(...), ...)

# For a NeuralBucket
bucket = create_neural_bucket(...)
nn_params = get_nn_initial_params(bucket)
# Returns: (flux_net = ComponentVector(...), state_net = ComponentVector(...), output_net = ComponentVector(...))
```
"""
function get_nn_initial_params(component::AbstractComponent; rng=Random.default_rng())
    nn_names = get_nn_names(component)
    isempty(nn_names) && return nothing

    nns_nt = _get_nn_params(component, nn_names, rng)
    return isempty(nns_nt) ? nothing : nns_nt
end

# Export main functions
export sort_fluxes, sort_components, expand_component_params
export get_default_states, extract_variables
export _as_componentvector
export get_initial_params, get_nn_initial_params
export @hydroflux_for
