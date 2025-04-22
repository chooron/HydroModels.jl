"""
    HydroModel{N} <: AbstractModel

Represents a complete hydrological model that integrates multiple components for simulating water processes.

# Fields
- `components::Vector{<:AbstractComponent}`: Hydrological computation components (e.g., `HydroBucket`, `UnitHydrograph`).
- `infos::NamedTuple`: Metadata including variable names (inputs, outputs, states, params, nns) aggregated from components.
- `_varindices::AbstractVector{AbstractVector{Int}}`: Input variable indices for each component, determined automatically during construction to manage data flow.
- `_outputindices::AbstractVector{Int}`: Indices used to select and order the final model outputs from all internal variables.

# Constructor
```julia
HydroModel(; components::Vector{<:AbstractComponent}, name::Union{Symbol,Nothing}=nothing)
```
- `components`: Required. A vector containing the different model components.
- `name`: Optional symbol to identify the model. If `nothing`, a unique name is generated.

# Usage Example
```julia
# Assuming bucket1 and routing1 are predefined AbstractComponent instances
model = HydroModel(
    name = :simple_catchment,
    components = [bucket1, routing1]
)

# Simulate the model
# output = model(input_data, parameters)
```

# Notes
- The type parameter `N` encodes the model `name` at the type level for potential dispatch optimizations.
- Component connections and variable routing (`_varindices`, `_outputindices`) are handled automatically.
- The constructed `model` object is callable to perform simulations (see the function call operator documentation for details).
"""
struct HydroModel{N} <: AbstractModel
    "hydrological computation elements"
    components::Vector{<:AbstractComponent}
    "meta data of hydrological model"
    infos::NamedTuple
    "input variables index for each components"
    _varindices::AbstractVector{AbstractVector{Int}}
    "output variables index for sort output variables"
    _outputindices::AbstractVector{Int}

    function HydroModel(;
        components::Vector{C},
        name::Union{Symbol,Nothing}=nothing,
    ) where {C<:AbstractComponent}
        inputs, outputs, states = get_vars(components)
        params = reduce(union, get_params.(components))
        nns = reduce(union, get_nns.(components))
        input_idx, output_idx = _prepare_indices(components, inputs, vcat(states, outputs))
        infos = (; inputs=inputs, outputs=outputs, states=states, params=params, nns=nns)
        model_name = isnothing(name) ? Symbol("##model#", hash(infos)) : name
        new{model_name}(components, infos, input_idx, output_idx)
    end
end

"""
    _prepare_indices(components::Vector{<:AbstractComponent}, 
                    input_names::Vector{Symbol}, 
                    vcat_names::Vector{Symbol}) -> Tuple{Vector{Vector{Int}}, Vector{Int}}

Prepare input and output variable indices for connecting components in a hydrological model.

# Arguments
- `components::Vector{<:AbstractComponent}`: Vector of model components to process
- `input_names::Vector{Symbol}`: Initial list of input variable names
- `vcat_names::Vector{Symbol}`: Concatenated list of all variable names (inputs, states, outputs)

# Returns
- `Tuple{Vector{Vector{Int}}, Vector{Int}}`: A tuple containing:
  - First element: Vector of input indices for each component
  - Second element: Vector of output indices for all components
"""
function _prepare_indices(components::Vector{<:AbstractComponent}, inputs::Vector{Num}, vars::Vector{Num})
    input_idx, output_idx = Vector{Int}[], Vector{Int}()
    input_names, var_names = tosymbol.(inputs), tosymbol.(vars)
    for component in components
        #* extract input index
        tmp_input_idx = map((nm) -> findfirst(varnm -> varnm == nm, input_names), get_input_names(component))
        push!(input_idx, tmp_input_idx)
        #* extract output index
        tmp_cpt_vcat_names = vcat(get_state_names(component), get_output_names(component))
        input_names = vcat(input_names, tmp_cpt_vcat_names)
    end
    for name in var_names
        push!(output_idx, findfirst(varnm -> varnm == name, input_names))
    end
    return input_idx, output_idx
end

"""
    @hydromodel name begin
        component1
        component2
        ...
    end

Creates a HydroModel with the specified name and components.

# Arguments
- `name`: Symbol for the model name
- Components can be:
  - HydroBucket instances
  - Flux definitions (using @hydroflux, @neuralflux, etc.)
  - Other model components

# Example
```julia
@hydromodel :model1 begin
    bucket1
    @neuralflux :flux1 [y] ~ chain([x1, x2])
    bucket2
end
```
"""
macro hydromodel(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after model name"
    components = filter(x -> !(x isa LineNumberNode), expr.args)
    return esc(:(HydroModel(; name=$(name), components=[$(components...)])))
end

"""
    (model::HydroModel)(input::AbstractArray, params::ComponentVector; kwargs...)

Runs the hydrological model simulation using the provided input data and parameters.

# Arguments
- `input::AbstractArray{T,D}`: Input data. Shape can be `(variables, timesteps)` (D=2, single-node) or `(variables, nodes, timesteps)` (D=3, multi-node).
- `params::ComponentVector`: Model parameters, structured according to the components' requirements.

# Keyword Arguments
- `initstates::ComponentVector`: Optional initial states for stateful components. Defaults to component-defined defaults.
- `config::Union{NamedTuple, Vector{<:NamedTuple}}`: Optional configuration for components (e.g., ODE solver, interpolation). Can be a single `NamedTuple` applied to all, or a `Vector` for component-specific settings.
- `kwargs...`: Other keyword arguments passed down to individual components if applicable.

# Returns
- `AbstractArray`: Model output, with shape matching the input's time (and node, if D=3) dimensions: `(output_variables, timesteps)` or `(output_variables, nodes, timesteps)`.

# Example
```julia
# Assuming 'model' is a constructed HydroModel instance
# input_data might be Matrix or 3D Array
# parameters is a ComponentVector
output = model(input_data, parameters; initstates=initial_conditions)
```

# Notes
- The function orchestrates data flow through the model's components sequentially.
- Input variable dimension (`size(input, 1)`) must match `get_input_names(model)`.
- State evolution and output collection are handled internally.
"""
function (model::HydroModel{N})(
    input::AbstractArray{T,D},
    params::ComponentVector;
    kwargs...,
) where {N,T<:Number,D}
    config = get(kwargs, :config, NamedTuple())
    initstates = get(kwargs, :initstates, get_default_states(model, input))
    comp_configs = config isa NamedTuple ? fill(config, length(model.components)) : config
    @assert D in (2, 3) "input array dimension must be 2 or 3"
    @assert length(comp_configs) == length(model.components) "component configs length must be equal to components length"
    @assert size(input, 1) == length(get_input_names(model)) "input variables length must be equal to input variables length"
    outputs = input
    for (idx_, comp_, config_) in zip(model._varindices, model.components, comp_configs)
        tmp_outputs = comp_(outputs[idx_, ntuple(_ -> Colon(), D-1)...], params; initstates=initstates, config_...)
        outputs = cat(outputs, tmp_outputs, dims=1)
    end
    return outputs[model._outputindices, ntuple(_ -> Colon(), D-1)...]
end