"""
    HydroBucket{S, MS, FF, OF, HT, NT} <: AbstractHydroBucket

A container that organizes and executes a collection of hydrological flux components (`fluxes` and `dfluxes`).

It automatically compiles the components into efficient, callable functions for both single-node and multi-node simulations.

$(FIELDS)
"""
struct HydroBucket{S,MS,FF,OF,HT,NT} <: AbstractHydroBucket
    "bucket name"
    name::Symbol
    "Generated function for calculating all hydrological fluxes. (Supports single-node data, multi-nodes data)"
    flux_func::FF
    "Generated function for ordinary differential equations (ODE) calculations, or nothing if no ODE calculations are needed. (Supports single-node data)"
    ode_func::OF
    "nodes type"
    hru_types::HT
    "Metadata about the bucket, including input, output, state, parameter, and neural network names."
    infos::NT

    function HydroBucket(;
        name::Optional{Symbol}=nothing,
        fluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux}=StateFlux[],
        hru_types::Vector{Int}=Int[]
    )
        inputs, outputs, states = get_var_names(vcat(fluxes, dfluxes))
        params = reduce(union, get_param_names.(vcat(fluxes, dfluxes)))
        nns = reduce(union, get_nn_names.(fluxes))
        infos = HydroModelCore.HydroInfos(
            inputs=inputs, states=states, outputs=outputs,
            params=params, nns=nns
        )
        flux_func, ode_func = build_bucket_func(fluxes, dfluxes, infos, length(hru_types) > 1)
        bucket_name = isnothing(name) ? Symbol("##bucket#", hash(infos)) : name
        hasstates, ismul = !isempty(states), length(hru_types) > 1
        return new{hasstates,ismul,typeof(flux_func),typeof(ode_func),typeof(hru_types),typeof(infos)}(
            bucket_name, flux_func, ode_func, hru_types, infos
        )
    end
end

"""
    @hydrobucket [name] begin ... end

A macro to conveniently create a `HydroBucket`.

# Usage
The macro takes an optional name and a `begin...end` block containing the component definitions.

```julia
@hydrobucket :my_bucket begin
    fluxes = begin
        # ... list of HydroFlux, NeuralFlux, etc.
        flux1
        flux2
    end
    dfluxes = begin
        # ... list of StateFlux
        state_flux1
    end
    hru_types = [1, 1, 2, 2] # Optional
end
```
"""
macro hydrobucket(args...)
    name = length(args) == 1 ? nothing : args[1]
    expr = length(args) == 1 ? args[1] : args[2]

    @assert Meta.isexpr(expr, :block) "Expected a begin...end block after bucket name"
    fluxes_expr, dfluxes_expr, hru_types_val = nothing, nothing, Int[]
    for assign in filter(x -> !(x isa LineNumberNode), expr.args)
        @assert Meta.isexpr(assign, :(=)) "Expected assignments in the form 'fluxes = begin...end'"
        lhs, rhs = assign.args
        if lhs == :fluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'fluxes' to be defined in a begin...end block"
            fluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :dfluxes
            @assert Meta.isexpr(rhs, :block) "Expected 'dfluxes' to be defined in a begin...end block"
            dfluxes_expr = Expr(:vect, filter(x -> !(x isa LineNumberNode), rhs.args)...)
        elseif lhs == :hru_types
            hru_types_val = rhs
        else
            error("Unknown assignment: $lhs. Expected 'fluxes', 'dfluxes', or 'hru_types'")
        end
    end

    kwargs = [:($:(name=$name)), :($:(fluxes=$fluxes_expr)), :($:(hru_types=$hru_types_val))]
    if !isnothing(dfluxes_expr)
        push!(kwargs, :($:(dfluxes=$dfluxes_expr)))
    end
    return esc(:(HydroBucket(; $(kwargs...))))
end

"""
$(TYPEDSIGNATURES)

Re-creates a `HydroBucket` with a new set of `hru_types`.

Note: This is not an in-place operation and returns a new `HydroBucket` instance.
"""
function set_hru_types!(bucket::HydroBucket, hru_types::Vector{Int})::HydroBucket
    bucket = HydroBucket(bucket.name, bucket.flux_func, bucket.ode_func, hru_types, bucket.infos)
end

build_bucket_func(
    fluxes::Vector{<:AbstractHydroFlux}, dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos, multiply::Bool
) = begin
    _build_bucket_func(fluxes, dfluxes, infos, Val(multiply))
end


function _build_bucket_func(
    fluxes::Vector{<:AbstractHydroFlux}, dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos, ::Val{false}
)
    nn_fluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)
    define_calls_1 = [
        generate_var_assignments(vars=infos.inputs, target=:inputs, dims=0)...,
        generate_var_assignments(vars=infos.states, target=:states, dims=0)...,
        generate_param_assignments(params=infos.params)...,
        generate_nn_assignments(nnfluxes=nn_fluxes)...
    ]
    define_calls_2 = [
        generate_var_assignments(vars=infos.inputs, target=:inputs, dims=1)...,
        generate_var_assignments(vars=infos.states, target=:states, dims=1)...,
        generate_param_assignments(params=infos.params)...,
        generate_nn_assignments(nnfluxes=nn_fluxes)...
    ]

    flux_func_expr = :(function (inputs, states, pas)
        Base.@_inline_meta
        $(define_calls_2...)
        $(generate_compute_calls(fluxes=fluxes, dims=1)...)
        return [$((infos.outputs)...)]
    end)
    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)

    if !isempty(infos.states)
        diff_func_expr = :(function (inputs, states, pas)
            Base.@_inline_meta
            $(define_calls_1...)
            $(generate_compute_calls(fluxes=fluxes, dims=0)...)
            $(generate_states_expression(dfluxes=dfluxes, broadcast=false))
        end)
        generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
        return generated_flux_func, generated_diff_func
    else
        return generated_flux_func, (_) -> nothing
    end
end

function _build_bucket_func(
    fluxes::Vector{<:AbstractHydroFlux}, dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos, ::Val{true}
)
    nn_fluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)

    define_calls_1 = [
        generate_var_assignments(vars=infos.inputs, target=:inputs, dims=1)...,
        generate_var_assignments(vars=infos.states, target=:states, dims=1)...,
        generate_param_assignments(params=infos.params)...,
        generate_nn_assignments(nnfluxes=nn_fluxes)...
    ]
    define_calls_2 = [
        generate_var_assignments(vars=infos.inputs, target=:inputs, dims=2)...,
        generate_var_assignments(vars=infos.states, target=:states, dims=2)...,
        generate_param_assignments(params=infos.params)...,
        generate_nn_assignments(nnfluxes=nn_fluxes)...
    ]

    multi_flux_func_expr = :(function (inputs, states, pas)
        Base.@_inline_meta
        $(define_calls_2...)
        $(generate_compute_calls(fluxes=fluxes, dims=1)...)
        return [$((infos.outputs)...)]
    end)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)

    if !isempty(infos.states)
        multi_diff_func_expr = :(function (inputs, states, pas)
            Base.@_inline_meta
            $(define_calls_1...)
            $(generate_compute_calls(fluxes=fluxes, dims=1)...)
            $(generate_states_expression(dfluxes=dfluxes, broadcast=true))
        end)
        generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
        return generated_multi_flux_func, generated_multi_diff_func
    else
        return generated_multi_flux_func, (_) -> nothing
    end
end


"""
    (bucket::HydroBucket)(input, params; kwargs...)

Executes the `HydroBucket` model. This is the functor implementation for `HydroBucket`.

The method dispatches based on the input data dimensions (2D for single-node, 3D for multi-node) and whether the bucket contains state variables. It solves the ODEs for state variables (if any) and then computes the output fluxes.

Common `kwargs` include `solver`, `interp`, `timeidx`, and `initstates`.
"""
function (bucket::HydroBucket{true,false})(
    input::AbstractArray{T,2}, params::ComponentVector;
    solver::S=ManualSolver(mutable=true),
    interp::I=DirectInterpolation,
    timeidx::AbstractVector=collect(1:size(input, 2)),
    initstates::AbstractArray=zeros(eltype(params), length(get_state_names(bucket))),
    kwargs...
) where {T,S,I}
    param_vec, params_axes = Vector(params), getaxes(params)
    itpfuncs = interp(input, timeidx)
    solved_states = solver(
        (u, p, t) -> bucket.ode_func(itpfuncs(t), u, ComponentVector(p, params_axes)),
        param_vec, Vector(initstates), timeidx
    )
    flux_output = bucket.flux_func(input, solved_states, params)
    vcat(solved_states, stack(flux_output, dims=1))
end

function (bucket::HydroBucket{true,true})(
    input::AbstractArray{T,3}, params::ComponentVector;
    solver::S=ManualSolver(mutable=true),
    interp::I=DirectInterpolation,
    timeidx::AbstractVector=collect(1:size(input, 3)),
    initstates::AbstractArray=zeros(eltype(params), length(get_state_names(bucket)) * size(input, 2)),
    device=solver.dev,
    kwargs...
) where {T,S,I}
    input_dims, num_nodes, time_len = size(input)
    new_params = expand_component_params(params, get_param_names(bucket), bucket.hru_types)
    params_vec, params_axes = Vector(new_params) |> device, getaxes(new_params)

    num_states = length(get_state_names(bucket))
    itpfuncs = interp(reshape(input, input_dims * num_nodes, :), timeidx)
    solved_states = solver(
        (u, p, t) -> bucket.ode_func(
            reshape(itpfuncs(t), input_dims, num_nodes),
            reshape(u, num_nodes, num_states)', ComponentVector(p, params_axes),
        ),
        params_vec, initstates, timeidx
    )
    solved_states_reshape = permutedims(reshape(solved_states, num_nodes, num_states, time_len), (2, 1, 3))
    output = bucket.flux_func(input, solved_states_reshape, new_params)
    cat(solved_states_reshape, stack(output, dims=1), dims=1)
end

(bucket::HydroBucket{false,MS})(input::AbstractArray{T,2}, params::ComponentVector; kwargs...) where {T,MS} = begin
    stack(bucket.flux_func(input, nothing, params), dims=1)
end

(bucket::HydroBucket{false,MS})(input::AbstractArray{T,3}, params::ComponentVector; kwargs...) where {T,MS} = begin
    new_params = expand_component_params(params, get_param_names(bucket), bucket.hru_types)
    stack(bucket.flux_func(input, nothing, new_params), dims=1)
end