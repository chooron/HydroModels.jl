"""
$(SIGNATURES)

Simplify a Julia expression. For example, convert `:((+)(x, (tanh)(y)))` to `:(x + tanh(y))`.
"""
@inline simplify_expr(expr) = Meta.parse(string(expr))

"""
$(SIGNATURES)

Generate variable assignment expressions from a target array.
"""
@inline generate_var_assignments(; vars, target, dims=0, prefix="") = begin
    [:($(Symbol(prefix, i)) = $(target)[$idx, ntuple(_ -> Colon(), $dims)...]) for (idx, i) in enumerate(vars)]
end

"""
$(SIGNATURES)

Generate parameter assignment expressions.
"""
@inline generate_param_assignments(; params, target=:pas) = [:($p = $(target).params.$p) for p in params]

"""
$(SIGNATURES)

Generate neural network assignment expressions.
"""
@inline generate_nn_assignments(; nnfluxes, target=:pas) = [:($nn = $(target).nns.$nn) for nn in [get_nn_names(f)[1] for f in nnfluxes]]


"""
    generate_flux_expression(flux::Union{AbstractHydroFlux, AbstractNeuralFlux}, ::Val{dims})

Generate flux computation expressions, dispatching on flux type and dimensionality (`dims`).
"""
@inline generate_flux_expression(flux::AbstractHydroFlux, ::Val{0}) = [:($nm = $(toexpr(expr))) for (nm, expr) in zip(get_output_names(flux), flux.exprs)]
@inline generate_flux_expression(flux::AbstractHydroFlux, ::Val{1}) = [:($nm = @. $(simplify_expr(toexpr(expr)))) for (nm, expr) in zip(get_output_names(flux), flux.exprs)]
@inline generate_flux_expression(flux::AbstractHydroFlux, ::Val{2}) = generate_flux_expression(flux, Val(1))
@inline generate_flux_expression(flux::AbstractNeuralFlux, ::Val{0}) = begin
    nn_names = get_nn_names(flux)[1]
    nn_inputs, nn_outputs = Symbol(nn_names, "_input"), Symbol(nn_names, "_output")
    [
        :($(nn_inputs) = [$(flux.infos.inputs...)] |> $(flux.norm_func)),
        :($(nn_outputs) = $(flux.chain_func)($(nn_inputs), $(nn_names))),
        [:($(nm) = $(nn_outputs)[$i]) for (i, nm) in enumerate(get_output_names(flux))]...
    ]
end
@inline generate_flux_expression(flux::AbstractNeuralFlux, ::Val{1}) = begin
    nn_names = get_nn_names(flux)[1]
    nn_inputs, nn_outputs = Symbol(nn_names, "_input"), Symbol(nn_names, "_output")
    [
        :($(nn_inputs) = stack([$(flux.infos.inputs...)], dims=1) |> $(flux.norm_func)),
        :($(nn_outputs) = $(flux.chain_func)($(nn_inputs), $(nn_names))),
        [:($(nm) = $(nn_outputs)[$i, :]) for (i, nm) in enumerate(get_output_names(flux))]...
    ]
end

@inline generate_flux_expression(flux::AbstractNeuralFlux, ::Val{2}) = begin
    nn_names = get_nn_names(flux)[1]
    nn_inputs, nn_outputs = Symbol(nn_names, "_input"), Symbol(nn_names, "_output")
    [
        :($(nn_inputs) = stack([$(flux.infos.inputs...)], dims=1) |> $(flux.norm_func)),
        :($(nn_outputs) = $(flux.chain_func)($(nn_inputs), $(nn_names))),
        [:($(nm) = $(nn_outputs)[$i, :, :]) for (i, nm) in enumerate(get_output_names(flux))]...
    ]
end

"""
$(SIGNATURES)

Generate a vector of computation expressions for a list of fluxes.
"""
@inline generate_compute_calls(; fluxes, dims=0) = vcat(map(f -> generate_flux_expression(f, Val(dims)), fluxes)...)

"""
$(SIGNATURES)

Generate the `return` expression for state differentials, with an option for broadcasting.
"""
@inline generate_states_expression(; dfluxes, broadcast=false) = _generate_states_expression(dfluxes, Val(broadcast))
@inline _generate_states_expression(dfluxes, ::Val{false}) = :(return [$(map(expr -> toexpr(expr), reduce(vcat, get_exprs.(dfluxes)))...)])
@inline _generate_states_expression(dfluxes, ::Val{true}) = :(return vcat($(map(expr -> :(@. $(simplify_expr(toexpr(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)))

"""
$(TYPEDSIGNATURES)

Topologically sort a vector of flux components based on their variable dependencies.

This function builds a directed graph representing the dependencies between flux components
(where edges connect input variables to output variables) and performs a topological sort
to determine the correct calculation order.

# Examples
```julia
fluxes = [flux1, flux2, flux3]
sorted_fluxes = sort_fluxes(fluxes)
# sorted_fluxes now contains the components in the correct calculation order
```
"""
function sort_fluxes(fluxes::AbstractVector{<:AbstractComponent})
    input_names = reduce(union, get_input_names.(fluxes))
    output_names = reduce(union, get_output_names.(fluxes))
    input_names = setdiff(input_names, output_names)
    output_names = setdiff(output_names, input_names)

    # Build a named tuple mapping output names to their corresponding flux instances
    fluxes_ntp = reduce(merge, map(fluxes) do flux
        tmp_output_names = get_output_names(flux)
        NamedTuple{Tuple(tmp_output_names)}(repeat([flux], length(tmp_output_names)))
    end)

    # Construct a directed calculation graph
    var_names = vcat(input_names, output_names)
    var_names_ntp = NamedTuple{Tuple(var_names)}(1:length(var_names))
    digraph = SimpleDiGraph(length(var_names))
    for flux in fluxes
        tmp_input_names, tmp_output_names = get_input_names(flux), get_output_names(flux)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end

    # Sort fluxes based on the topological order of the directed graph
    sorted_fluxes = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_flux = fluxes_ntp[tmp_var_nm]
            if !(tmp_flux in sorted_fluxes)
                push!(sorted_fluxes, tmp_flux)
            end
        end
    end
    sorted_fluxes
end

"""
$(TYPEDSIGNATURES)

Topologically sort a vector of hydrological components based on their variable dependencies.

This function builds a directed graph of dependencies between components (inputs to outputs/states)
and uses a topological sort to find the correct execution order.

# Examples
```julia
components = [bucket1, flux1, route1]
sorted_components = sort_components(components)
# sorted_components now contains the components in the correct calculation order
```
"""
function sort_components(components::AbstractVector{<:AbstractComponent})
    input_names, output_names, state_names = get_var_names(components)
    components_ntp = reduce(merge, map(components) do component
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        tmp_output_state_names = vcat(tmp_output_names, tmp_state_names)
        NamedTuple{Tuple(tmp_output_state_names)}(repeat([component], length(tmp_output_state_names)))
    end)
    var_names = reduce(union, [input_names, output_names, state_names])
    var_names_ntp = NamedTuple{Tuple(var_names)}(collect(1:length(var_names)))
    digraph = SimpleDiGraph(length(var_names))
    for component in components
        tmp_input_names, tmp_output_names, tmp_state_names = get_var_names(component)
        tmp_output_names = vcat(tmp_output_names, tmp_state_names)
        for ipnm in tmp_input_names
            for opnm in tmp_output_names
                add_edge!(digraph, var_names_ntp[ipnm], var_names_ntp[opnm])
            end
        end
    end
    sorted_components = AbstractComponent[]
    for idx in topological_sort(digraph)
        tmp_var_nm = var_names[idx]
        if (tmp_var_nm in output_names)
            tmp_component = components_ntp[tmp_var_nm]
            if !(tmp_component in sorted_components)
                push!(sorted_components, tmp_component)
            end
        end
    end
    sorted_components
end


"""
$(TYPEDSIGNATURES)

Expand specified parameters of a `ComponentVector` using a given index.

This is useful for spatial modeling, where a parameter vector needs to be expanded to match the grid size.
"""
function expand_component_params(pas::ComponentVector, pnames::AbstractVector, ptyidx::AbstractVector)
    params_ntp = NamedTuple(pas[:params])
    expand_params = (; params=NamedTuple{Tuple(pnames)}([params_ntp[p][ptyidx] for p in pnames]))
    return merge(pas |> NamedTuple, expand_params) |> ComponentVector
end

"""
$(TYPEDSIGNATURES)

Get a `ComponentVector` of default (zeroed) states for a given component.

The shape of the states depends on the input array dimensions:
- 2D input `(time, features)`: Returns a vector of zeros for each state.
- 3D input `(time, locations, features)`: Returns a `NamedTuple` of state vectors for each location.
"""
function get_default_states(component::AbstractComponent, input::AbstractArray{T,2}) where {T}
    state_names = get_state_names(component)
    states_axes = Axis(NamedTuple{Tuple(state_names)}(1:length(state_names)))
    return ComponentVector(zeros(T, length(state_names)), states_axes)
end

function get_default_states(component::AbstractComponent, input::AbstractArray{T,3}) where {T}
    state_names = get_state_names(component)
    return ComponentVector(NamedTuple{Tuple(state_names)}(fill(zeros(T, size(input, 2)), length(state_names))))
end


"""
$(SIGNATURES)

Extract all variable symbols from a Julia expression.
"""
function extract_variables(expr)
    function _extract_vars!(expr::Any, vars::Set{Symbol})
    end

    function _extract_vars!(s::Symbol, vars::Set{Symbol})
        push!(vars, s)
    end

    function _extract_vars!(ex::Expr, vars::Set{Symbol})
        if ex.head == :call
            for arg in ex.args[2:end]
                _extract_vars!(arg, vars)
            end
        elseif ex.head == :.
            _extract_vars!(ex.args[1], vars)
        elseif ex.head == :(=) || ex.head == :function || ex.head == :macro
        else
            for arg in ex.args
                _extract_vars!(arg, vars)
            end
        end
    end
    vars = Set{Symbol}()
    _extract_vars!(expr, vars)
    return vars
end

"""
$(SIGNATURES)

Check if all variables within an expression are defined in the current module. Throws an error for undefined variables.
"""
function check_defined_variables(expr)
    for var_name in extract_variables(expr)
        if !isdefined(__module__, var_name)
            expr_str = string(expr)
            return error("Undefined variable '$(var_name)' detected in expression: `$(expr_str)`")
        end
    end
    return true
end


macro hydroflux_for(args...)
    name = length(args) == 1 ? nothing : args[1]
    eqs_expr = length(args) == 1 ? args[1] : args[2]

    if !Meta.isexpr(eqs_expr, :for)
        return :(error("@hydroflux_for macro must be used with a 'for' loop expression."))
    end

    loop_var = eqs_expr.args[1].args[1]
    range_expr = eqs_expr.args[1].args[2]
    range_val = if Meta.isexpr(range_expr, :call) && range_expr.args[1] == :(:)
        if length(range_expr.args) == 3
            range_expr.args[2]:range_expr.args[3]
        elseif length(range_expr.args) == 4
            range_expr.args[2]:range_expr.args[3]:range_expr.args[4]
        end
    else
        eval(range_expr)
    end

    loop_body = eqs_expr.args[2]
    if !Meta.isexpr(loop_body, :block)
        loop_body = Expr(:block, loop_body)
    end

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

    for var_name in extract_variables(vect_eqs_expr)
        if !@isdefined(var_name)
            expr_str = string(vect_eqs_expr)
            return :(error("Undefined variable '", $(string(var_name)), "' detected in expression: `", $expr_str, "`"))
        end
    end

    return esc(:(HydroFlux(exprs=$vect_eqs_expr, name=$name)))
end