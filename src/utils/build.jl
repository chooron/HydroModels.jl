"""
    simplify_expr(expr)

Simplifies an expression by converting it to a string and then parsing it back to an expression.
eg: convert `:((+)(x, (tanh)(y)))` to `:(x + tanh(y))`
"""
@inline simplify_expr(expr) = Meta.parse(string(expr))

"""
    extract_names(infos)

Extracts various names from a NamedTuple, standardizing the naming extraction process.
"""
function extract_names(infos)
    input_names = length(infos.inputs) == 0 ? [] : tosymbol.(infos.inputs)
    output_names = length(infos.outputs) == 0 ? [] : tosymbol.(infos.outputs)
    state_names = length(infos.states) == 0 ? [] : tosymbol.(infos.states)
    param_names = length(infos.params) == 0 ? [] : tosymbol.(infos.params)

    return (input_names=input_names, output_names=output_names,
        state_names=state_names, param_names=param_names)
end

"""
    generate_input_assignments(names; dims=0, target="inputs", prefix="")

Generates variable assignment statements, reducing repetitive assignment expressions in different functions.

# Arguments
- `names`: List of variable names
- `dims`: Specifies dimension access, 0 for one dimension, 1 for two dimensions [:, :], 2 for three dimensions [:, :, :]
- `target`: Target variable name to assign from
- `prefix`: Variable name prefix, used to distinguish variables with the same name in different contexts

# Returns
- Array of assignment statements
"""
function generate_input_assignments(names; dims=0, prefix="")
    if dims == 0
        return [:($(Symbol(prefix, i)) = inputs[$idx]) for (idx, i) in enumerate(names)]
    elseif dims == 1
        return [:($(Symbol(prefix, i)) = inputs[$idx, :]) for (idx, i) in enumerate(names)]
    elseif dims == 2
        return [:($(Symbol(prefix, i)) = inputs[$idx, :, :]) for (idx, i) in enumerate(names)]
    end
end

function generate_state_assignments(names; dims=0, prefix="", reshape=false)
    if dims == 0
        return [:($(Symbol(prefix, i)) = states[$idx]) for (idx, i) in enumerate(names)]
    elseif dims == 1
        if reshape
            reshape_expr = [:(new_states = reshape(states, length($names), :))]
            return vcat(reshape_expr, [:($(Symbol(prefix, i)) = new_states[$idx, :]) for (idx, i) in enumerate(names)])
        else
            return [:($(Symbol(prefix, i)) = states[$idx, :]) for (idx, i) in enumerate(names)]
        end
    elseif dims == 2
        return [:($(Symbol(prefix, i)) = states[$idx, :, :]) for (idx, i) in enumerate(names)]
    end
end

"""
    generate_param_assignments(param_names, target=:pas)

Generates parameter assignment statements.
"""
@inline generate_param_assignments(param_names, target=:pas) = [:($p = $(target).params.$p) for p in param_names]

"""
    generate_nn_assignments(fluxes, target=:pas)

Generates neural network assignment statements.
"""
@inline generate_nn_assignments(nn_fluxes, target=:pas) = [:($nn = $(target).nns.$nn) for nn in [get_nn_names(f)[1] for f in nn_fluxes]]

"""
    generate_compute_calls(fluxes, mode)

Generates computation call statements based on mode, dispatching to specialized functions.

# Arguments
- `fluxes`: Array of flux components
- `mode`: Computation mode, can be :state, :flux, :multi_state, :multi_flux
"""
@inline generate_compute_calls(fluxes, mode::Symbol) = generate_compute_calls(fluxes, Val(mode))

"""
    generate_compute_calls(fluxes, ::Val{:state})

Generates computation calls for state mode (scalar inputs and outputs).
"""
function generate_compute_calls(fluxes, ::Val{:state})
    compute_calls = []

    for f in fluxes
        input_names, output_names, _ = get_var_names(f)
        if f isa AbstractNeuralFlux
            nn_names = get_nn_names(f)[1]
            append!(compute_calls, [:($(f.infos[:nn_inputs]) = [$(input_names...)])])
            push!(compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(nn_names))))
            append!(compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i]) for (i, nm) in enumerate(output_names)])
        else
            # Process regular flux
            append!(compute_calls, [:($nm = $(toexpr(expr))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    return compute_calls
end

"""
    generate_compute_calls(fluxes, ::Val{:flux})

Generates computation calls for flux mode (vectorized inputs and outputs).
"""
function generate_compute_calls(fluxes, ::Val{:flux})
    compute_calls = []

    for f in fluxes
        input_names, output_names, _ = get_var_names(f)
        if f isa AbstractNeuralFlux
            nn_names = get_nn_names(f)[1]
            push!(compute_calls, :($(f.infos[:nn_inputs]) = stack([$(input_names...)], dims=1)))
            push!(compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(nn_names))))
            append!(compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :]) for (i, nm) in enumerate(output_names)])
        else
            # Process regular flux
            append!(compute_calls, [:($nm = @. $(simplify_expr(toexpr(expr)))) for (nm, expr) in zip(output_names, f.exprs)])
        end
    end

    return compute_calls
end

"""
    generate_compute_calls(fluxes, ::Val{:multi_state})

Generates computation calls for multi-state mode (batch processing for state calculations).
"""
function generate_compute_calls(fluxes, ::Val{:multi_state})
    compute_calls = []

    for f in fluxes
        input_names, output_names, _ = get_var_names(f)
        if f isa AbstractNeuralFlux
            nn_names = get_nn_names(f)[1]
            push!(compute_calls, :($(f.infos[:nn_inputs]) = stack([$(input_names...)], dims=1)))
            push!(compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(nn_names))))
            append!(compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :]) for (i, nm) in enumerate(output_names)])
        else
            # Process regular flux
            append!(compute_calls, [:($nm = @. $(simplify_expr(toexpr(expr)))) for (nm, expr) in zip(output_names, f.exprs)])
        end
    end

    return compute_calls
end

"""
    generate_compute_calls(fluxes, ::Val{:multi_flux})

Generates computation calls for multi-flux mode (multi-dimensional batch processing).
"""
function generate_compute_calls(fluxes, ::Val{:multi_flux})
    compute_calls = []

    for f in fluxes
        input_names, output_names, _ = get_var_names(f)
        if f isa AbstractNeuralFlux
            nn_names = get_nn_names(f)[1]
            push!(compute_calls, :($(f.infos[:nn_inputs]) = stack([$(input_names...)], dims=1)))
            push!(compute_calls, :($(f.infos[:nn_outputs]) = stack($(f.func).(eachslice($(f.infos[:nn_inputs]), dims=2), Ref($(nn_names))), dims=2)))
            append!(compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :, :]) for (i, nm) in enumerate(output_names)])
        else
            append!(compute_calls, [:($nm = @. $(simplify_expr(toexpr(expr)))) for (nm, expr) in zip(output_names, f.exprs)])
        end
    end

    return compute_calls
end

"""
    generate_return_expression(output_names, dfluxes=nothing; mode=:normal)

Generates return statements based on mode, dispatching to specialized functions.

# Arguments
- `output_names`: Output variable names
- `dfluxes`: State differential component array
- `mode`: Return mode, can be :normal, :state, :multi_state, :route

# Returns
- Return expression
"""
@inline function generate_return_expression(; output_names=nothing, dfluxes=nothing, mode=:normal)
    generate_return_expression(output_names, dfluxes, Val(mode))
end

"""
    generate_return_expression(output_names, dfluxes, ::Val{:normal})

Generates standard return expression with output names.
"""
function generate_return_expression(output_names, ::Any, ::Val{:normal})
    return :(return [$(output_names...)])
end

"""
    generate_return_expression(output_names, dfluxes, ::Val{:state})

Generates return expression for state differentials.
"""
function generate_return_expression(::Any, dfluxes, ::Val{:state})
    return :(return [$(map(expr -> :($(toexpr(expr))), reduce(vcat, get_exprs.(dfluxes)))...)])
end

"""
    generate_return_expression(output_names, dfluxes, ::Val{:multi_state})

Generates return expression for multi-state differentials.
"""
function generate_return_expression(::Any, dfluxes, ::Val{:multi_state})
    return :(return reduce(vcat, [$(map(expr -> :(@. $(simplify_expr(toexpr(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)]))
end

"""
    generate_return_expression(output_names, dfluxes, ::Val{:route})

Generates return expression for routing calculations.
"""
function generate_return_expression(output_names, dfluxes, ::Val{:output_state})
    return :(return $(output_names...), $(map(expr -> :(@. $(simplify_expr(toexpr(expr)))), reduce(vcat, get_exprs.(dfluxes)))...))
end

"""
    generate_return_expression(output_names, dfluxes, ::Val{:nnlayer})

Generates return expression for neural network layer.
"""
function generate_return_expression(output_names, dfluxes, ::Val{:nnlayer})
    dfluxes_exprs = reduce(vcat, get_exprs.(dfluxes))
    return :(return [$(output_names...)], [$(map(expr -> :(@. $(simplify_expr(toexpr(expr)))), dfluxes_exprs)...)])
end

"""
    build_flux_func(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num})

Generates a runtime function for flux calculations based on symbolic expressions.
"""
function build_flux_func(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num})
    names = (
        input_names=length(inputs) == 0 ? [] : tosymbol.(inputs),
        output_names=length(outputs) == 0 ? [] : tosymbol.(outputs),
        param_names=length(params) == 0 ? [] : tosymbol.(params)
    )

    flux_exprs = map(expr -> :(@. $(simplify_expr(toexpr(expr)))), exprs)
    input_assign_calls = generate_input_assignments(names.input_names)
    params_assign_calls = generate_param_assignments(names.param_names)
    compute_calls = [:($o = $expr) for (o, expr) in zip(names.output_names, flux_exprs)]
    return_expr = :(return [$(names.output_names...)])

    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, pas)
        $(meta_exprs...)
        $(input_assign_calls...)
        $(params_assign_calls...)
        $(compute_calls...)
        $(return_expr)
    end)
    return @RuntimeGeneratedFunction(flux_func_expr)
end

"""
    build_single_bucket_func(fluxes, dfluxes, infos)

Builds runtime functions for flux calculations and state differentials in a hydrological model element.
"""
function build_single_bucket_func(fluxes::Vector{<:AbstractFlux}, dfluxes::Vector{<:AbstractStateFlux}, infos::NamedTuple)
    names = extract_names(infos)

    # Define variable assignments
    input_define_calls_1 = generate_input_assignments(names.input_names, dims=0)
    state_define_calls_1 = generate_state_assignments(names.state_names, dims=0)
    input_define_calls_2 = generate_input_assignments(names.input_names, dims=1)
    state_define_calls_2 = generate_state_assignments(names.state_names, dims=1)
    params_assign_calls = generate_param_assignments(names.param_names)
    nn_params_assign_calls = generate_nn_assignments(filter(f -> f isa AbstractNeuralFlux, fluxes))
    define_calls_1 = [input_define_calls_1..., state_define_calls_1..., params_assign_calls..., nn_params_assign_calls...]
    define_calls_2 = [input_define_calls_2..., state_define_calls_2..., params_assign_calls..., nn_params_assign_calls...]

    # Generate computation calls
    flux_compute_calls = generate_compute_calls(fluxes, :flux)

    # Return expression
    return_flux = generate_return_expression(output_names=names.output_names, mode=:normal)

    # Create function expression
    flux_func_expr = :(function (inputs, states, pas)
        $([:(Base.@_inline_meta)]...)
        $(define_calls_2...)
        $(flux_compute_calls...)
        $(return_flux)
    end)
    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)

    # If state variables exist, create differential function
    if !isempty(names.state_names)
        state_compute_calls = generate_compute_calls(fluxes, :state)
        return_state = generate_return_expression(dfluxes=dfluxes, mode=:state)
        diff_func_expr = :(function (inputs, states, pas)
            $([:(Base.@_inline_meta)]...)
            $(define_calls_1...)
            $(state_compute_calls...)
            $(return_state)
        end)
        generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
        return generated_flux_func, generated_diff_func
    else
        return generated_flux_func, nothing
    end
end

"""
    build_multi_bucket_func(fluxes, dfluxes, infos)

Builds runtime functions for multi-dimensional batch processing in a hydrological model.
"""
function build_multi_bucket_func(fluxes::Vector{<:AbstractFlux}, dfluxes::Vector{<:AbstractStateFlux}, infos::NamedTuple)
    names = extract_names(infos)

    # Define variable assignments for different dimensions
    input_define_calls_1 = generate_input_assignments(names.input_names, dims=1)
    state_define_calls_1 = generate_state_assignments(names.state_names, dims=1, reshape=true)

    input_define_calls_2 = generate_input_assignments(names.input_names, dims=2)
    state_define_calls_2 = generate_state_assignments(names.state_names, dims=2)

    params_assign_calls = generate_param_assignments(names.param_names)
    nn_params_assign_calls = generate_nn_assignments(filter(f -> f isa AbstractNeuralFlux, fluxes))

    define_calls_1 = [input_define_calls_1..., state_define_calls_1..., params_assign_calls..., nn_params_assign_calls...]
    define_calls_2 = [input_define_calls_2..., state_define_calls_2..., params_assign_calls..., nn_params_assign_calls...]

    # Generate computation calls
    multi_flux_compute_calls = generate_compute_calls(fluxes, :multi_flux)

    # Return expression
    return_flux = generate_return_expression(output_names=names.output_names, mode=:normal)

    # Create multi-dimensional flux function
    multi_flux_func_expr = :(function (inputs, states, pas)
        $([:(Base.@_inline_meta)]...)
        $(define_calls_2...)
        $(multi_flux_compute_calls...)
        $(return_flux)
    end)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)

    # If state variables exist, create multi-dimensional differential function
    if !isempty(names.state_names)
        multi_state_compute_calls = generate_compute_calls(fluxes, :multi_state)
        return_multi_state = generate_return_expression(dfluxes=dfluxes, mode=:multi_state)
        multi_diff_func_expr = :(function (inputs, states, pas)
            $([:(Base.@_inline_meta)]...)
            $(define_calls_1...)
            $(multi_state_compute_calls...)
            $(return_multi_state)
        end)
        generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
        return generated_multi_flux_func, generated_multi_diff_func
    else
        return generated_multi_flux_func, nothing
    end
end

"""
    build_route_func(fluxes, dfluxes, infos)

Builds routing calculation functions for a hydrological model.
"""
function build_route_func(fluxes::AbstractVector{<:AbstractFlux}, dfluxes::AbstractVector{<:AbstractStateFlux}, infos::NamedTuple)
    names = extract_names(infos)

    # Define variable assignments for different dimensions
    input_define_calls_1 = generate_input_assignments(names.input_names, dims=1)
    state_define_calls_1 = generate_state_assignments(names.state_names, dims=1, reshape=true)

    input_define_calls_2 = generate_input_assignments(names.input_names, dims=2)
    state_define_calls_2 = generate_state_assignments(names.state_names, dims=2)

    params_assign_calls = generate_param_assignments(names.param_names)
    nn_params_assign_calls = generate_nn_assignments(filter(f -> f isa AbstractNeuralFlux, fluxes))

    define_calls_1 = [input_define_calls_1..., state_define_calls_1..., params_assign_calls..., nn_params_assign_calls...]
    define_calls_2 = [input_define_calls_2..., state_define_calls_2..., params_assign_calls..., nn_params_assign_calls...]

    # Generate computation calls
    multi_flux_compute_calls = generate_compute_calls(fluxes, :multi_flux)

    # Return expression
    return_flux = generate_return_expression(output_names=names.output_names, mode=:normal)

    # Create multi-dimensional flux function
    multi_flux_func_expr = :(function (inputs, states, pas)
        $([:(Base.@_inline_meta)]...)
        $(define_calls_2...)
        $(multi_flux_compute_calls...)
        $(return_flux)
    end)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)

    # If state variables exist, create multi-dimensional differential function
    multi_state_compute_calls = generate_compute_calls(fluxes, :multi_state)
    return_multi_state = generate_return_expression(output_names=names.output_names, dfluxes=dfluxes, mode=:output_state)
    multi_diff_func_expr = :(function (inputs, states, pas)
        $([:(Base.@_inline_meta)]...)
        $(define_calls_1...)
        $(multi_state_compute_calls...)
        $(return_multi_state)
    end)
    generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
    return generated_multi_flux_func, generated_multi_diff_func
end


"""
    build_nnlayer_func(fluxes, dfluxes, infos)

Builds a neural network layer function for a hydrological model.
"""
function build_nnlayer_func(fluxes::Vector{<:AbstractHydroFlux}, dfluxes::Vector{<:AbstractStateFlux}, infos::NamedTuple)
    names = extract_names(infos)

    # Define variable assignments
    input_define_calls = generate_input_assignments(names.input_names, dims=1)
    state_define_calls = generate_state_assignments(names.state_names, dims=1)
    params_assign_calls = generate_param_assignments(names.param_names)
    nn_params_assign_calls = [:($(nflux.infos[:nns][1]) = pas.nns.$(nflux.infos[:nns][1])) for nflux in filter(f -> f isa AbstractNeuralFlux, fluxes)]
    define_calls = [input_define_calls..., state_define_calls..., params_assign_calls..., nn_params_assign_calls...]

    # Generate computation calls
    compute_calls = generate_compute_calls(fluxes, :flux)

    # Return expression
    return_flux = generate_return_expression(output_names=names.output_names, dfluxes=dfluxes, mode=:output_state)

    # Create function expression
    func_expr = :(function (inputs, pas, states)
        $([:(Base.@_inline_meta)]...)
        $(define_calls...)
        $(compute_calls...)
        $(return_flux)
    end)

    return @RuntimeGeneratedFunction(func_expr)
end

"""
    build_uh_func(uh_pairs, params, max_lag)

Builds unit hydrograph functions.
"""
function build_uh_func(uh_pairs::AbstractVector{<:Pair}, params::AbstractVector, max_lag::Number)
    conditions_rev = vcat([0], reverse(first.(uh_pairs)))
    values_rev = reverse(last.(uh_pairs))
    param_names = tosymbol.(params)
    params_assign_calls = generate_param_assignments(param_names)

    values_exprs = map(eachindex(values_rev)) do i
        :(
            if $(toexpr(conditions_rev[i])) < t < $(toexpr(conditions_rev[i+1]))
                return $(toexpr(values_rev[i]))
            end
        )
    end
    default_return = :(return 0.0)

    uh_func_expr = :(function (t, pas)
        $(params_assign_calls...)
        $(values_exprs...)
        $(default_return)
    end)

    max_lag_expr = :(function (pas)
        $(params_assign_calls...)
        ceil($(toexpr(max_lag)))
    end)

    return @RuntimeGeneratedFunction(uh_func_expr), @RuntimeGeneratedFunction(max_lag_expr)
end
