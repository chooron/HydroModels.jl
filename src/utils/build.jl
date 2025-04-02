"""
    build_flux_func(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num})

Generates a runtime function that computes flux calculations based on symbolic expressions.

# Arguments
- `inputs::Vector{Num}`: Vector of symbolic input variables that will be provided at runtime
- `outputs::Vector{Num}`: Vector of symbolic output variables that will be computed
- `params::Vector{Num}`: Vector of symbolic parameters used in the flux calculations
- `exprs::Vector{Num}`: Vector of symbolic expressions defining how outputs are computed from inputs and parameters

# Returns
- A runtime-generated function with signature `(inputs, pas)` where:
  - `inputs`: Vector of input values corresponding to the symbolic inputs
  - `pas`: A parameter struct containing fields matching the parameter names

# Details
The function generates code that:
1. Assigns input values from the input vector to local variables
2. Retrieves parameter values from the parameter struct
3. Computes outputs using the provided expressions
4. Returns a vector of computed outputs
"""
function build_flux_func(inputs::Vector{Num}, outputs::Vector{Num}, params::Vector{Num}, exprs::Vector{Num})
    input_names, output_names = Symbolics.tosymbol.(inputs), Symbolics.tosymbol.(outputs)
    param_names = Symbolics.tosymbol.(params)
    flux_exprs = toexprv2.(unwrap.(exprs))
    input_assign_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    compute_calls = [:($o = $expr) for (o, expr) in zip(output_names, flux_exprs)]
    return_calls = :(return [$(output_names...)])
    flux_func_expr = :(function (inputs, pas)
        $(input_assign_calls...)
        $(params_assign_calls...)
        $(compute_calls...)
        $(return_calls)
    end)
    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    return generated_flux_func
end

"""
    build_ele_func(fluxes::Vector{<:AbstractFlux}, dfluxes::Vector{<:AbstractStateFlux}, meta::ComponentVector)

Builds runtime-generated functions for both flux calculations and state differentials in a hydrological model element.

# Arguments
- `fluxes::Vector{<:AbstractFlux}`: Vector of flux components that define the element's behavior
- `dfluxes::Vector{<:AbstractStateFlux}`: Vector of state differential components that define state changes
- `meta::ComponentVector`: Metadata containing:
  - `inputs`: Input variable names
  - `outputs`: Output variable names
  - `states`: State variable names
  - `params`: Parameter names

# Returns
A tuple containing:
- First element: `Vector{Function}` with two functions:
  1. Regular flux function `(inputs, states, pas) -> outputs`
  2. Multi-dimensional flux function for batch processing
- Second element: Either `nothing` (if no states) or `Vector{Function}` with two functions:
  1. State differential function `(inputs, states, pas) -> dstates`
  2. Multi-dimensional state differential function for batch processing

# Details
The function generates four types of runtime functions:
1. Single-sample flux computation
2. Multi-sample flux computation (batched)
3. Single-sample state differential computation (if states exist)
4. Multi-sample state differential computation (if states exist)

For neural network fluxes (`AbstractNeuralFlux`), the function handles:
- Input tensor preparation
- Neural network forward passes
- Output tensor reshaping

For regular fluxes, it directly computes using provided expressions.
```
"""
function build_ele_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    infos::NamedTuple,
)
    input_names, output_names = infos.inputs, infos.outputs
    state_names, param_names = infos.states, infos.params

    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(state_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($(nflux.infos[:nns][1]) = pas.nns.$(nflux.infos[:nns][1])) for nflux in filter(f -> f isa AbstractNeuralFlux, fluxes)]
    define_calls = reduce(vcat, [input_define_calls, state_define_calls, params_assign_calls, nn_params_assign_calls])

    # varibles definitions expressions
    state_compute_calls, multi_state_compute_calls, flux_compute_calls, multi_flux_compute_calls = [], [], [], []
    for f in fluxes
        if f isa AbstractNeuralFlux
            append!(state_compute_calls, [:($(f.infos[:nn_inputs]) = [$(get_input_names(f)...)])])
            push!(state_compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(f.infos[:nns][1]))))
            append!(state_compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i]) for (i, nm) in enumerate(get_output_names(f))])

            push!(multi_state_compute_calls, :($(f.infos[:nn_inputs]) = stack([$(get_input_names(f)...)], dims=1)))
            push!(multi_state_compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(f.infos[:nns][1]))))
            append!(multi_state_compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :]) for (i, nm) in enumerate(get_output_names(f))])

            push!(flux_compute_calls, :($(f.infos[:nn_inputs]) = stack([$(get_input_names(f)...)], dims=1)))
            push!(flux_compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(f.infos[:nns][1]))))
            append!(flux_compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :]) for (i, nm) in enumerate(get_output_names(f))])

            push!(multi_flux_compute_calls, :($(f.infos[:nn_inputs]) = stack([$(get_input_names(f)...)], dims=1)))
            push!(multi_flux_compute_calls, :($(f.infos[:nn_outputs]) = stack($(f.func).(eachslice($(f.infos[:nn_inputs]), dims=2), Ref($(f.infos[:nns][1]))), dims=2)))
            append!(multi_flux_compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :, :]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($nm = $(toexpr(expr))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(multi_state_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(multi_flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    # Create return expressions with concrete values
    return_flux = :(return [$(output_names...)])
    return_state = :(return [$(map(expr -> :($(toexpr(expr))), reduce(vcat, get_exprs.(dfluxes)))...)])
    return_multi_state = :(return stack([$(map(expr -> :($(toexprv2(unwrap(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)], dims=1))

    # Create fcuntion expression
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(flux_compute_calls...)
        $(return_flux)
    end)

    multi_flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(multi_flux_compute_calls...)
        $(return_flux)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)

    if length(state_names) > 0
        diff_func_expr = :(function (inputs, states, pas)
            $(meta_exprs...)
            $(define_calls...)
            $(state_compute_calls...)
            $(return_state)
        end)

        multi_diff_func_expr = :(function (inputs, states, pas)
            $(meta_exprs...)
            $(define_calls...)
            $(multi_state_compute_calls...)
            $(return_multi_state)
        end)

        generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
        generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
        return [generated_flux_func, generated_multi_flux_func], [generated_diff_func, generated_multi_diff_func]
    else
        return [generated_flux_func, generated_multi_flux_func], nothing
    end
end

"""
    build_route_func(fluxes::AbstractVector{<:AbstractFlux}, dfluxes::AbstractVector{<:AbstractStateFlux}, meta::ComponentVector)

Builds runtime-generated functions for routing calculations in a hydrological model, handling both flux computations and state updates with special support for outflow tracking.

# Arguments
- `fluxes::AbstractVector{<:AbstractFlux}`: Vector of flux components that define the routing behavior
- `dfluxes::AbstractVector{<:AbstractStateFlux}`: Vector of state differential components that define state changes
- `meta::ComponentVector`: Metadata containing:
  - `inputs`: Input variable names
  - `outputs`: Output variable names
  - `states`: State variable names
  - `params`: Parameter names

# Returns
A tuple containing two runtime-generated functions:
1. `flux_func(inputs, states, pas) -> outputs`: Computes routing fluxes
2. `diff_func(inputs, states, pas) -> (dstates, outflows)`: Computes state changes and tracks outflows

# Details
The function specializes in routing calculations by:
1. Handling both regular and neural network-based flux components
2. Supporting batch processing for neural network operations
3. Managing tensor transformations for multi-dimensional routing
4. Tracking outflows separately from other state changes

For neural network fluxes (`AbstractNeuralFlux`), the function:
- Prepares input tensors with appropriate dimensions
- Handles batched neural network forward passes
- Reshapes outputs for compatibility with the routing system

For regular fluxes:
- Directly computes expressions for both states and fluxes
- Maintains dimensional consistency with neural network outputs
"""
function build_route_func(
    fluxes::AbstractVector{<:AbstractFlux},
    dfluxes::AbstractVector{<:AbstractStateFlux},
    infos::NamedTuple,
)
    input_names, output_names = tosymbol.(infos.inputs), tosymbol.(infos.outputs)
    state_names, param_names = tosymbol.(infos.states), tosymbol.(infos.params)

    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(state_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($nn = pas.nns.$nn) for nn in [nflux.infos[:nns][1] for nflux in filter(f -> f isa AbstractNeuralFlux, fluxes)]]
    define_calls = reduce(vcat, [input_define_calls, state_define_calls, params_assign_calls, nn_params_assign_calls])
    state_compute_calls, flux_compute_calls, = [], []
    for f in fluxes
        if f isa AbstractNeuralFlux
            push!(state_compute_calls, :($(f.infos[:nn_inputs]) = stack([$(get_input_names(f)...)], dims=1)))
            push!(state_compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(f.infos[:nns]))))
            append!(state_compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :]) for (i, nm) in enumerate(get_output_names(f))])

            push!(flux_compute_calls, :($(f.infos[:nn_inputs]) = stack([$(get_input_names(f)...)], dims=1)))
            push!(flux_compute_calls, :($(f.infos[:nn_outputs]) = stack($(f.func).(eachslice($(f.infos[:nn_inputs]), dims=2), Ref($(f.infos[:nns][1]))), dims=2)))
            append!(flux_compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :, :]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($(nm) = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($(nm) = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    dfluxes_outflows = reduce(vcat, [dflux.infos.outflows for dflux in dfluxes])
    return_state = :(return [$(map(expr -> :($(toexprv2(unwrap(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)], [$(dfluxes_outflows...)])
    # Create function expression
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(flux_compute_calls...)
        $(:(return [$(output_names...)]))
    end)

    diff_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(state_compute_calls...)
        $(return_state)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
    return generated_flux_func, generated_diff_func
end
