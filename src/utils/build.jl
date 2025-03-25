function build_flux_func(
    inputs::Vector{Num},
    outputs::Vector{Num},
    params::Vector{Num},
    exprs::Vector{Num},
)
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

# function build_multi_flux_func(
#     inputs::Vector{Num},
#     outputs::Vector{Num},
#     params::Vector{Num},
#     exprs::Vector{Num},
# )
#     input_names, output_names = Symbolics.tosymbol.(inputs), Symbolics.tosymbol.(outputs)
#     param_names = Symbolics.tosymbol.(params)
#     flux_exprs = toexprv2.(unwrap.(exprs))
#     input_assign_calls = [:($i = inputs[$idx, :, :]) for (idx, i) in enumerate(input_names)]
#     params_assign_calls = [:($p = pas.params.$p) for p in param_names]
#     compute_calls = [:($o = $expr) for (o, expr) in zip(output_names, flux_exprs)]
#     return_calls = :(return [$(output_names...)])
#     flux_func_expr = :(function (inputs, pas)
#         $(input_assign_calls...)
#         $(params_assign_calls...)
#         $(compute_calls...)
#         $(return_calls)
#     end)
#     generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
#     return generated_flux_func
# end

function build_single_ele_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    meta::ComponentVector,
)
    input_names, output_names = tosymbol.(meta.inputs), tosymbol.(meta.outputs)
    state_names, param_names = tosymbol.(meta.states), tosymbol.(meta.params)
    nfluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)

    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($nn = pas.nns.$nn) for nn in [nflux.nninfos[:params] for nflux in nfluxes]]

    # varibles definitions expressions
    diff_define_calls = reduce(vcat, [
        [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)],
        [:($s = states[$idx]) for (idx, s) in enumerate(state_names)],
        [:($(f.nninfos[:inputs]) = [$(get_input_names(f)...)]) for f in nfluxes],
        params_assign_calls, nn_params_assign_calls,
    ])

    flux_define_calls = reduce(vcat, [
        [:($i = inputs[$idx, :]) for (idx, i) in enumerate(input_names)],
        [:($s = states[$idx, :]) for (idx, s) in enumerate(state_names)],
        [:($(f.nninfos[:inputs]) = permutedims(reduce(hcat, [$(get_input_names(f)...)]))) for f in nfluxes],
        params_assign_calls, nn_params_assign_calls,
    ])

    state_compute_calls, flux_compute_calls = [], []
    for f in fluxes
        if f isa AbstractNeuralFlux
            push!(state_compute_calls, :($(f.nninfos[:outputs]) = f.func($(f.nninfos[:inputs]), $(f.nninfos[:params]))))
            append!(state_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[i]) for (i, nm) in enumerate(get_output_names(f))])
            push!(flux_compute_calls, :($(f.nninfos[:outputs]) = f.func($(f.nninfos[:inputs]), $(f.nninfos[:params]))))
            append!(flux_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[i]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($nm = $(toexpr(expr))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    # Create return expressions with concrete values
    return_flux = :(return [$(output_names...)])
    return_state = :(return [$(map(expr -> :($(toexpr(expr))), reduce(vcat, get_exprs.(dfluxes)))...)])

    # Add meta information for optimization
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(flux_define_calls...)
        $(flux_compute_calls...)
        $(return_flux)
    end)

    diff_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(diff_define_calls...)
        $(state_compute_calls...)
        $(return_state)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
    return generated_flux_func, generated_diff_func
end


function build_multi_ele_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    meta::ComponentVector,
)
    input_names, output_names = tosymbol.(meta.inputs), tosymbol.(meta.outputs)
    state_names, param_names = tosymbol.(meta.states), tosymbol.(meta.params)
    nfluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)

    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($nn = pas.nns.$nn) for nn in [nflux.nninfos[:params] for nflux in nfluxes]]

    # varibles definitions expressions
    diff_define_calls = reduce(vcat, [
        [:($i = inputs[$idx, :]) for (idx, i) in enumerate(input_names)],
        [:($s = states[$idx, :]) for (idx, s) in enumerate(state_names)],
        [:($(f.nninfos[:inputs]) = permutedims(reduce(hcat, [$(get_input_names(f)...)]))) for f in nfluxes],
        params_assign_calls, nn_params_assign_calls,
    ])

    flux_define_calls = reduce(vcat, [
        [:($i = inputs[$idx, :, :]) for (idx, i) in enumerate(input_names)],
        [:($s = states[$idx, :, :]) for (idx, s) in enumerate(state_names)],
        [:($(f.nninfos[:inputs]) = permutedims(reduce(hcat, [$(get_input_names(f)...)]))) for f in nfluxes],
        params_assign_calls, nn_params_assign_calls,
    ])

    state_compute_calls, flux_compute_calls = [], []
    for f in fluxes
        if f isa AbstractNeuralFlux
            push!(state_compute_calls, :($(f.nninfos[:outputs]) = f.func($(f.nninfos[:inputs]), $(f.nninfos[:params]))))
            append!(state_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[i]) for (i, nm) in enumerate(get_output_names(f))])
            push!(flux_compute_calls, :($(f.nninfos[:outputs]) = f.func($(f.nninfos[:inputs]), $(f.nninfos[:params]))))
            append!(flux_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[i]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    # Create return expressions with concrete values
    return_flux = :(return [$(output_names...)])
    return_state = :(return [$(map(expr -> :($(toexprv2(unwrap(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)])

    # Add meta information for optimization
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(flux_define_calls...)
        $(flux_compute_calls...)
        $(return_flux)
    end)

    diff_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(diff_define_calls...)
        $(state_compute_calls...)
        $(return_state)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
    return generated_flux_func, generated_diff_func
end

function build_route_func(
    func::AbstractFlux,
    states::Vector{Num},
)
    inputs = setdiff(get_input_vars(func), states)
    func_params = get_param_vars(func)
    func_nns = get_nn_vars(func)
    assign_list = []
    output_list = []
    func_nns_bounds = []

    if func isa AbstractNeuralFlux
        #* For NeuralFlux, use nn_input to match the input variable
        push!(assign_list, Assignment(func.nninfos[:inputs], MakeArray(func.inputs, Vector)))
        #* For NeuralFlux, use nn_output to match the calculation result of nn expr
        push!(assign_list, Assignment(func.nninfos[:outputs], get_exprs(func)[1]))
        #* According to the output variable name, match each index of nn_output
        for (idx, output) in enumerate(get_output_vars(func))
            push!(assign_list, Assignment(output, func.nninfos[:outputs][idx]))
            push!(output_list, output)
        end
        funcs_nns_len = length.(funcs_nns)
        start_indices = [1; cumsum(funcs_nns_len)[1:end-1] .+ 1]
        end_indices = cumsum(funcs_nns_len)
        func_nns_bounds = [start:stop for (start, stop) in zip(start_indices, end_indices)]
    else
        #* According to the output variable name, match each result of the flux exprs
        for (output, expr) in zip(get_output_vars(func), get_exprs(func))
            push!(assign_list, Assignment(output, expr))
            push!(output_list, output)
        end
    end

    func_args = [
        DestructuredArgs(inputs, :inputs, inbounds=true),
        #* argument 1: Function calculation parameters
        DestructuredArgs(states, :states, inbounds=true),
        #* argument 2: Function calculation parameters
        DestructuredArgs(func_params, :params, inbounds=true),
        #* argument 3: Function neuralnetwork parameters
        DestructuredArgs(func_nns, :nns, inds=func_nns_bounds, inbounds=true),
    ]

    outputs_arr = MakeArray(output_list, Vector)
    func_bodies = Let(assign_list, outputs_arr, false)
    other_setting = [Expr(:meta, :view), Expr(:meta, :inline)]
    route_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], func_bodies, other_setting))
    )
    route_func
end
