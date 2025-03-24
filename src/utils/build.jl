function build_flux_func(
    inputs::Vector{Num},
    outputs::Vector{Num},
    params::Vector{Num},
    exprs::Vector{Num},
)
    input_names = Symbolics.tosymbol.(inputs)
    output_names = Symbolics.tosymbol.(outputs)
    param_names = Symbolics.tosymbol.(params)
    flux_exprs = toexpr.(exprs)
    input_assign_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    params_assign_calls = [:($p = params[$idx]) for (idx, p) in enumerate(param_names)]
    compute_calls = [:($target = @. $expr) for (target, expr) in zip(output_names, flux_exprs)]
    return_calls = :(return [$(output_names...)])
    flux_func_expr = :(function (inputs, params)
        $(input_assign_calls...)
        $(params_assign_calls...)
        $(compute_calls...)
        $(return_calls)
    end)
    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    return generated_flux_func
end

function build_ele_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    meta::ComponentVector,
)
    input_names, output_names = tosymbol.(meta.inputs), tosymbol.(meta.outputs)
    state_names, param_names = tosymbol.(meta.states), tosymbol.(meta.params)
    nfluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)
    funcs_nns_bounds = []
    #* prepare nn bounds
    if !isempty(nfluxes)
        funcs_nns_len = [length(nflux.nninfos[:params]) for nflux in nflux]
        start_indices = [1; cumsum(funcs_nns_len)[1:end-1] .+ 1]
        end_indices = cumsum(funcs_nns_len)
        funcs_nns_bounds = [start:stop for (start, stop) in zip(start_indices, end_indices)]
    end

    params_assign_calls = [:($p = params[$idx]) for (idx, p) in enumerate(param_names)]
    nn_params_assign_calls = [:($p = nn_params[$idx]) for (idx, p) in zip(funcs_nns_bounds, [nflux.nninfos[:params] for nflux in nfluxes])]

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

    state_compute_calls = []
    flux_compute_calls = []
    for f in fluxes
        if f isa AbstractNeuralFlux
            push!(state_compute_calls, :($(f.nninfos[:outputs]) = f.func($(f.nninfos[:inputs]), $(f.nninfos[:params]))))
            append!(state_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[i]) for (i, nm) in enumerate(get_output_names(f))])
            push!(flux_compute_calls, :($(f.nninfos[:outputs]) = f.func($(f.nninfos[:inputs]), $(f.nninfos[:params]))))
            append!(flux_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[i]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($nm = $(toexpr(expr))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr), LazyState()))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    # Create return expressions with concrete values
    return_flux = :(return [$(output_names...)])
    return_state = :(return [$(map(expr -> :($(toexpr(expr))), reduce(vcat, [dflux.exprs for dflux in dfluxes]))...)])

    # Add meta information for optimization
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, params, nn_params)
        $(meta_exprs...)
        $(flux_define_calls...)
        $(flux_compute_calls...)
        $(return_flux)
    end)

    diff_func_expr = :(function (inputs, states, params, nn_params)
        $(meta_exprs...)
        $(diff_define_calls...)
        $(state_compute_calls...)
        $(return_state)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
    return generated_flux_func, generated_diff_func
end


#* 构建bucket的函数
function build_ele_func_old(
    funcs::Vector{<:AbstractFlux},
    dfuncs::Vector{<:AbstractStateFlux},
    meta::ComponentVector,
)
    #* prepare variables for function building
    funcs_inputs, funcs_states, funcs_params = meta.inputs, meta.states, meta.params
    #* prepare Assignment and Output Variables
    assign_list, output_list, funcs_nns = Assignment[], Num[], []
    for func in funcs
        if func isa AbstractNeuralFlux
            push!(funcs_nns, func.nninfos[:nns])
            #* For NeuralFlux, use nn_input to match the input variable
            push!(assign_list, Assignment(func.nninfos[:inputs], MakeArray(get_input_vars(func), Vector)))
            #* For NeuralFlux, use nn_output to match the calculation result of nn expr
            push!(assign_list, Assignment(func.nninfos[:outputs], get_exprs(func)[1]))
            #* According to the output variable name, match each index of nn_output
            for (idx, output) in enumerate(get_output_vars(func))
                push!(assign_list, Assignment(output, func.nninfos[:outputs][idx]))
                push!(output_list, output)
            end
        else
            #* According to the output variable name, match each result of the flux exprs
            for (output, expr) in zip(get_output_vars(func), get_exprs(func))
                push!(assign_list, Assignment(output, expr))
                push!(output_list, output)
            end
        end
    end

    funcs_nns_bounds = []
    #* prepare nn bounds
    if !isempty(funcs_nns)
        funcs_nns_len = [length(funcs_nns[k]) for k in keys(funcs_nns)]
        start_indices = [1; cumsum(funcs_nns_len)[1:end-1] .+ 1]
        end_indices = cumsum(funcs_nns_len)
        funcs_nns_bounds = [start:stop for (start, stop) in zip(start_indices, end_indices)]
    end

    #* convert flux output to array (include states and outputs)
    flux_output_array = MakeArray(vcat(funcs_states, output_list), Vector)

    #* Set the input argument of ODE Function
    func_args = [
        #* argument 1: Function input variables
        DestructuredArgs(funcs_inputs, :inputs, inbounds=true),
        #* argument 2: Function state variables
        DestructuredArgs(funcs_states, :states, inbounds=true),
        #* argument 2: Function calculation parameters
        DestructuredArgs(funcs_params, :params, inbounds=true),
        #* argument 3: Function neuralnetwork parameters
        DestructuredArgs(funcs_nns, :nns, inds=funcs_nns_bounds, inbounds=true),
    ]

    #* Construct Flux Function: Func(args, kwargs, body), where body represents the matching formula between each variable and expression
    generated_flux_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, flux_output_array, false)))
    )

    if !isempty(dfuncs)
        #* convert diff state output to array
        diffst_output_array = MakeArray(reduce(vcat, get_exprs.(dfuncs)), Vector)
        generated_diff_func = @RuntimeGeneratedFunction(
            toexpr(Func(func_args, [], Let(assign_list, diffst_output_array, false)))
        )
    else
        generated_diff_func = nothing
    end
    generated_flux_func, generated_diff_func
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
