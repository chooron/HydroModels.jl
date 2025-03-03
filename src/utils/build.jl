function build_flux_func(
    inputs::Vector{Num},
    outputs::Vector{Num},
    params::Vector{Num},
    exprs::Vector{Num},
)
    assign_list = Assignment.(outputs, exprs)
    outputs_arr = MakeArray(outputs, Vector)
    func_args = [
        #* argument 1: Function calculation parameters
        DestructuredArgs(inputs),
        #* argument 2: Function neuralnetwork parameters
        DestructuredArgs(params),
    ]
    call_func = @RuntimeGeneratedFunction(
        toexpr(Func(func_args, [], Let(assign_list, outputs_arr, false)))
    )
    call_func
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

#* 构建bucket的函数
function build_ele_func(
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

