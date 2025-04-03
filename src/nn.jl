struct HydroNNLayer{N} <: AbstractNNLayer
    "hydrological fluxes"
    fluxes::Vector{<:AbstractHydroFlux}
    "hydrological state derivatives"
    dfluxes::Vector{<:AbstractStateFlux}
    "layer function"
    layer_func::Function
    "meta data of hydrological nn layer"
    infos::NamedTuple

    function HydroNNLayer(;
        name::Union{Symbol,Nothing}=nothing,
        fluxes::Vector{<:AbstractHydroFlux},
        dfluxes::Vector{<:AbstractStateFlux}=StateFlux[],
        sort_fluxes::Bool=false,
    )
        #* sort the fluxes if needed
        fluxes = sort_fluxes ? sort_fluxes(fluxes) : fluxes
        #* Extract all variable names of fluxes and dfluxes
        input_names, output_names, state_names = get_vars(fluxes, dfluxes)
        param_names = reduce(union, get_param_names.(vcat(fluxes, dfluxes)))
        nn_names = reduce(union, get_nn_names.(fluxes))
        infos = (; inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        layer_func = _build_layer_func(fluxes, dfluxes, infos)
        layer_name = isnothing(name) ? Symbol("##layer#", hash(infos)) : name
        return new{layer_name}(fluxes, dfluxes, layer_func, infos)
    end
end

function _build_layer_func(fluxes::Vector{<:AbstractHydroFlux}, dfluxes::Vector{<:AbstractStateFlux}, infos::NamedTuple)
    input_names = length(infos.inputs) == 0 ? [] : tosymbol.(infos.inputs)
    output_names = length(infos.outputs) == 0 ? [] : tosymbol.(infos.outputs)
    state_names = length(infos.states) == 0 ? [] : tosymbol.(infos.states)
    param_names = length(infos.params) == 0 ? [] : tosymbol.(infos.params)

    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(state_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($(nflux.infos[:nns][1]) = pas.nns.$(nflux.infos[:nns][1])) for nflux in filter(f -> f isa AbstractNeuralFlux, fluxes)]
    define_calls = reduce(vcat, [input_define_calls, state_define_calls, params_assign_calls, nn_params_assign_calls])

    # varibles definitions expressions
    compute_calls = []
    for f in fluxes
        if f isa AbstractNeuralFlux
            push!(compute_calls, :($(f.infos[:nn_inputs]) = stack([$(get_input_names(f)...)], dims=1)))
            push!(compute_calls, :($(f.infos[:nn_outputs]) = $(f.func)($(f.infos[:nn_inputs]), $(f.infos[:nns][1]))))
            append!(compute_calls, [:($(nm) = $(f.infos[:nn_outputs])[$i, :]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    # Create return expressions with concrete values
    return_flux = :(return tuple($(output_names...)), tuple($(map(expr -> :($(toexprv2(unwrap(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)))

    # Create fcuntion expression
    meta_exprs = [:(Base.@_inline_meta)]

    func_expr = :(function (inputs, pas, states)
        $(meta_exprs...)
        $(define_calls...)
        $(compute_calls...)
        $(return_flux)
    end)


    println(func_expr)

    return @RuntimeGeneratedFunction(func_expr)
end

function LuxCore.initialparameters(rng::AbstractRNG, layer::HydroNNLayer)
    init_params = NamedTuple{Tuple(get_param_names(layer))}(map(get_params(layer)) do p
        lb, ub = getbounds(p)
        rand(rng, layer.hidden_dims) .* (ub - lb) .+ lb
    end)
    init_nns = NamedTuple{Tuple(get_nn_names(layer))}(map(get_nn_names(layer)) do nn
        glorot_normal(rng, length(nn))
    end)
    return ComponentVector(params=init_params, nns=init_nns)
end

function LuxCore.initialstates(rng::AbstractRNG, layer::HydroNNLayer)
    init_states = NamedTuple{Tuple(get_state_names(layer))}(map(get_states(layer)) do s
        zeros(layer.hidden_dims)
    end)
    return init_states
end

function (layer::HydroNNLayer)(input::AbstractArray{T,2}, params::ComponentVector, states::NamedTuple) where {T}
    return layer.layer_func(eachslice(input, dims=1), params, [states[nm] for nm in get_state_names(layer)])
end

struct HydroNNModel{N} <: AbstractNNModel
    recur_layers::Vector{<:AbstractNNLayer}
    "fc-layer的几种情况: (1) 接入一个dense然后整合不同节点的计算结果,(2) 接入一个UH然后求和每个节点的汇流结果"
    fc_layer::AbstractLuxLayer
    recur_op::Function
    fc_idx::Int
    infos::NamedTuple

    function HydroNNModel(;
        recur_layers::Vector{<:AbstractNNLayer},
        fc_layer::AbstractLuxLayer,
        name::Union{Symbol,Nothing}=nothing,
        fc_variables::Union{AbstractArray{<:Symbol},Nothing}=nothing,
    )
        input_names, output_names, state_names = get_var_names(recur_layers)
        param_names = reduce(union, get_param_names.(recur_layers))
        nn_names = reduce(union, get_nn_names.(recur_layers))
        fc_idx = findfirst(x -> x == fc_variables, output_names)

        infos = (; inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        recur_op = _build_recur_op(recur_layers, infos)
        model_name = isnothing(name) ? Symbol("##model#", hash(infos)) : name
        new{model_name}(recur_layers, fc_layer, recur_op, fc_idx, infos)
    end
end

function _build_recur_op(layers::Vector{<:AbstractNNLayer}, infos::NamedTuple)
    input_names = length(infos.inputs) == 0 ? [] : tosymbol.(infos.inputs)
    state_names = length(infos.states) == 0 ? [] : tosymbol.(infos.states)
    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(state_names)]
    state_update_names = Symbol.(state_names, :_)
    compute_calls = map(layers) do layer
        i, s, o = get_var_names(layer)
        s′ = Symbol.(s, :_)
        :(((($o...), ($s′...)), st) = $layer(([($i...)], [($s...)]), params, st))
    end
    return_calls = :((stack([$infos.outputs...], dims=1), ($state_update_names...)), st)
    return quote
        function (inputs, states, params, st)
            $input_define_calls
            $state_define_calls
            $(compute_calls...)
            $return_calls
        end
    end
end

function (nn::HydroNNModel)(input::AbstractArray{T,3}, params::ComponentVector, st::NamedTuple; initstates::NamedTuple) where {T}
    function recur_op(::Nothing, input)
        (out, carry), st′ = nn.recur_op((input, initstates), params, st)
        return [out], carry, st′
    end
    function recur_op((outputs, carry, state), input)
        (out, carry′), st′ = nn.recur_op((input, carry), params, state)
        return vcat(outputs, [out]), carry′, st′
    end
    results = foldl_init(recur_op, eachslice(input, dims=3))
    output = stack(first(results), dims=3)
    fc_output, st′ = LuxCore.apply(nn.fc_layer, output[end, :, :], params.fc, last(results))
    return fc_output, st′
end