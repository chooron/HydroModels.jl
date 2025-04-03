struct HydroNNLayer{N} <: AbstractHydroNNLayer
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
        input_names, output_names, state_names = get_var_names(fluxes, dfluxes)
        param_names = reduce(union, get_param_names.(vcat(fluxes, dfluxes)))
        nn_names = reduce(union, get_nn_names.(fluxes))
        infos = (; inputs=input_names, outputs=output_names, states=state_names, params=param_names, nns=nn_names)
        #* Construct a function for ordinary differential calculation based on dfunc and funcs
        layer_func = build_layer_func(fluxes, dfluxes, infos)
        layer_name = isnothing(name) ? Symbol("##layer#", hash(infos)) : name
        return new{layer_name,!isempty(state_names)}(fluxes, dfluxes, layer_func, infos)
    end
end

function LuxCore.initialparameters(rng::AbstractRNG, layer::HydroNNLayer)
    init_params = NamedTuple{Tuple{get_param_names(layer)}}(map(x -> rand(rng, x), get_param_names(layer)))
    init_nns = NamedTuple{Tuple{get_nn_names(layer)}}(map(x -> rand(rng, x), get_nn_names(layer)))
    return ComponentVector(params=init_params, nns=init_nns)
end

function LuxCore.initialstates(rng::AbstractRNG, layer::HydroNNLayer)
    init_states = NamedTuple{Tuple{get_state_names(layer)}}(map(x -> rand(rng, x), get_state_names(layer)))
    return init_states
end

function (layer::HydroNNLayer)(input::AbstractArray{T,2}, params::ComponentVector, states::NamedTuple)
    return layer.layer_func(eachslice(input, dims=1), params, [states[nm] for nm in get_state_names(layer)])
end

struct HydroNNModel{N} <: AbstractHydroNNModel
    recur_layers::Vector{<:AbstractHydroNNLayer}
    "fc-layer的几种情况: (1) 接入一个dense然后整合不同节点的计算结果,(2) 接入一个UH然后求和每个节点的汇流结果"
    fc_layer::AbstractLuxLayer
    recur_op::Function
    fc_idx::Int
    infos::NamedTuple

    function HydroNNModel(;
        recur_layers::Vector{<:AbstractHydroNNLayer},
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

function _build_recur_op(layers::Vector{<:AbstractHydroNNLayer}, infos::NamedTuple)
    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(infos.inputs)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(infos.states)]
    state_update_names = Symbol.(infos.states, :_)
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

function (nn::HydroNNModel)(input::AbstractArray{T,2}, params::ComponentVector, st::NamedTuple; initstates::NamedTuple)
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