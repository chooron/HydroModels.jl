struct NeuralWrapper2{N,M} <: AbstractNeuralWrapper
    model
    meta::M

    function NeuralWrapper2(fluxes::Pair, model::N; model_name::Symbol, name::Union{Symbol,Nothing}=nothing, rng::R=StableRNG(42)) where {N,R}
        ps, st = Lux.setup(rng, model)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func(x, p) = LuxCore.apply(model, x, ComponentVector(p, ps_axes), st)[1]
        nn_var = first(@parameters $(model_name)[1:length(ps)])
        meta = ComponentVector(inputs=fluxes[1], outputs=fluxes[2], nns=NamedTuple{Tuple([model_name])}([nn_var]))
        name = isnothing(name) ? Symbol("##wrapper#", hash(meta)) : name
        return new{name,typeof(meta)}(model, meta)
    end
end

function (wrapper::NeuralWrapper2)(input::Array{T,2}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nns][get_nn_names(wrapper)]
    output = wrapper.func(input, nn_ps)
    return output
end

function (wrapper::NeuralWrapper2)(input::Array{T,3}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nns][get_nn_names(wrapper)]
    output = wrapper.func.(Matrix.(eachslice(input, dims=2)), Ref(nn_ps))
    return permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), output), (1, 3, 2))
end

struct NeuralWrapper3{N,M} <: AbstractNeuralWrapper
    model::N
    meta::M

    function NeuralWrapper3(fluxes::Pair, model::N;model_name::Symbol, name::Union{Symbol,Nothing}=nothing, rng::R=StableRNG(42)) where {N,R}
        ps, st = Lux.setup(rng, model)
        ps_axes = getaxes(ComponentVector(ps))
        nn_func(x, p) = LuxCore.apply(model, x, ComponentVector(p, ps_axes), st)[1]
        meta = ComponentVector(inputs=fluxes[1], outputs=fluxes[2], nns=NamedTuple{Tuple([model_name])}([nn_var]))
        name = isnothing(name) ? Symbol("##wrapper#", hash(meta)) : name
        return new{name,typeof(meta)}(model, meta)
    end
end

function (wrapper::NeuralWrapper3)(input::Array{T,3}, pas::ComponentVector; kwargs...) where {T}
    nn_ps = pas[:nns][get_nn_names(wrapper)]
    output = wrapper.func(input, Ref(nn_ps))
    return output
end
