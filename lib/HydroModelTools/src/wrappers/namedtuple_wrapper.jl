struct NamedTuplePreprocessor{N} <: AbstractDataPreprocessor
    infos::NamedTuple

    function NamedTuplePreprocessor(component::C; name::Union{Symbol,Nothing}=nothing) where {C}
        wrapper_name = isnothing(name) ? Symbol("##ntp-pre#", get_name(component)) : name
        new{wrapper_name}(component.infos)
    end
end

function (adapter::NamedTuplePreprocessor)(input::NamedTuple, params, kwargs)
    input_matrix = Matrix(permutedims(reduce(hcat, [input[k] for k in get_input_names(adapter)])))
    return (input_matrix, params, kwargs)
end

function (adapter::NamedTuplePreprocessor)(input::Vector{<:NamedTuple}, params, kwargs)
    input_mats = [reduce(hcat, collect(input[i][k] for k in get_input_names(adapter))) for i in eachindex(input)]
    input_arr = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), input_mats), (2, 3, 1))
    return (input_arr, params, kwargs)
end

struct NamedTuplePostprocessor{N} <: AbstractDataPostprocessor
    infos::NamedTuple

    function NamedTuplePostprocessor(component::C; name::Union{Symbol,Nothing}=nothing) where {C}
        wrapper_name = isnothing(name) ? Symbol("##ntp-post#", get_name(component)) : name
        new{wrapper_name}(component.infos)
    end
end

function (adapter::NamedTuplePostprocessor)(output::AbstractArray)
    output_names_tuple = Tuple(vcat(get_state_names(adapter), get_output_names(adapter)))
    if length(size(output)) == 2
        return NamedTuple{output_names_tuple}(eachslice(output, dims=1))
    elseif length(size(output)) == 3
        return [NamedTuple{output_names_tuple}(eachslice(output_, dims=1)) for output_ in eachslice(output, dims=2)]
    else
        throw(ArgumentError("Output must be a 2D or 3D array"))
    end
end