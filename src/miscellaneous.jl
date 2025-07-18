using Symbolics

struct FunctionElement <: AbstractElement
    inputs::Vector{Num}
    func::Function
end

function (ele::FunctionElement)(input::AbstractArray, params::ComponentVector; kwargs...) where {T}
    ele.func(input, params; kwargs...)
end

function SummationElement(inputs::Vector{Num})
    return FunctionElement(inputs, (input, params) -> sum(input, dims=2))
end

"""
# GroupHydroFlux.jl
构建全是由HydroFlux构建的组合Flux
"""
function GroupHydroFlux(
    pairs::Vector{Pair{Num,HF}};
    name::Union{Symbol,Nothing}=nothing
) where HF<:AbstractHydroFlux
    weights = [f[1] for f in pairs]  # Extract the fluxes from the pairs
    fluxes = [f[2] for f in pairs]  # Extract the fluxes from the pairs

    @assert all(f -> f isa HydroFlux, fluxes) "All elements in fluxes must be of type AbstractHydroFlux."
    @assert all(length(f.infos.outputs) == 1 for f in fluxes) "Each flux only has one output variable."
    @assert allequal([f.infos.outputs[1] for f in fluxes]) "All fluxes must have the same output names."

    # build group expression
    outputs = fluxes[1].infos.outputs
    group_exprs = reduce(+, map(zip(fluxes, weights)) do (f, w)
        f.exprs[1] * w
    end)

    group_inputs = unique(reduce(vcat, [f.infos.inputs for f in fluxes]))
    weight_params = unique(reduce(vcat, map(weights) do weight
        HydroModels.Num.([w for w in HydroModels.get_variables(weight)])
    end))
    group_params = vcat(unique(reduce(vcat, [f.infos.params for f in fluxes])), weight_params)
    flux_name = isnothing(name) ? Symbol("##group_flux#", hash(group_params)) : name

    return HydroFlux(group_inputs, outputs, group_params; exprs=[group_exprs], name=flux_name)
end

macro grouphydroflux(args...)
    local name_arg, block_expr
    if length(args) == 1
        name_arg = :nothing
        block_expr = args[1]
    elseif length(args) == 2
        name_arg = QuoteNode(args[1])
        block_expr = args[2]
    else
        return :(error("@grouphydroflux expects an optional name and a begin...end block."))
    end
    if !(block_expr isa Expr && block_expr.head == :block)
        return :(error("The last argument to @grouphydroflux must be a begin...end block."))
    end
    pair_exprs = filter(x -> !(x isa LineNumberNode), block_expr.args)
    pairs_vector_expr = Expr(:vect, pair_exprs...)
    return esc(quote
        GroupHydroFlux($pairs_vector_expr; name=$name_arg)
    end)
end

"""
    GroupStateFluxes
    stateflux可能会涉及多个状态变量: 
"""
function GroupStateFlux(
    pairs::Vector{<:Pair};
    name::Union{Vector{Symbol},Nothing}=nothing
)
    weights = [f[1] for f in pairs]
    fluxes_vec = [f[2] for f in pairs]
    max_flux_num = maximum(length.(fluxes_vec))

    group_fluxes_vec = map(1:max_flux_num) do i
        [fluxes[i] for fluxes in fluxes_vec if i <= length(fluxes)]
    end

    states = [fluxes[1].infos.states for fluxes in group_fluxes_vec]
    group_inputs = [unique(reduce(vcat, [f.infos.inputs for f in fluxes])) for fluxes in group_fluxes_vec]
    group_params = [unique(reduce(vcat, [f.infos.params for f in fluxes])) for fluxes in group_fluxes_vec]

    group_exprs = []
    for i in eachindex(group_fluxes_vec)
        push!(group_exprs, reduce(+, map(zip(group_fluxes_vec[i], weights)) do (f, w)
            f.exprs[1] * w
        end))
    end

    return map(eachindex(group_fluxes_vec)) do i
        flux_name = isnothing(name) ? Symbol("##group_flux#", hash(group_exprs[i])) : name
        StateFlux(group_inputs[i], states[i][1], group_params[i]; expr=group_exprs[i], name=flux_name)
    end
end

macro groupstateflux(args...)
    local name_arg, block_expr
    if length(args) == 1
        name_arg = :nothing
        block_expr = args[1]
    elseif length(args) == 2
        name_arg = QuoteNode(args[1])
        block_expr = args[2]
    else
        return :(error("@groupstateflux expects an optional name and a begin...end block."))
    end
    if !(block_expr isa Expr && block_expr.head == :block)
        return :(error("The last argument to @groupstateflux must be a begin...end block."))
    end
    pair_exprs = filter(x -> !(x isa LineNumberNode), block_expr.args)
    pairs_vector_expr = Expr(:vect, pair_exprs...)
    return esc(quote
        GroupStateFlux($pairs_vector_expr; name=$name_arg)
    end)
end