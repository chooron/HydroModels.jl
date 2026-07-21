"""
    $(TYPEDEF)

Stores attribute information for a hydrological component.

# Fields
$(FIELDS)
"""
struct HydroInfos{IS,OS,SS,PS,NS} <: AbstractInfos
    "Input variable names."
    inputs::IS
    "Output variable names."
    outputs::OS
    "State variable names."
    states::SS
    "Parameter names."
    params::PS
    "Neural network names."
    nns::NS

    function HydroInfos(;
        inputs::AbstractVector{Symbol}=Symbol[],
        outputs::AbstractVector{Symbol}=Symbol[],
        states::AbstractVector{Symbol}=Symbol[],
        params::AbstractVector{Symbol}=Symbol[],
        nns::AbstractVector{Symbol}=Symbol[]
    )
        new{typeof(inputs),typeof(outputs),typeof(states),typeof(params),typeof(nns)}(
            inputs, outputs, states, params, nns
        )
    end
end


# ================================================================================================
# 基础访问器（使用 Val 分派优化）
# ================================================================================================

"""
    $(SIGNATURES)

Get component name.
"""
@inline get_name(cpt::AbstractComponent) = cpt.name

"""
    $(SIGNATURES)

Get input names from `HydroInfos` or `AbstractComponent`.
"""
@inline get_input_names(infos::HydroInfos) = infos.inputs
@inline get_input_names(cpt::AbstractComponent) = get_input_names(cpt.infos)

"""
    $(SIGNATURES)

Get output names from `HydroInfos` or `AbstractComponent`.
"""
@inline get_output_names(infos::HydroInfos) = infos.outputs
@inline get_output_names(cpt::AbstractComponent) = get_output_names(cpt.infos)

"""
    $(SIGNATURES)

Get state names from `HydroInfos` or `AbstractComponent`.
"""
@inline get_state_names(infos::HydroInfos) = infos.states
@inline get_state_names(cpt::AbstractComponent) = get_state_names(cpt.infos)

"""
    $(SIGNATURES)

Get parameter names from `HydroInfos` or `AbstractComponent`.
"""
@inline get_param_names(infos::HydroInfos) = infos.params
@inline get_param_names(cpt::AbstractComponent) = get_param_names(cpt.infos)

"""
    $(SIGNATURES)

Get neural network names from `HydroInfos` or `AbstractComponent`.
"""
@inline get_nn_names(infos::HydroInfos) = infos.nns
@inline get_nn_names(cpt::AbstractComponent) = get_nn_names(cpt.infos)


# ================================================================================================
# 组合访问器
# ================================================================================================

"""
    $(SIGNATURES)

Get input, output, and state names from `HydroInfos` or `AbstractComponent`.
Returns a tuple `(inputs, outputs, states)`.
"""
@inline get_var_names(infos::HydroInfos) = (get_input_names(infos), get_output_names(infos), get_state_names(infos))
@inline get_var_names(cpt::AbstractComponent) = get_var_names(cpt.infos)

"""
    $(SIGNATURES)

Get all names (inputs, outputs, states, params, nns) from `HydroInfos` or `AbstractComponent`.
Returns a tuple `(inputs, outputs, states, params, nns)`.
"""
@inline get_all_names(infos::HydroInfos) = (infos.inputs, infos.outputs, infos.states, infos.params, infos.nns)
@inline get_all_names(cpt::AbstractComponent) = get_all_names(cpt.infos)

"""
    $(SIGNATURES)

Get unique variable names from a collection of components.
Returns a tuple `(inputs, outputs, states)`.

# Performance
This function is optimized to avoid unnecessary allocations when possible.
"""
function get_var_names(components::CT) where {CT}
    inputs, outputs = Vector{Symbol}(), Vector{Symbol}()
    states = reduce(union, get_state_names.(components); init=Symbol[])

    for comp in components
        tmp_inputs = get_input_names(comp)
        tmp_outputs = get_output_names(comp)

        # 只添加不在 outputs 中的 inputs
        tmp_inputs = setdiff(tmp_inputs, outputs)
        union!(inputs, tmp_inputs)
        union!(outputs, tmp_outputs)
    end

    # 从 inputs 中移除 states
    setdiff!(inputs, states)
    return inputs, outputs, states
end


# ================================================================================================
# 表达式获取器（多态分派）
# ================================================================================================

"""
    $(SIGNATURES)

Get expressions from a flux component. Dispatches based on component type.
"""
@inline get_exprs(cpt::AbstractFlux) = cpt.exprs
@inline get_exprs(cpt::AbstractNeuralFlux) = get_output_names(cpt) ~ cpt.chain


# ================================================================================================
# 数量统计器（使用 Val 分派优化）
# ================================================================================================

"""
    $(SIGNATURES)

Count the number of a specific attribute type in `HydroInfos`.
Type is specified using `Val{:symbol}` where symbol is one of `:inputs`, `:outputs`, `:states`, `:params`, `:nns`.
"""
@inline count_attrs(infos::HydroInfos, ::Val{:inputs}) = length(infos.inputs)
@inline count_attrs(infos::HydroInfos, ::Val{:outputs}) = length(infos.outputs)
@inline count_attrs(infos::HydroInfos, ::Val{:states}) = length(infos.states)
@inline count_attrs(infos::HydroInfos, ::Val{:params}) = length(infos.params)
@inline count_attrs(infos::HydroInfos, ::Val{:nns}) = length(infos.nns)

"""
    $(SIGNATURES)

Convenience wrapper for counting attributes.
"""
@inline count_inputs(x) = count_attrs(get_infos(x), Val(:inputs))
@inline count_outputs(x) = count_attrs(get_infos(x), Val(:outputs))
@inline count_states(x) = count_attrs(get_infos(x), Val(:states))
@inline count_params(x) = count_attrs(get_infos(x), Val(:params))
@inline count_nns(x) = count_attrs(get_infos(x), Val(:nns))

@inline get_infos(infos::HydroInfos) = infos
@inline get_infos(cpt::AbstractComponent) = cpt.infos


# ================================================================================================
# 存在性检查器（优化查询性能）
# ================================================================================================

"""
    $(SIGNATURES)

Check if a variable name exists in the specified attribute category.
"""
@inline has_input(infos::HydroInfos, name::Symbol) = name in infos.inputs
@inline has_output(infos::HydroInfos, name::Symbol) = name in infos.outputs
@inline has_state(infos::HydroInfos, name::Symbol) = name in infos.states
@inline has_param(infos::HydroInfos, name::Symbol) = name in infos.params
@inline has_nn(infos::HydroInfos, name::Symbol) = name in infos.nns

@inline has_input(cpt::AbstractComponent, name::Symbol) = has_input(cpt.infos, name)
@inline has_output(cpt::AbstractComponent, name::Symbol) = has_output(cpt.infos, name)
@inline has_state(cpt::AbstractComponent, name::Symbol) = has_state(cpt.infos, name)
@inline has_param(cpt::AbstractComponent, name::Symbol) = has_param(cpt.infos, name)
@inline has_nn(cpt::AbstractComponent, name::Symbol) = has_nn(cpt.infos, name)

"""
    $(SIGNATURES)

Check if any variable names exist in the specified attribute category.
"""
@inline has_any_input(infos::HydroInfos, names::AbstractVector{Symbol}) = any(n -> has_input(infos, n), names)
@inline has_any_output(infos::HydroInfos, names::AbstractVector{Symbol}) = any(n -> has_output(infos, n), names)
@inline has_any_state(infos::HydroInfos, names::AbstractVector{Symbol}) = any(n -> has_state(infos, n), names)
@inline has_any_param(infos::HydroInfos, names::AbstractVector{Symbol}) = any(n -> has_param(infos, n), names)
@inline has_any_nn(infos::HydroInfos, names::AbstractVector{Symbol}) = any(n -> has_nn(infos, n), names)


# ================================================================================================
# 空值检查器
# ================================================================================================

"""
    $(SIGNATURES)

Check if a specific attribute category is empty.
"""
@inline is_empty(infos::HydroInfos, ::Val{:inputs}) = isempty(infos.inputs)
@inline is_empty(infos::HydroInfos, ::Val{:outputs}) = isempty(infos.outputs)
@inline is_empty(infos::HydroInfos, ::Val{:states}) = isempty(infos.states)
@inline is_empty(infos::HydroInfos, ::Val{:params}) = isempty(infos.params)
@inline is_empty(infos::HydroInfos, ::Val{:nns}) = isempty(infos.nns)

"""
    $(SIGNATURES)

Convenience wrappers for checking if attribute categories are empty.
"""
@inline has_inputs(x) = !is_empty(get_infos(x), Val(:inputs))
@inline has_outputs(x) = !is_empty(get_infos(x), Val(:outputs))
@inline has_states(x) = !is_empty(get_infos(x), Val(:states))
@inline has_params(x) = !is_empty(get_infos(x), Val(:params))
@inline has_nns(x) = !is_empty(get_infos(x), Val(:nns))


# ================================================================================================
# 集合操作辅助函数
# ================================================================================================

"""
    $(SIGNATURES)

Find intersection of names between two `HydroInfos` objects for a specific attribute type.
"""
@inline intersect_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:inputs}) = intersect(info1.inputs, info2.inputs)
@inline intersect_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:outputs}) = intersect(info1.outputs, info2.outputs)
@inline intersect_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:states}) = intersect(info1.states, info2.states)
@inline intersect_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:params}) = intersect(info1.params, info2.params)
@inline intersect_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:nns}) = intersect(info1.nns, info2.nns)

"""
    $(SIGNATURES)

Find union of names between two `HydroInfos` objects for a specific attribute type.
"""
@inline union_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:inputs}) = union(info1.inputs, info2.inputs)
@inline union_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:outputs}) = union(info1.outputs, info2.outputs)
@inline union_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:states}) = union(info1.states, info2.states)
@inline union_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:params}) = union(info1.params, info2.params)
@inline union_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:nns}) = union(info1.nns, info2.nns)

"""
    $(SIGNATURES)

Find difference of names between two `HydroInfos` objects for a specific attribute type.
"""
@inline setdiff_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:inputs}) = setdiff(info1.inputs, info2.inputs)
@inline setdiff_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:outputs}) = setdiff(info1.outputs, info2.outputs)
@inline setdiff_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:states}) = setdiff(info1.states, info2.states)
@inline setdiff_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:params}) = setdiff(info1.params, info2.params)
@inline setdiff_attrs(info1::HydroInfos, info2::HydroInfos, ::Val{:nns}) = setdiff(info1.nns, info2.nns)


# ================================================================================================
# 批量操作辅助函数
# ================================================================================================

"""
    $(SIGNATURES)

Collect all unique names of a specific attribute type from a collection of components.
"""
@inline function collect_unique_attrs(components, ::Val{:inputs})
    return reduce(union, get_input_names.(components); init=Symbol[])
end

@inline function collect_unique_attrs(components, ::Val{:outputs})
    return reduce(union, get_output_names.(components); init=Symbol[])
end

@inline function collect_unique_attrs(components, ::Val{:states})
    return reduce(union, get_state_names.(components); init=Symbol[])
end

@inline function collect_unique_attrs(components, ::Val{:params})
    return reduce(union, get_param_names.(components); init=Symbol[])
end

@inline function collect_unique_attrs(components, ::Val{:nns})
    return reduce(union, get_nn_names.(components); init=Symbol[])
end

"""
    $(SIGNATURES)

Convenience wrappers for collecting unique attributes from components.
"""
@inline collect_all_inputs(components) = collect_unique_attrs(components, Val(:inputs))
@inline collect_all_outputs(components) = collect_unique_attrs(components, Val(:outputs))
@inline collect_all_states(components) = collect_unique_attrs(components, Val(:states))
@inline collect_all_params(components) = collect_unique_attrs(components, Val(:params))
@inline collect_all_nns(components) = collect_unique_attrs(components, Val(:nns))
