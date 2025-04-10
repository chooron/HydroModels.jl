
function check_input(component::AbstractComponent, input::AbstractArray{<:Number,2}, timeidx::AbstractVector)
    @assert size(input, 1) == length(get_input_names(component)) "Input variables in component '$(component.meta.name)' do not match required dimensions. Expected $(length(get_input_names(component))) variables ($(get_input_names(component))), got $(size(input, 1)) variables"
    @assert size(input, 2) == length(timeidx) "Time steps in component '$(component.meta.name)' do not match required length. Expected $(length(timeidx)) steps, got $(size(input, 2)) steps"
end

function check_input(component::AbstractComponent, input::AbstractArray{<:Number,3}, timeidx::AbstractVector)
    @assert size(input, 1) == length(get_input_names(component)) "Input variables in component '$(component.meta.name)' do not match required dimensions. Expected $(length(get_input_names(component))) variables ($(get_input_names(component))), got $(size(input, 1)) variables"
    @assert size(input, 3) == length(timeidx) "Time steps in component '$(component.meta.name)' do not match required length. Expected $(length(timeidx)) steps, got $(size(input, 3)) steps"
end

function check_ptypes(component::AbstractComponent, input::AbstractArray{<:Number,3}, ptypes::AbstractVector{Symbol})
    @assert length(ptypes) == size(input, 2) "Number of parameter types mismatch In $(get_name(component)). Expected $(size(input, 2)) parameter types, got $(length(ptypes))."
end

function check_stypes(component::AbstractComponent, input::AbstractArray{<:Number,3}, stypes::AbstractVector{Symbol})
    @assert length(stypes) == size(input, 2) "Number of state types mismatch In $(get_name(component)). Expected $(size(input, 2)) state types, got $(length(stypes))."
end

function check_pas(component::AbstractComponent, pas::ComponentVector)
    check_parameters(component, pas)
    check_initstates(component, pas)
    check_nns(component, pas)
end

function check_pas(component::AbstractComponent, pas::ComponentVector, ptypes::AbstractVector{Symbol}, stypes::AbstractVector{Symbol})
    check_parameters(component, pas, ptypes)
    check_initstates(component, pas, stypes)
    check_nns(component, pas)
end

function check_parameters(component::AbstractComponent, pas::ComponentVector)
    param_names = get_param_names(component)
    cpt_name = get_name(component)
    for param_name in param_names
        @assert(param_name in keys(pas[:params]),
            "Parameter '$(param_name)' in component '$(cpt_name)' is required but not found in pas[:params]. Available parameters: $(keys(pas[:params]))"
        )
    end
end

function check_parameters(component::AbstractComponent, pas::ComponentVector, ptypes::AbstractVector{Symbol})
    param_names = get_param_names(component)
    cpt_name = get_name(component)
    for ptype in ptypes
        tmp_ptype_params_keys = keys(pas[:params][ptype])
        for param_name in param_names
            @assert(param_name in tmp_ptype_params_keys,
                "Parameter '$(param_name)' in component '$(cpt_name)' is required but not found in parameter type '$(ptype)'. Available parameters: $(tmp_ptype_params_keys)"
            )
        end
    end
end

function check_initstates(component::AbstractComponent, pas::ComponentVector)
    state_names = get_state_names(component)
    cpt_name = get_name(component)
    for state_name in state_names
        tmp_ptype_initstates_keys = keys(pas[:initstates])
        @assert(state_name in tmp_ptype_initstates_keys,
            "Initial state '$(state_name)' in component '$(cpt_name)' is required but not found in parameter type '$(ptype)'. Available states: $(tmp_ptype_initstates_keys)"
        )
    end
end

function check_initstates(component::AbstractComponent, pas::ComponentVector, stypes::AbstractVector{Symbol})
    state_names = get_state_names(component)
    cpt_name = get_name(component)
    for stype in stypes
        tmp_ptype_initstates_keys = keys(pas[:initstates][stype])
        for state_name in state_names
            @assert(state_name in tmp_ptype_initstates_keys,
                "Initial state '$(state_name)' in component '$(cpt_name)' is required but not found in state type '$(stype)'. Available states: $(tmp_ptype_initstates_keys)"
            )
        end
    end
end

function check_nns(component::AbstractComponent, pas::ComponentVector)
    nn_names = get_nn_names(component)
    cpt_name = get_name(component)
    if !isempty(nn_names)
        for nn_name in nn_names
            @assert nn_name in keys(pas[:nn]) "Neural network parameters '$(nn_name)' in component '$(cpt_name)' is required but not found in pas[:nn]. Available networks: $(keys(pas[:nn]))"
        end
    end
end

# 处理pas，避免在后续计算中反复转换
function process_pas(c::AbstractComponent, pas::ComponentVector)
    check_pas(c, pas)
    pas_type = eltype(pas)
    new_pas = ComponentVector(
        params=haskey(pas, :params) ? Vector(pas[:params][HydroModelCore.get_param_names(c)]) : Vector{pas_type}(undef, 0),
        initstates=haskey(pas, :initstates) ? Vector(pas[:initstates][HydroModelCore.get_state_names(c)]) : Vector{pas_type}(undef, 0),
        nns=haskey(pas, :nns) ? Vector(pas[:nns][HydroModelCore.get_nn_names(c)]) : Vector{pas_type}(undef, 0),
    )
    return new_pas
end

function process_pas(c::AbstractComponent, pas::ComponentVector, ptypes::AbstractVector{Symbol}, stypes::AbstractVector{Symbol})
    check_pas(c, pas, ptypes, stypes)
    pas_type = eltype(pas)
    params_names, state_names = get_param_names(c), get_state_names(c)

    unq_ptypes = unique(ptypes)
    unq_stypes = unique(stypes)

    ptype_indices = [findfirst(isequal(ptype), unq_ptypes) for ptype in ptypes]
    stype_indices = [findfirst(isequal(stype), unq_stypes) for stype in stypes]

    params_mat = isempty(params_names) ? Matrix{pas_type}(undef, 0, 0) : reduce(hcat, [Vector(pas[:params][ptype][params_names]) for ptype in unq_ptypes])
    initstates_mat = isempty(state_names) ? Matrix{pas_type}(undef, 0, 0) : reduce(hcat, [Vector(pas[:initstates][stype][state_names]) for stype in unq_stypes])

    new_pas = ComponentVector(
        params=params_mat,
        initstates=initstates_mat,
        nns=haskey(pas, :nns) ? Vector(pas[:nns][HydroModelCore.get_nn_names(c)]) : Vector{pas_type}(undef, 0),
    )
    return new_pas, ptype_indices, stype_indices
end

function process_pas(c::AbstractComponent, pas::ComponentVector, ptypes::AbstractVector{AbstractVector{Symbol}}, stypes::AbstractVector{AbstractVector{Symbol}})
    @assert length(c.components) == length(ptypes) == length(stypes) "ptypes and stypes length must be equal to components length"
    pas_type = eltype(pas)
    params_names, state_names = get_param_names(c), get_state_names(c)

    unq_ptypes = unique.(ptypes)
    unq_stypes = unique.(stypes)

    max_len_ptype = maximum(length.(unq_ptypes))
    max_len_stype = maximum(length.(unq_stypes))

    params_mat = zeros(pas_type, length(params_names), max_len_ptype)
    initstates_mat = zeros(pas_type, length(state_names), max_len_stype)

    for i in 1:length(c.components)
        check_pas(c.components[i], pas, ptypes[i], stypes[i])
        tmp_params_names, tmp_state_names = get_param_names(c.components[i]), get_state_names(c.components[i])
        tmp_params_indices = [findfirst(isequal(tmp_param_name), params_names) for tmp_param_name in tmp_params_names]
        tmp_state_indices = [findfirst(isequal(tmp_state_name), state_names) for tmp_state_name in tmp_state_names]

        if !isempty(tmp_params_names)
            params_mat[tmp_params_indices, 1:length(unq_ptypes[i])] .= [Vector(pas[:params][unq_ptype][tmp_params_names]) for unq_ptype in unq_ptypes[i]]
        end
        if !isempty(tmp_state_names)
            initstates_mat[tmp_params_indices, 1:length(unq_stypes[i])] .= [Vector(pas[:initstates][unq_stype][tmp_state_names]) for unq_stype in unq_stypes[i]]
        end
    end

    new_pas = ComponentVector(
        params=params_mat,
        initstates=initstates_mat,
        nns=haskey(pas, :nns) ? Vector(pas[:nns][HydroModelCore.get_nn_names(c)]) : Vector{pas_type}(undef, 0),
    )

    ptype_indices = map(unq_ptypes) do unq_ptype
        [findfirst(isequal(ptype), ptypes) for ptype in unq_ptype]
    end
    stype_indices = map(unq_stypes) do unq_stype
        [findfirst(isequal(stype), stypes) for stype in unq_stype]
    end

    return new_pas, ptype_indices, stype_indices
end