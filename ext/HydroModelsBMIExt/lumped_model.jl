@kwdef mutable struct LumpedModelBMI{C}
    # Component
    component::C

    # BMI required
    start_time::Date
    end_time::Date
    dt::Number
    time_unit::String = "days"
    time::Int = 1

    # Internal state
    state::ComponentVector # current state
    index::Int # current index

    # Independent of model
    forcing::AbstractArray
    pas::ComponentVector
    initstates::ComponentVector
    config::NamedTuple
end

function BMI.initialize(
    component::C, forcing::AbstractArray, pas::ComponentVector,
    initstates::ComponentVector, timeidx::AbstractVector{<:Date}, config::NamedTuple
) where {C<:AbstractComponent}
    check(component=component, input=forcing, pas=pas, initstates=initstates, timeidx=timeidx)
    return LumpedModelBMI{C}(component, forcing, pas, initstates, config, initstates, timeidx[1])
end

function BMI.update(model::LumpedModelBMI)
    output = model.component(
        model.forcing, model.pas,
        initstates=ComponentVector(model.state),
        timeidx=model.timeidx[model.index:model.index+1],
        config=model.config
    )
    model_state_names = get_state_names(model.component)
    model.state = NamedTuple{model_state_names}(eachslice(output, dims=1)[1:length(model_state_names)])
    model.index += 1
end

function BMI.update_until(model::LumpedModelBMI, time::Date)
    end_time = ceil(Int, (time - model.start_time) / model.dt)
    @assert end_time >= model.time
    output = model.component(
        model.forcing, model.pas,
        initstates=ComponentVector(model.state),
        timeidx=model.timeidx[model.index:end_time],
        config=model.config
    )
    model_state_names = get_state_names(model.component)
    model_output_names = get_output_names(model.component)
    latest_output = output[ntuple((_) -> Colon(), length(model_state_names))..., end:end]
    model.state = NamedTuple{Tuple(model_state_names)}(eachslice(latest_output, dims=1))
    model.index = end_time
    return NamedTuple{Tuple(vcat(model_output_names, model_state_names))}(eachslice(output, dims=1))
end

function BMI.finalize(model::LumpedModelBMI)
    output = BMI.update_until(model, model.end_time)
    # todo clear memory
    return output
end

get_component_name(model::LumpedModelBMI) = get_name(model.component)

get_input_item_count(model::LumpedModelBMI) = length(get_input_names(model.component))

get_output_item_count(model::LumpedModelBMI) = length(get_output_names(model.component))

get_input_var_names(model::LumpedModelBMI) = string.(get_input_names(model.component))

get_output_var_names(model::LumpedModelBMI) = string.(vcat(get_state_names(model.component), get_output_names(model.component)))

get_var_grid(::LumpedModelBMI, ::String) = 0

get_var_type(model::LumpedModelBMI, ::String) = repr(eltype(model.forcing))

function get_var_units(model::LumpedModelBMI, name::String)
    if name in get_output_names(model.component)
        return get_unit(model.component.infos.outputs[name])
    elseif name in get_input_names(model.component)
        return get_unit(model.component.infos.inputs[name])
    elseif name in get_state_names(model.component)
        return get_unit(model.component.infos.states[name])
    else
        error("Variable $name not found in model")
    end
end

get_var_itemsize(model::LumpedModelBMI, ::String) = sizeof(eltype(model.forcing))

get_var_nbytes(model::LumpedModelBMI, name::String) = get_var_itemsize(model, name) * 10000 # todo add get total timesteps

get_var_location(model::LumpedModelBMI, name::String) = "node"

get_current_time(model::LumpedModelBMI) = model.time

get_start_time(model::LumpedModelBMI) = model.start_time

get_end_time(model::LumpedModelBMI) = model.end_time

get_time_units(model::LumpedModelBMI) = model.time_unit

get_time_step(model::LumpedModelBMI) = model.dt

function get_value(model::LumpedModelBMI, name::String, dest::AbstractArray)
    dest .= copy(get_value_ptr(model, name))
    return dest
end

function get_value_ptr(model::LumpedModelBMI, name::String)
    if name in get_output_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
        return model.component.outputs[var_index, ntuple((_) -> Colon(), size(model.component.outputs) - 1)]
    elseif name in get_input_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
        return model.component.inputs[var_index, ntuple((_) -> Colon(), size(model.component.inputs) - 1)]
    elseif name in get_state_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_state_names(model.component))
        return model.component.states[var_index, ntuple((_) -> Colon(), size(model.component.states) - 1)]
    else
        error("Variable $name not found in model")
    end
end

function get_value_at_indices(model::LumpedModelBMI, name::String, dest::AbstractArray, indices::AbstractArray{Int})
    # 获取指定索引处的变量值
    values = get_value_ptr(model, name)
    dest .= values[ntuple((_) -> Colon(), size(values) - 1), indices]
    return dest
end

function set_value(model::LumpedModelBMI, name::String, src::AbstractArray)
    if name in get_output_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
        model.component.outputs[var_index, ntuple((_) -> Colon(), size(model.component.outputs) - 1)] .= src
    elseif name in get_input_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
        model.component.inputs[var_index, ntuple((_) -> Colon(), size(model.component.inputs) - 1)] .= src
    elseif name in get_state_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_state_names(model.component))
        model.component.states[var_index, ntuple((_) -> Colon(), size(model.component.states) - 1)] .= src
    else
        error("Variable $name not found in model")
    end
end

function set_value_at_indices(model::LumpedModelBMI, name::String, indices::AbstractArray{Int}, src::AbstractArray)
    if name in get_output_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
        model.component.outputs[var_index, ntuple((_) -> Colon(), size(model.component.outputs) - 2), indices] .= src
    elseif name in get_input_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
        model.component.inputs[var_index, ntuple((_) -> Colon(), size(model.component.inputs) - 2), indices] .= src
    elseif name in get_state_names(model.component)
        var_index = findfirst(x -> string(x) == name, get_state_names(model.component))
        model.component.states[var_index, ntuple((_) -> Colon(), size(model.component.states) - 2), indices] .= src
    else
        error("Variable $name not found in model")
    end
end

# Grid information
get_grid_rank(::LumpedModelBMI, ::Int) = 2 # fixed in LumpedModelBMI

get_grid_size(::LumpedModelBMI, ::Int) = 1 # fixed in LumpedModelBMI

get_grid_type(::LumpedModelBMI, ::Int) = "uniform_rectilinear"

get_grid_shape(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = (1, 1) # fixed in LumpedModelBMI

get_grid_spacing(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = (0, 0)

get_grid_origin(model::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = (model.lat, model.lon)

get_grid_x(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = 1

get_grid_y(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = 1

get_grid_z(::LumpedModelBMI, ::Int, ::AbstractVector{Float64}) = 1

get_grid_node_count(model::LumpedModelBMI, grid::Int) = 1 # fixed in LumpedModelBMI

get_grid_edge_count(model::LumpedModelBMI, grid::Int) = 1 # fixed in LumpedModelBMI

get_grid_face_count(::LumpedModelBMI, ::Int) = 1 # fixed in LumpedModelBMI

function get_grid_edge_nodes(model::LumpedModelBMI, grid::Int, edge_nodes::AbstractVector{Int})
    # 对于一维时间序列，每条边连接相邻的两个时间点
    n_edges = get_grid_edge_count(model, grid)
    for i in 1:n_edges
        edge_nodes[2*i-1] = i        # 边的起点
        edge_nodes[2*i] = i + 1      # 边的终点
    end
    return edge_nodes
end

get_grid_face_edges(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = 1

get_grid_face_nodes(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = 1

get_grid_nodes_per_face(::LumpedModelBMI, ::Int, ::AbstractVector{Int}) = 1