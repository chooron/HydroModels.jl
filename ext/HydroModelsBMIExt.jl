module HydroModelsBMIExt

using HydroModels
using Dates
using ComponentArrays
import BasicModelInterface as BMI

@kwdef mutable struct HydroModelsBMI{C}
    component::C

    # BMI required
    start_time::Date
    end_time::Date
    dt::Number
    time_unit::String = "days"
    time::Int = 1
    lat::Vector{Float} = []
    lon::Vector{Float} = []

    # Internal state
    state::ComponentVector # current state
    index::Int # current index

    # Independent of model
    forcing::AbstractArray
    pas::ComponentVector
    initstates::ComponentVector
    config::NamedTuple
end

"""
    initialize(::Type{Model}, [config_file::String])::Model

Perform startup tasks for the model.

Perform all tasks that take place before entering the model's time
loop, including opening files and initializing the model state. Model
inputs are read from a text-based configuration file, specified by
`config_file`.

# Arguments
- `Model`: the type of your model, only used for dispatch
- `config_file`: String, optional
    The path to the model configuration file.

# Returns
The initialized model.

# Example
If your model struct is named `MyModel`, you can implement this
in the following way:

```julia
function BMI.initialize(::Type{MyModel}, config_file)::MyModel
    ...
    m = MyModel(...)
    return m
end
```

# Notes
Models should be refactored, if necessary, to use a
configuration file. CSDMS does not impose any constraint on
how configuration files are formatted, although YAML is
recommended. A template of a model's configuration file
with placeholder values is used by the BMI.
"""
function BMI.initialize(
    component::C,
    forcing::AbstractArray,
    pas::ComponentVector,
    initstates::ComponentVector,
    timeidx::AbstractVector{<:Date},
    config::NamedTuple,
    lat::Vector{Number},
    lon::Vector{Number}
) where {C<:AbstractComponent}
    check(component=component, input=forcing, pas=pas, initstates=initstates, timeidx=timeidx)
    return HydroModelsBMI{C}(component, forcing, pas, initstates, config, initstates, timeidx[1])
end

"""
    update(model)::Nothing

Advance model state by one time step.

Perform all tasks that take place within one pass through the model's
time loop. This typically includes incrementing all of the model's
state variables. If the model's state variables don't change in time,
then they can be computed by the `initialize` method and this
method can return with no action.
"""
function BMI.update(model::HydroModelsBMI)
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


"""
    update_until(model, time::Date)::Nothing

Advance model state until the given time.

The given `time` must be a model time later than the current model time.
"""

function BMI.update_until(model::HydroModelsBMI, time::Date)
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


"""
    finalize(model)::Nothing

Perform tear-down tasks for the model.

Perform all tasks that take place after exiting the model's time
loop. This typically includes deallocating memory, closing files and
printing reports.
"""
function BMI.finalize(model::HydroModelsBMI)
    output = BMI.update_until(model, model.end_time)
    #todo clear memory
    return output
end

"""
    get_component_name(model)::String

Name of the component.
"""
function get_component_name(model::HydroModelsBMI)
    return get_name(model.component)
end

"""
    get_input_item_count(model)::Int

Count of a model's input variables.
"""
function get_input_item_count(model::HydroModelsBMI)
    return length(get_input_names(model.component))
end

"""
    get_output_item_count(model)::Int

Count of a model's output variables.
"""
function get_output_item_count(model::HydroModelsBMI)
    return length(get_output_names(model.component))
end

"""
    get_input_var_names(model)::Vector{String}

List of a model's input variables.

Input variable names must be CSDMS Standard Names, also known
as *long variable names*.

# Notes
Standard Names enable the CSDMS framework to determine whether
an input variable in one model is equivalent to, or compatible
with, an output variable in another model. This allows the
framework to automatically connect components.

Standard Names do not have to be used within the model.
"""
function get_input_var_names(model::HydroModelsBMI)
    return string.(get_input_names(model.component))
end

"""
    get_output_var_names(model)::Vector{String}

List of a model's output variables.

Output variable names must be CSDMS Standard Names, also known
as *long variable names*.
"""
function get_output_var_names(model::HydroModelsBMI)
    return string.(get_output_names(model.component))
end

"""
    get_var_grid(model, name)::Int

Get grid identifier integer for the given variable.

The `name` can be an input or output variable name, a CSDMS Standard Name.
"""
function get_var_grid(model::HydroModelsBMI, name::String)
    # 默认所有变量使用同一个网格
    return 0
end

"""
    get_var_type(model, name)::String

Get data type of the given variable.

The `name` can be an input or output variable name, a CSDMS Standard Name.
"""
function get_var_type(model::HydroModelsBMI, name::String)
    # 假设所有变量都是Float64类型
    return eltype(model.forcing)
end

"""
    get_var_units(model, name)::String

Get units of the given variable.

Standard unit names, in lower case, should be used, such as
``meters`` or ``seconds``. Standard abbreviations, like ``m`` for
meters, are also supported. For variables with compound units,
each unit name is separated by a single space, with exponents
other than 1 placed immediately after the name, as in ``m s-1``
for velocity, ``W m-2`` for an energy flux, or ``km2`` for an
area.

The `name` can be an input or output variable name, a CSDMS Standard Name.

# Notes
CSDMS uses the [UDUNITS](http://www.unidata.ucar.edu/software/udunits)
standard from Unidata.
"""
function get_var_units(model::HydroModelsBMI, name::String)
    # 这里需要根据具体变量返回对应的单位
    # 可以从组件的元数据中获取
    return "mm/day"  # 默认单位，实际应用中应该从组件元数据获取
end

"""
    get_var_itemsize(model, name)::Int

Get memory use for each array element in bytes.

The `name` can be an input or output variable name, a CSDMS Standard Name.
"""
function get_var_itemsize(model::HydroModelsBMI, name::String)
    # 假设所有变量都是Float64类型
    return sizeof(eltype(model.forcing))
end

"""
    get_var_nbytes(model, name)::Int

Get size, in bytes, of the given variable.

The `name` can be an input or output variable name, a CSDMS Standard Name.
"""
function get_var_nbytes(model::HydroModelsBMI, name::String)
    # 计算变量总字节数 = 元素大小 × 元素数量
    # 这里简化处理，假设每个变量都是一个时间序列数组
    return get_var_itemsize(model, name) * model.component.config.n_timesteps
end

"""
    get_var_location(model, name)::String

Get the grid element type that the a given variable is defined on.

The grid topology can be composed of *nodes*, *edges*, and *faces*.

- node: A point that has a coordinate pair or triplet: the most
    basic element of the topology.
- edge: A line or curve bounded by two nodes.
- face: A plane or surface enclosed by a set of edges. In a 2D
    horizontal application one may consider the word "polygon",
    but in the hierarchy of elements the word "face" is most common.

The `name` can be an input or output variable name, a CSDMS Standard Name.

# Returns
The grid location on which the variable is defined. Must be one of
`"node"`, `"edge"`, or `"face"`.

# Notes
CSDMS uses the [ugrid conventions](http://ugrid-conventions.github.io/ugrid-conventions)
to define unstructured grids.
"""
function get_var_location(model::HydroModelsBMI, name::String)
    # 大多数水文模型变量位于节点上
    return "node"
end

"""
    get_current_time(model)::Float64

Current time of the model.
"""
function get_current_time(model::HydroModelsBMI)
    return model.time
end

"""
    get_start_time(model)::Float64

Start time of the model.
"""
function get_start_time(model::HydroModelsBMI)
    return model.start_time
end

"""
    get_end_time(model)::Float64

End time of the model.
"""
function get_end_time(model::HydroModelsBMI)
    return model.end_time
end

"""
    get_time_units(model)::String

Time units of the model; e.g., `days` or `s`.

# Notes
CSDMS uses the [UDUNITS](http://www.unidata.ucar.edu/software/udunits)
standard from Unidata.
"""
function get_time_units(model::HydroModelsBMI)
    return model.time_unit
end

"""
    get_time_step(model)::Float64

Current time step of the model.
"""
"""    
    get_time_step(model)::Float64

Current time step of the model.
"""
function get_time_step(model::HydroModelsBMI)
    return model.dt
end

"""
    get_value(model, name, dest)::DenseVector

Get a copy of values of the given variable.

This is a getter for the model, used to access the model's
current state. It returns a *copy* of a model variable, with
the return type, size and rank dependent on the variable.

# Arguments
- `name`: An input or output variable name, a CSDMS Standard Name.
- `dest`: An array into which to place the values.

# Returns
The same array that was passed as an input buffer, `dest`.
"""
function get_value(model::HydroModelsBMI, name::String, dest::AbstractArray)
    # 从模型中获取变量值
    # 这里需要根据变量名从模型状态中获取对应的值
    # 简化实现，假设所有变量都可以从组件的状态中获取
    var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
    if !isnothing(var_index)
        # 如果是输出变量
        copyto!(dest, model.component.states[var_index, :])
    else
        # 如果是输入变量
        var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
        if !isnothing(var_index)
            copyto!(dest, model.component.inputs[var_index, :])
        else
            error("Variable $name not found in model")
        end
    end
    return dest
end

"""
    get_value_ptr(model, name)::DenseVector

Get a reference to values of the given variable.

This is a getter for the model, used to access the model's
current state. It returns a reference to a model variable,
with the return type, size and rank dependent on the variable.

The `name` can be an input or output variable name, a CSDMS Standard Name.
"""
function get_value_ptr(model::HydroModelsBMI, name::String)
    # Julia不支持直接返回指针，但可以返回对数组的引用
    # 这里简化实现，返回对应变量的视图
    var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
    if !isnothing(var_index)
        # 如果是输出变量
        return @view model.component.states[var_index, :]
    else
        # 如果是输入变量
        var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
        if !isnothing(var_index)
            return @view model.component.inputs[var_index, :]
        else
            error("Variable $name not found in model")
        end
    end
end

"""
    get_value_at_indices(model, name, dest, inds)::DenseVector

Get values at particular indices.

# Arguments
- `name`: An input or output variable name, a CSDMS Standard Name.
- `dest`: An array into which to place the values.
- `inds`: The indices into the variable array.

# Returns
The same array that was passed as an input buffer, `dest`.
"""
function get_value_at_indices(model::HydroModelsBMI, name::String, dest::AbstractArray, indices::AbstractArray{Int})
    # 获取指定索引处的变量值
    values = get_value_ptr(model, name)
    for (i, idx) in enumerate(indices)
        dest[i] = values[idx]
    end
    return dest
end

"""
    set_value(model, name::String, value)::Nothing

Specify a new value for a model variable.

This is the setter for the model, used to change the model's
current state. It accepts, through *value*, a new value for a
model variable, with the type, size and rank of *value*
dependent on the variable.

# Arguments
- `name`: An input or output variable name, a CSDMS Standard Name.
- `value`: The new value for the specified variable.
"""
function set_value(model::HydroModelsBMI, name::String, src::AbstractArray)
    # 设置变量值
    # 这里需要根据变量名设置模型中对应的值
    var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
    if !isnothing(var_index)
        # 如果是输入变量
        copyto!(model.component.inputs[var_index, :], src)
    else
        # 如果是状态变量
        var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
        if !isnothing(var_index)
            copyto!(model.component.states[var_index, :], src)
        else
            error("Variable $name not found in model or is not settable")
        end
    end
    return nothing
end

"""
    set_value_at_indices(model, name::String, inds::DenseVector{Int}, value)::Nothing

Specify a new value for a model variable at particular indices.

# Arguments
- `name`: An input or output variable name, a CSDMS Standard Name.
- `inds`: The indices into the variable array.
- `value`: The new value for the specified variable.
"""
function set_value_at_indices(model::HydroModelsBMI, name::String, indices::AbstractArray{Int}, src::AbstractArray)
    # 在指定索引处设置变量值
    var_index = findfirst(x -> string(x) == name, get_input_names(model.component))
    if !isnothing(var_index)
        # 如果是输入变量
        for (i, idx) in enumerate(indices)
            model.component.inputs[var_index, idx] = src[i]
        end
    else
        # 如果是状态变量
        var_index = findfirst(x -> string(x) == name, get_output_names(model.component))
        if !isnothing(var_index)
            for (i, idx) in enumerate(indices)
                model.component.states[var_index, idx] = src[i]
            end
        else
            error("Variable $name not found in model or is not settable")
        end
    end
    return nothing
end

# Grid information

"""
    get_grid_rank(model, grid::Int)::Int

Get number of dimensions of the computational grid.
"""
function get_grid_rank(model::HydroModelsBMI, grid::Int)
    # 水文模型通常使用1D网格（时间序列）
    return 1
end

"""
    get_grid_size(model, grid::Int)::Int

Get the total number of elements in the computational grid.
"""
function get_grid_size(model::HydroModelsBMI, grid::Int)
    # 返回网格中的元素数量，对于时间序列，就是时间步数
    return model.component.config.n_timesteps
end

"""
    get_grid_type(model, grid::Int)::String

Get the grid type as a string.
"""
function get_grid_type(model::HydroModelsBMI, grid::Int)
    # 水文模型通常使用一维均匀网格
    return "uniform_rectilinear"
end

# Uniform rectilinear

"""
    get_grid_shape(model, grid::Int, shape::DenseVector{Int})::DenseVector{Int}

Get dimensions of the computational grid.

Returns the filled `shape` array.
"""
function get_grid_shape(model::HydroModelsBMI, grid::Int, shape::AbstractVector{Int})
    # 对于一维网格，形状就是时间步数
    shape[1] = model.component.config.n_timesteps
    return shape
end

"""
    get_grid_spacing(model, grid::Int, spacing::DenseVector{Float64})::DenseVector{Float64}

Get distance between nodes of the computational grid.

# Arguments
- `grid`: A grid identifier.
- `spacing`: An array to hold the spacing between grid rows and columns.

Returns the filled `spacing` array.
"""
function get_grid_spacing(model::HydroModelsBMI, grid::Int, spacing::AbstractVector{Float64})
    # 对于时间序列，间距就是时间步长
    spacing[1] = model.component.config.step_size
    return spacing
end

"""
    get_grid_origin(model, grid::Int, origin::DenseVector{Float64})::DenseVector{Float64}

Get coordinates for the lower-left corner of the computational grid.

# Arguments
- `grid`: A grid identifier.
- `origin`: An array to hold the coordinates of the lower-left corner of the grid.

Returns the filled `origin` array.
"""
function get_grid_origin(model::HydroModelsBMI, grid::Int, origin::AbstractVector{Float64})
    # 对于时间序列，原点就是起始时间
    origin[1] = model.component.config.start_time
    return origin
end

# Non-uniform rectilinear, curvilinear

"""
    get_grid_x(model, grid::Int, x::DenseVector{Float64})::DenseVector{Float64}

Get coordinates of grid nodes in the x direction.

# Arguments
- `grid`: A grid identifier.
- `x`: An array to hold the x-coordinates of the grid nodes.

Returns the filled `x` array.
"""
function get_grid_x(model::HydroModelsBMI, grid::Int, x::AbstractVector{Float64})
    # 对于时间序列，x坐标就是时间点
    for i in 1:length(x)
        x[i] = model.component.config.start_time + (i - 1) * model.component.config.step_size
    end
    return x
end

"""
    get_grid_y(model, grid::Int, y::DenseVector{Float64})::DenseVector{Float64}

Get coordinates of grid nodes in the y direction.

# Arguments
- `grid`: A grid identifier.
- `y`: An array to hold the y-coordinates of the grid nodes.

Returns the filled `y` array.
"""
function get_grid_y(model::HydroModelsBMI, grid::Int, y::AbstractVector{Float64})
    # 对于一维时间序列，y坐标不适用
    # 返回全零数组
    fill!(y, 0.0)
    return y
end

"""
    get_grid_z(model, grid::Int, z::DenseVector{Float64})::DenseVector{Float64}

Get coordinates of grid nodes in the z direction.

# Arguments
- `grid`: A grid identifier.
- `z`: An array to hold the z-coordinates of the grid nodes.

Returns the filled `z` array.
"""
function get_grid_z(model::HydroModelsBMI, grid::Int, z::AbstractVector{Float64})
    # 对于一维时间序列，z坐标不适用
    # 返回全零数组
    fill!(z, 0.0)
    return z
end

"""
    get_grid_node_count(model, grid::Int)::Int

Get the number of nodes in the grid.
"""
function get_grid_node_count(model::HydroModelsBMI, grid::Int)
    # 对于一维时间序列，节点数就是时间步数
    return model.component.config.n_timesteps
end

"""
    get_grid_edge_count(model, grid::Int)::Int

Get the number of edges in the grid.
"""
function get_grid_edge_count(model::HydroModelsBMI, grid::Int)
    # 对于一维时间序列，边数是节点数减1
    return model.component.config.n_timesteps - 1
end

"""
    get_grid_face_count(model, grid::Int)::Int

Get the number of faces in the grid.
"""
function get_grid_face_count(model::HydroModelsBMI, grid::Int)
    # 对于一维时间序列，没有面
    return 0
end

"""
    get_grid_edge_nodes(model, grid::Int, edge_nodes::DenseVector{Int})::DenseVector{Int}

Get the edge-node connectivity.

# Arguments
- `grid`: A grid identifier.
- `edge_nodes`: An array of integers to place the edge-node connectivity.
    For each edge, connectivity is given as node at edge tail,
    followed by node at edge head. Therefore this array must be twice
    the number of nodes long.

Returns the filled `edge_nodes` array.
"""
function get_grid_edge_nodes(model::HydroModelsBMI, grid::Int, edge_nodes::AbstractVector{Int})
    # 对于一维时间序列，每条边连接相邻的两个时间点
    n_edges = get_grid_edge_count(model, grid)
    for i in 1:n_edges
        edge_nodes[2*i-1] = i        # 边的起点
        edge_nodes[2*i] = i + 1      # 边的终点
    end
    return edge_nodes
end

"""
    get_grid_face_edges(model, grid::Int, face_edges::DenseVector{Int})::DenseVector{Int}

Get the face-edge connectivity.

# Arguments
- `grid`: A grid identifier.
- `face_edges`: An array of integers in which to place the face-edge connectivity.

Returns the filled `face_edges` array.
"""
function get_grid_face_edges(model::HydroModelsBMI, grid::Int, face_edges::AbstractVector{Int})
    # 对于一维时间序列，没有面，返回空数组
    return face_edges
end

"""
    get_grid_face_nodes(model, grid::Int, face_nodes::DenseVector{Int})::DenseVector{Int}

Get the face-node connectivity.

# Arguments
- `grid`: A grid identifier.
- `face_nodes`: An array of integers in which to place the face-node connectivity.
    For each face, the nodes (listed in a counter-clockwise direction) that form the
    boundary of the face.

Returns the filled `face_nodes` array.
"""
function get_grid_face_nodes(model::HydroModelsBMI, grid::Int, face_nodes::AbstractVector{Int})
    # 对于一维时间序列，没有面，返回空数组
    return face_nodes
end

"""
    get_grid_nodes_per_face(model, grid::Int, nodes_per_face::DenseVector{Int})::DenseVector{Int}

Get the number of nodes for each face.

# Arguments
- `grid`: A grid identifier.
- `nodes_per_face`: An array in which to place the number of edges per face.

Returns the filled `nodes_per_face` array.
"""
function get_grid_nodes_per_face(model::HydroModelsBMI, grid::Int, nodes_per_face::AbstractVector{Int})
    # 对于一维时间序列，没有面，返回空数组
    return nodes_per_face
end
end