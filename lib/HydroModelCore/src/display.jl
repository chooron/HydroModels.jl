# ================================================================================================
# 显示配置和辅助函数
# ================================================================================================

"""
Color scheme for component display.
"""
const DISPLAY_COLORS = (
    header = :light_blue,
    name = :light_black,
    label_input = :light_green,
    label_output = :light_yellow,
    label_state = :blue,
    label_param = :light_magenta,
    label_nn = :light_cyan,
    label_expr = :cyan,
    value = :white
)

"""
$(SIGNATURES)

Print a styled header for component display.
"""
@inline function print_component_header(io::IO, type_name::String, component_name::String="")
    printstyled(io, "┌ ", color=DISPLAY_COLORS.header, bold=true)
    printstyled(io, type_name, color=DISPLAY_COLORS.header, bold=true)
    if !isempty(component_name)
        printstyled(io, "{$component_name}", color=DISPLAY_COLORS.name)
    end
    println(io)
end

"""
$(SIGNATURES)

Print a labeled field in component display.
"""
@inline function print_field(io::IO, label::String, values::AbstractVector{Symbol}, color::Symbol=:white; prefix::String="│ ")
    print(io, prefix)
    printstyled(io, label, color=color)
    println(io, "[", join(values, ", "), "]")
end

"""
$(SIGNATURES)

Print expressions with labels.
"""
@inline function print_expressions(io::IO, labels::AbstractVector{Symbol}, exprs, label_color::Symbol=:yellow)
    print(io, "│ ")
    printstyled(io, "Expressions:", color=DISPLAY_COLORS.label_expr)
    println(io)
    for (label, expr) in zip(labels, exprs)
        print(io, "│   ")
        printstyled(io, "$label = ", color=label_color)
        println(io, expr)
    end
end

"""
$(SIGNATURES)

Print a footer for component display.
"""
@inline function print_component_footer(io::IO)
    printstyled(io, "└─", color=DISPLAY_COLORS.header)
    println(io)
end

"""
$(SIGNATURES)

Print compact representation of a component.
"""
@inline function print_compact(io::IO, type_name::String, fields::NamedTuple)
    print(io, type_name, "(")
    field_strs = String[]
    for (key, val) in pairs(fields)
        if val isa AbstractVector
            push!(field_strs, "$key=[$(join(val, ","))]")
        else
            push!(field_strs, "$key=$val")
        end
    end
    print(io, join(field_strs, ", "))
    print(io, ")")
end


# ================================================================================================
# AbstractHydroFlux 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractHydroFlux`.
"""
function Base.show(io::IO, flux::AbstractHydroFlux)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "HydroFlux", (
            inputs=get_input_names(flux),
            outputs=get_output_names(flux),
            params=get_param_names(flux)
        ))
    else
        print_component_header(io, "HydroFlux", flux.name)
        print_field(io, "Inputs:  ", get_input_names(flux), DISPLAY_COLORS.label_input)
        print_field(io, "Outputs: ", get_output_names(flux), DISPLAY_COLORS.label_output)
        print_field(io, "Params:  ", get_param_names(flux), DISPLAY_COLORS.label_param)

        if !isempty(flux.exprs)
            print_expressions(io, get_output_names(flux), flux.exprs)
        end

        print_component_footer(io)
    end
end


# ================================================================================================
# AbstractStateFlux 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractStateFlux`.
"""
function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "StateFlux", (
            inputs=get_input_names(flux),
            states=get_state_names(flux),
            params=get_param_names(flux)
        ))
    else
        print_component_header(io, "StateFlux", flux.name)
        print_field(io, "Inputs:  ", get_input_names(flux), DISPLAY_COLORS.label_input)
        print_field(io, "States:  ", get_state_names(flux), DISPLAY_COLORS.label_state)
        print_field(io, "Params:  ", get_param_names(flux), DISPLAY_COLORS.label_param)

        if !isempty(get_exprs(flux))
            print_expressions(io, get_state_names(flux), get_exprs(flux), DISPLAY_COLORS.label_state)
        end

        print_component_footer(io)
    end
end


# ================================================================================================
# AbstractNeuralFlux 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractNeuralFlux`.
"""
function Base.show(io::IO, flux::AbstractNeuralFlux)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "NeuralFlux", (
            inputs=get_input_names(flux),
            outputs=get_output_names(flux),
            nns=get_nn_names(flux)
        ))
    else
        print_component_header(io, "NeuralFlux", flux.name)
        print_field(io, "Inputs:  ", get_input_names(flux), DISPLAY_COLORS.label_input)
        print_field(io, "Outputs: ", get_output_names(flux), DISPLAY_COLORS.label_output)
        print_field(io, "NNs:     ", get_nn_names(flux), DISPLAY_COLORS.label_nn)

        print(io, "│ ")
        printstyled(io, "Expressions:", color=DISPLAY_COLORS.label_expr)
        println(io)
        print(io, "│   ")
        printstyled(io, "$(get_output_names(flux)) = ", color=DISPLAY_COLORS.label_output)
        println(io, "$(flux.chain)($(get_input_names(flux)))")

        print_component_footer(io)
    end
end


# ================================================================================================
# AbstractHydrograph 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractHydrograph`.
"""
function Base.show(io::IO, uh::AbstractHydrograph)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "UnitHydroFlux", (
            inputs=get_input_names(uh),
            outputs=get_output_names(uh),
            params=get_param_names(uh)
        ))
    else
        print_component_header(io, "UnitHydroFlux")
        print_field(io, "Inputs:       ", get_input_names(uh), DISPLAY_COLORS.label_input)
        print_field(io, "Outputs:      ", get_output_names(uh), DISPLAY_COLORS.label_output)
        print_field(io, "Parameters:   ", get_param_names(uh), DISPLAY_COLORS.label_param)
        print_component_footer(io)
    end
end


# ================================================================================================
# AbstractBucket 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractBucket`.
"""
function Base.show(io::IO, ele::AbstractBucket)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "HydroBucket", (
            inputs=get_input_names(ele),
            states=get_state_names(ele),
            outputs=get_output_names(ele),
            params=get_param_names(ele),
            nns=get_nn_names(ele)
        ))
    else
        print_component_header(io, "HydroBucket", ele.name)
        print_field(io, "Inputs:  ", get_input_names(ele), DISPLAY_COLORS.label_input)
        print_field(io, "States:  ", get_state_names(ele), DISPLAY_COLORS.label_state)
        print_field(io, "Outputs: ", get_output_names(ele), DISPLAY_COLORS.label_output)
        print_field(io, "Params:  ", get_param_names(ele), DISPLAY_COLORS.label_param)
        print_field(io, "NNs:     ", get_nn_names(ele), DISPLAY_COLORS.label_nn, prefix="└─")
    end
end


# ================================================================================================
# AbstractHydroRoute 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractHydroRoute`.
"""
function Base.show(io::IO, route::AbstractHydroRoute)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "HydroRoute", (
            inputs=get_input_names(route),
            states=get_state_names(route),
            outputs=get_output_names(route),
            params=get_param_names(route),
            nns=get_nn_names(route)
        ))
    else
        print_component_header(io, "HydroRoute", route.name)
        print_field(io, "Inputs:  ", get_input_names(route), DISPLAY_COLORS.label_input)
        print_field(io, "States:  ", get_state_names(route), DISPLAY_COLORS.label_state)
        print_field(io, "Outputs: ", get_output_names(route), DISPLAY_COLORS.label_output)
        print_field(io, "Params:  ", get_param_names(route), DISPLAY_COLORS.label_param)
        print_field(io, "NNs:     ", get_nn_names(route), DISPLAY_COLORS.label_nn, prefix="└─")
    end
end


# ================================================================================================
# AbstractModel 显示
# ================================================================================================

"""
$(SIGNATURES)

Custom pretty-printing for `AbstractModel`.
"""
function Base.show(io::IO, model::AbstractModel)
    compact = get(io, :compact, false)

    if compact
        print_compact(io, "HydroModel", (
            name=model.name,
            components=length(model.components)
        ))
    else
        # 统计各类组件
        fluxes_in_model = filter(x -> x isa AbstractFlux, model.components)
        buckets_in_model = filter(x -> x isa AbstractBucket, model.components)
        routes_in_model = filter(x -> x isa AbstractRoute, model.components)
        uh_in_model = filter(x -> x isa AbstractHydrograph, model.components)

        # 打印头部
        printstyled(io, "┌ ", color=DISPLAY_COLORS.header, bold=true)
        printstyled(io, "HydroModel: ", color=DISPLAY_COLORS.header, bold=true)
        printstyled(io, model.name, color=:white, bold=true)
        println(io)

        # 打印组件列表
        print(io, "│ ")
        printstyled(io, "Components: ", color=DISPLAY_COLORS.label_expr)
        println(io, join(map(c -> c.name, model.components), ", "))

        # 打印模型变量信息
        print_field(io, "Inputs:  ", get_input_names(model), DISPLAY_COLORS.label_input)
        print_field(io, "States:  ", get_state_names(model), DISPLAY_COLORS.label_state)
        print_field(io, "Outputs: ", get_output_names(model), DISPLAY_COLORS.label_output)
        print_field(io, "Params:  ", get_param_names(model), DISPLAY_COLORS.label_param)
        print_field(io, "NNs:     ", get_nn_names(model), DISPLAY_COLORS.label_nn)

        # 打印组件统计
        print(io, "│ ")
        printstyled(io, "Summary:", color=:white, bold=true)
        println(io)

        print(io, "│   ")
        printstyled(io, "Fluxes:          ", color=DISPLAY_COLORS.label_expr)
        println(io, length(fluxes_in_model), " flux", length(fluxes_in_model) == 1 ? "" : "es")

        print(io, "│   ")
        printstyled(io, "Buckets:         ", color=DISPLAY_COLORS.label_expr)
        println(io, length(buckets_in_model), " bucket", length(buckets_in_model) == 1 ? "" : "s")

        print(io, "│   ")
        printstyled(io, "Routes:          ", color=DISPLAY_COLORS.label_expr)
        println(io, length(routes_in_model), " route", length(routes_in_model) == 1 ? "" : "s")

        print(io, "│   ")
        printstyled(io, "UnitHydrographs: ", color=DISPLAY_COLORS.label_expr)
        println(io, length(uh_in_model), " uh", length(uh_in_model) == 1 ? "" : "s")

        print_component_footer(io)
    end
end


# ================================================================================================
# 额外的显示工具函数
# ================================================================================================

"""
$(SIGNATURES)

Print a summary table of multiple components.

# Arguments
- `components::AbstractVector{<:AbstractComponent}`: Components to summarize
- `io::IO=stdout`: Output stream
"""
function print_component_table(components::AbstractVector{<:AbstractComponent}; io::IO=stdout)
    println(io, "Component Summary Table")
    println(io, "=" ^ 80)

    # 打印表头
    @printf(io, "%-20s %-15s %-10s %-10s %-10s %-10s %-10s\n",
        "Name", "Type", "#Inputs", "#Outputs", "#States", "#Params", "#NNs")
    println(io, "-" ^ 80)

    # 打印每个组件
    for comp in components
        comp_type = split(string(typeof(comp)), ".")[end]
        @printf(io, "%-20s %-15s %-10d %-10d %-10d %-10d %-10d\n",
            get_name(comp),
            comp_type,
            length(get_input_names(comp)),
            length(get_output_names(comp)),
            length(get_state_names(comp)),
            length(get_param_names(comp)),
            length(get_nn_names(comp))
        )
    end

    println(io, "=" ^ 80)
end

"""
$(SIGNATURES)

Print detailed information about component variables.

# Arguments
- `component::AbstractComponent`: Component to inspect
- `io::IO=stdout`: Output stream
"""
function print_variable_details(component::AbstractComponent; io::IO=stdout)
    println(io, "Variable Details for Component: $(get_name(component))")
    println(io, "=" ^ 80)

    # 输入变量
    inputs = get_input_names(component)
    if !isempty(inputs)
        printstyled(io, "Input Variables ($(length(inputs))):\n", color=DISPLAY_COLORS.label_input, bold=true)
        for (i, var) in enumerate(inputs)
            println(io, "  $i. $var")
        end
        println(io)
    end

    # 输出变量
    outputs = get_output_names(component)
    if !isempty(outputs)
        printstyled(io, "Output Variables ($(length(outputs))):\n", color=DISPLAY_COLORS.label_output, bold=true)
        for (i, var) in enumerate(outputs)
            println(io, "  $i. $var")
        end
        println(io)
    end

    # 状态变量
    states = get_state_names(component)
    if !isempty(states)
        printstyled(io, "State Variables ($(length(states))):\n", color=DISPLAY_COLORS.label_state, bold=true)
        for (i, var) in enumerate(states)
            println(io, "  $i. $var")
        end
        println(io)
    end

    # 参数
    params = get_param_names(component)
    if !isempty(params)
        printstyled(io, "Parameters ($(length(params))):\n", color=DISPLAY_COLORS.label_param, bold=true)
        for (i, var) in enumerate(params)
            println(io, "  $i. $var")
        end
        println(io)
    end

    # 神经网络
    nns = get_nn_names(component)
    if !isempty(nns)
        printstyled(io, "Neural Networks ($(length(nns))):\n", color=DISPLAY_COLORS.label_nn, bold=true)
        for (i, var) in enumerate(nns)
            println(io, "  $i. $var")
        end
        println(io)
    end

    println(io, "=" ^ 80)
end
