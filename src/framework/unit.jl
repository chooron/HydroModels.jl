@kwdef mutable struct Unit{T,E} <: AbstractUnit where {T<:Number,E<:AbstractElement}
    id::String

    # parameters
    parameters::ComponentVector{T}

    # model structure
    elements::Vector{E} #::Vector{ODEElement}

    # attribute
    input_names::Set{Symbol}
    output_names::Set{Symbol}

    # inner variables
    fluxes::ComponentVector = ComponentVector()
end

function build_unit(; id::String, elements::Vector{E}) where {E<:AbstractElement}
    parameters = ComponentVector()
    input_names = Set{Symbol}()
    output_names = Set{Symbol}()
    for ele in elements
        parameters = ComponentVector(parameters; ele.parameters...)
        union!(input_names, ele.input_names)
        union!(output_names, ele.output_names)
    end
    Unit(
        id=id,
        parameters=parameters,
        elements=elements,
        input_names=input_names,
        output_names=output_names
    )
end

function update_fluxes!(unit::AbstractUnit; fluxes::ComponentVector)
    unit.fluxes = ComponentVector(unit.fluxes; fluxes...)
end

function get_fluxes(unit::AbstractUnit; flux_names::Set{Symbol})
    output = Dict{Symbol,Vector}()
    for flux_nm in flux_names
        output[flux_nm] = unit.fluxes[flux_nm]
    end
    return ComponentVector(; output...)
end

function get_init_states(sort_eles::Vector{E}) where {E<:AbstractElement}
    u_init = ComponentVector()
    for tmp_ele in sort_eles
        if isa(tmp_ele, ODEElement)
            u_init = ComponentVector(u_init; tmp_ele.init_states...)
        end
    end
    return u_init
end

function set_parameters!(unit::Unit; paraminfos::Vector{ParamInfo{T}}) where {T<:Number}
    for idx in topological_sort(unit.structure)
        tmp_ele = get_prop(unit.structure, idx, :ele)
        if isa(tmp_ele, ParameterizedElement) | isa(tmp_ele, StateParameterizedElement)
            set_parameters!(tmp_ele, paraminfos=paraminfos)
        end
        if isa(tmp_ele, StateElement) | isa(tmp_ele, StateParameterizedElement)
            set_states!(tmp_ele, paraminfos=paraminfos)
        end
    end
end

function get_output(unit::Unit; input::ComponentVector{T}, step::Bool=true, kwargs...) where {T<:Number}
    # 开始计算
    if step
        # * This function is calculated element by element
        # initialize unit fluxes
        unit.fluxes = input
        # traversal of the directed graph
        for tmp_ele in unit.elements
            tmp_fluxes = get_output(tmp_ele, input=unit.fluxes)
            update_fluxes!(unit, fluxes=tmp_fluxes)
        end
    else
        # * This function is calculated based on the whole Unit
        dt = 1
        xs = 1:dt:length(input[first(keys(input))])
        # fit interpolation functions
        itp = Dict(k => linear_interpolation(xs, input[k]) for k in keys(input))
        # 获取ODEsElement的所有state初始值
        u_init = get_init_states(unit.elements)
        # 定义整体的ode函数
        function ode_function!(du, u, p, t)
            # element input
            # 使用插值方法获取该时段下的输入值
            tmp_input = ComponentVector(; Dict(k => itp[k](t) for k in keys(itp))...)
            # 遍历Unit中所有的Element进行求解
            for tmp_ele in unit.elements
                # 判断是否为ODEsElement
                if isa(tmp_ele, ODEElement)
                    # 计算出各个flux，更新至tmp_input中
                    tmp_fluxes = get_fluxes(tmp_ele, state=u, input=tmp_input)
                    # 求解du并更新du
                    tmp_du = tmp_ele.get_du(tmp_fluxes, tmp_ele.parameters)
                    for k in tmp_ele.state_names
                        du[k] = tmp_du[k]
                    end
                else
                    # 计算出各个flux，更新至tmp_input中
                    tmp_fluxes = get_fluxes(tmp_ele, input=tmp_input)
                end
                tmp_input = ComponentVector(tmp_input; tmp_fluxes...)
            end
        end

        # *求解这个函数
        prob = ODEProblem(ode_function!, u_init, (xs[1], maximum(xs)))
        sol = solve(
            prob,
            get(kwargs, :alg, BS3()),
            saveat=xs,
            reltol=get(kwargs, :reltol, 1e-3),
            abstol=get(kwargs, :abstol, 1e-3),
            sensealg=get(kwargs, :sensealg, ForwardDiffSensitivity())
        )
        solved_u = sol.u
        solved_u_matrix = hcat(solved_u...)
        solved_u = ComponentVector(; Dict(nm => solved_u_matrix[idx, :] for (idx, nm) in enumerate(keys(solved_u[1])))...)
        for ele in unit.elements
            ele.states = solved_u[collect(ele.state_names)]
        end
        unit.fluxes = input
        # *带入求解的结果计算最终的输出结果
        for tmp_ele in unit.elements
            if isa(tmp_ele, ODEElement)
                tmp_fluxes = get_fluxes(tmp_ele, state=solved_u, input=unit.fluxes)
            else
                tmp_fluxes = get_fluxes(tmp_ele, input=unit.fluxes)
            end
            update_fluxes!(unit, fluxes=tmp_fluxes)
        end
    end
    return get_fluxes(unit, flux_names=unit.output_names)
end
