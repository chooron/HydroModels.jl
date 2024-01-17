# Element Methods
# get_id
function get_id(ele::BaseElement)::String
    return ele.id
end

# copy element
function copy_element(ele::ParameterizedElement)
    p = ele.parameters
    return typeof(ele)(ele.id, p)
end

function deepcopy_element(ele::ParameterizedElement)
    p = deepcopy(ele.parameters)
    return typeof(ele)(ele.id, p)
end

function copy_element(ele::StateElement)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, states)
end

function deepcopy_element(ele::StateElement)
    return copy_element(ele)
end

function copy_element(ele::StateElement)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, states)
end

function deepcopy_element(ele::StateElement)
    return copy_element(ele)
end

function copy_element(ele::Union{DiscElement,StateParameterizedElement})
    p = ele.parameters
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states)
end

function deepcopy_element(ele::Union{DiscElement,StateParameterizedElement})
    p = deepcopy(ele.parameters)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states)
end

function copy_element(ele::ODEsElement)
    p = ele.parameters
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states, ele.solver)
end

function deepcopy_element(ele::ODEsElement)
    p = deepcopy(ele.parameters)
    states = deepcopy(ele.states)
    return typeof(ele)(ele.id, p, states, ele.solver)
end


function get_parameters(ele::Union{ParameterizedElement,StateParameterizedElement}; names::Vector{String}=nothing)::Dict{String,Any}
    if isnothing(names)
        return ele.parameters
    else
        return Dict(name => ele.parameters[name] for name in names)
    end
end

function get_states(ele::Union{StateElement,StateParameterizedElement}; names::Vector{String}=nothing)::Dict{String,Vector{Number}}
    if isnothing(names)
        return ele.states
    else
        return Dict(name => ele.states[name] for name in names)
    end
end


function solve_prob(ele::ODEsElement, input::Dict{String,Any})::Vector{Number}
    dt = 1
    xs = 1:dt:length(input[keys(input)[1]])
    tspan = (times[1], times[end])

    # fit interpolation functions
    itp = Dict()
    for (key, value) in pairs(input)
        itp[key] = linear_interpolation(xs, value)
    end

    function func(u, p, t)
        # interpolate value by fitted functions
        tmp_input = Dict()
        for (key, value) in pairs(input)
            tmp_input[key] = value(t)
        end
        # return dt
        get_du(ele, u, tmp_input)
    end
    prob = ODEProblem(func, ele.init_states, tspan)
    sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6, saveat=dt)
    solved_u = sol.u
end


function get_output(ele::ODEsElement; input::Dict{String,Vector{Number}}, solved::Bool=true)::Dict{String,Vector{Number}}
    S = solve_prob(ele, input)
    fluxes = get_fluxes(ele, S, input)
end

# mutable struct BaseElement <: BaseElementType
#     id::String
# end

# function BaseElement(id::String)
#     BaseElement(id)
# end

# mutable struct ParameterizedElement <: ParameterizedElementType
#     id::String
#     parameters::Dict{String,Any}
# end

# function ParameterizedElement(id::String, parameters::Dict{String,Any})
#     ParameterizedElement(id, parameters)
# end


# mutable struct StateElement <: StateElementType
#     id::String
#     init_states::Dict{String,Number}

#     # 中间计算状态
#     states::Dict{String,Vector{Number}}
# end

# function StateElement(id::String, states::Dict{String,Number})
#     StateElement(id, states, nothing)
# end


# mutable struct StateParameterizedElement <: StateParameterizedElementType
#     id::String
#     parameters::Dict{String,Any}
#     init_states::Dict{String,Number}

#     # 中间计算状态
#     states::Dict{String,Vector{Number}}
# end

# function StateParameterizedElement(id::String, parameters::Dict{String,Any}, states::Dict{String,Number})
#     StateParameterizedElement(id, parameters, states, nothing)
# end


# function get_parameters(ele::Union{ParameterizedElement,StateParameterizedElement}; names::Vector{String}=nothing)
#     if isnothing(names)
#         return ele.parameters
#     else
#         return Dict(name => ele.parameters[name] for name in names)
#     end
# end


# function get_output(ele::BaseElement, input::Dict{String,Vector{Number}})
#     # 根据input输入到模块中获取fluxes
#     fluxes = get_fluxes(ele, input)
# end