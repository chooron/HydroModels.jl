get_input_names(func::AbstractFlux) = func.infos[:input]
get_input_names(ele::AbstractBucket) = ele.infos[:input]
get_input_names(route::AbstractRoute) = route.infos[:input]
get_input_names(unit::AbstractModel) = unit.infos[:input]

get_output_names(func::AbstractFlux) = func.infos[:output]
get_output_names(::AbstractStateFlux) = Symbol[]
get_output_names(ele::AbstractBucket) = ele.infos[:output]
get_output_names(route::AbstractRoute) = route.infos[:output]
get_output_names(unit::AbstractModel) = unit.infos[:output]

get_state_names(flux::AbstractFlux) = flux.infos[:state]
get_state_names(fluxes::Vector{<:AbstractFlux}) = reduce(union, get_state_names.(fluxes))
get_state_names(ele::AbstractBucket) = ele.infos[:state]
get_state_names(route::AbstractRoute) = route.infos[:state]
get_state_names(unit::AbstractModel) = unit.infos[:state]

get_param_names(func::AbstractFlux) = func.infos[:param]
get_param_names(::AbstractNeuralFlux) = Symbol[]
get_param_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_param_names.(funcs))
get_param_names(ele::AbstractBucket) = ele.infos[:param]
get_param_names(route::AbstractRoute) = route.infos[:param]
get_param_names(unit::AbstractModel) = unit.infos[:param]

get_nn_names(::AbstractFlux) = Symbol[]
get_nn_names(func::AbstractNeuralFlux) = func.infos[:nn]
get_nn_names(funcs::Vector{<:AbstractFlux}) = reduce(union, get_nn_names.(funcs))
get_nn_names(ele::AbstractBucket) = ele.infos[:nn]
get_nn_names(route::AbstractRoute) = route.infos[:nn]
get_nn_names(unit::AbstractModel) = unit.infos[:nn]

get_var_names(func::AbstractFlux) = get_input_names(func), get_output_names(func), get_state_names(func)
function get_var_names(funcs::Vector{<:AbstractFlux}, dfuncs::Vector{<:AbstractFlux})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = get_state_names(dfuncs)
    for func in vcat(funcs, dfuncs)
        #* 输入需要排除已有的输出变量，表明这个输入是中间计算得到的，此外state也不作为输入
        tmp_input_names = setdiff(setdiff(get_input_names(func), output_names), state_names)
        #* 排除一些输出，比如在flux中既作为输入又作为输出的变量，这时候他仅能代表输入
        tmp_output_names = setdiff(get_output_names(func), input_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
    end
    input_names, output_names, state_names
end
#* elements name utils
get_var_names(ele::AbstractBucket) = ele.infos[:input], ele.infos[:output], ele.infos[:state]
get_var_names(route::AbstractRoute) = route.infos[:input], route.infos[:output], route.infos[:state]

function get_var_names(components::Vector{<:AbstractComponent})
    input_names = Vector{Symbol}()
    output_names = Vector{Symbol}()
    state_names = Vector{Symbol}()
    for com in components
        tmp_input_names = get_input_names(com)
        tmp_output_names = get_output_names(com)
        tmp_state_names = get_state_names(com)
        tmp_input_names = setdiff(tmp_input_names, output_names)
        tmp_input_names = setdiff(tmp_input_names, state_names)
        union!(input_names, tmp_input_names)
        union!(output_names, tmp_output_names)
        union!(state_names, tmp_state_names)
    end
    input_names, output_names, state_names
end

get_var_names(unit::AbstractModel) = unit.infos[:var]