function Base.show(io::IO, flux::AbstractHydroFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "HydroFlux(")
        print(io, "inputs=[", join(get_input_names(flux), ","), "]")
        print(io, ", outputs=[", join(get_output_names(flux), ","), "]")
        print(io, ", params=[", join(get_param_names(flux), ","), "]")
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "HydroFlux", color=:light_blue, bold=true)
        printstyled(io, "{$(get_name(flux))}", color=:light_black)
        println(io)
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:  ", color=:light_green)
        println(io, "[", join(get_input_names(flux), ", "), "]")
        
        # Print outputs with a different color
        print(io, "│ ")
        printstyled(io, "Outputs: ", color=:light_yellow)
        println(io, "[", join(get_output_names(flux), ", "), "]")
        
        # Print parameters with a different color
        print(io, "│ ")
        printstyled(io, "Params:  ", color=:light_magenta)
        println(io, "[", join(get_param_names(flux), ", "), "]")
        
        if !isempty(flux.exprs)
            print(io, "│ ")
            printstyled(io, "Expressions:", color=:cyan)
            println(io)
            for (output, expr) in zip(get_output_names(flux), flux.exprs)
                print(io, "│   ")
                printstyled(io, "$output = ", color=:yellow)
                println(io, expr)
            end
        end
        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end

function Base.show(io::IO, flux::AbstractStateFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "StateFlux(")
        print(io, "inputs=[", join(get_input_names(flux), ","), "]")
        print(io, ", states=[", join(get_state_names(flux), ","), "]")
        print(io, ", params=[", join(get_param_names(flux), ","), "]")
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "StateFlux", color=:light_blue, bold=true)
        printstyled(io, "{$(get_name(flux))}", color=:light_black)
        println(io)
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:  ", color=:light_green)
        println(io, "[", join(get_input_names(flux), ", "), "]")
        
        # Print states with a different color
        print(io, "│ ")
        printstyled(io, "States:  ", color=:blue)
        println(io, "[", join(get_state_names(flux), ", "), "]")
        
        # Print parameters with a different color
        print(io, "│ ")
        printstyled(io, "Params:  ", color=:light_magenta)
        println(io, "[", join(get_param_names(flux), ", "), "]")
        
        if !isempty(get_exprs(flux))
            print(io, "│ ")
            printstyled(io, "Expressions:", color=:cyan)
            println(io)
            for (state, expr) in zip(get_state_names(flux), get_exprs(flux))
                print(io, "│   ")
                printstyled(io, "$state = ", color=:blue)
                println(io, expr)
            end
        end
        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end

function Base.show(io::IO, flux::AbstractNeuralFlux)
    compact = get(io, :compact, false)

    if compact
        print(io, "NeuralFlux(")
        print(io, "inputs=[", join(get_input_names(flux), ","), "]")
        print(io, ", outputs=[", join(get_output_names(flux), ","), "]")
        print(io, ", nns=[", join(get_nn_names(flux), ","), "]")
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "NeuralFlux", color=:light_blue, bold=true)
        printstyled(io, "{$(get_name(flux))}", color=:light_black)
        println(io)
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:  ", color=:light_green)
        println(io, "[", join(get_input_names(flux), ", "), "]")
        
        # Print outputs with a different color
        print(io, "│ ")
        printstyled(io, "Outputs: ", color=:light_yellow)
        println(io, "[", join(get_output_names(flux), ", "), "]")
        
        # Print neural networks with a different color
        print(io, "│ ")
        printstyled(io, "NNs:     ", color=:light_cyan)
        println(io, "[", join(get_nn_names(flux), ", "), "]")
        
        print(io, "│ ")
        printstyled(io, "Expressions:", color=:cyan)
        println(io)
        print(io, "│   ")
        printstyled(io, "$(get_output_names(flux)) = ", color=:yellow)
        println(io, "$(flux.chain)($(get_input_names(flux)))")

        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end

function Base.show(io::IO, uh::AbstractHydrograph)
    compact = get(io, :compact, false)
    if compact
        print(io, "UnitHydroFlux(")
        print(io, "inputs=[", join(get_input_names(uh), ","), "]")
        print(io, ", outputs=[", join(get_output_names(uh), ","), "]")
        print(io, ", params=[", join(get_param_names(uh), ","), "]")
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "UnitHydroFlux", color=:light_blue, bold=true)
        println(io)
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:       ", color=:light_green)
        println(io, "[", join(get_input_names(uh), ", "), "]")
        
        # Print outputs with a different color
        print(io, "│ ")
        printstyled(io, "Outputs:      ", color=:light_yellow)
        println(io, "[", join(get_output_names(uh), ", "), "]")
        
        # Print parameters with a different color
        print(io, "│ ")
        printstyled(io, "Parameters:   ", color=:light_magenta)
        println(io, "[", join(get_param_names(uh), ", "), "]")
        
        # todo 还没想好怎么打印这个
        # # Print UnitFunction with a different color
        # print(io, "│ ")
        # printstyled(io, "UnitFunction: ", color=:cyan)
        # println(io, typeof(uh.uhfunc).parameters[1])
        
        # # Print SolveType with a different color
        # print(io, "│ ")
        # printstyled(io, "SolveType:    ", color=:cyan)
        # println(io, typeof(uh).parameters[end])
        
        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end

function Base.show(io::IO, ele::AbstractBucket)
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroBucket(")
        print(io, "inputs=[", join(get_input_names(ele), ","), "]")
        print(io, ", states=[", join(get_state_names(ele), ","), "]")
        print(io, ", outputs=[", join(get_output_names(ele), ","), "]")
        print(io, ", params=[", join(get_param_names(ele), ","), "]")
        print(io, ", nns=[", join(get_nn_names(ele), ","), "]")
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "HydroBucket", color=:light_blue, bold=true)
        printstyled(io, "{$(get_name(ele))}", color=:light_black)
        println(io)
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:  ", color=:light_green)
        println(io, "[", join(get_input_names(ele), ", "), "]")
        
        # Print states with a different color
        print(io, "│ ")
        printstyled(io, "States:  ", color=:blue)
        println(io, "[", join(get_state_names(ele), ", "), "]")
        
        # Print outputs with a different color
        print(io, "│ ")
        printstyled(io, "Outputs: ", color=:light_yellow)
        println(io, "[", join(get_output_names(ele), ", "), "]")
        
        # Print parameters with a different color
        print(io, "│ ")
        printstyled(io, "Params:  ", color=:light_magenta)
        println(io, "[", join(get_param_names(ele), ", "), "]")
        
        # Print neural networks with a different color
        print(io, "│ ")
        printstyled(io, "NNs:     ", color=:light_cyan)
        println(io, "[", join(get_nn_names(ele), ", "), "]")

        # Print fluxes expressions
        println(io, "│")
        print(io, "│ ")
        printstyled(io, "Fluxes:", color=:yellow)
        println(io)
        for flux in ele.fluxes
            for (output, ex) in zip(get_outputs(flux), get_exprs(flux))
                print(io, "│   ")
                println(io, output ~ ex)
            end
        end

        # Print dfluxes expressions if any
        if !isempty(ele.dfluxes)
            println(io, "│")
            print(io, "│ ")
            printstyled(io, "State Fluxes:", color=:yellow)
            println(io)
            for dflux in ele.dfluxes
                for (st, ex) in zip(get_states(dflux), get_exprs(dflux))
                    print(io, "│   ")
                    println(io, st ~ ex)
                end
            end
        end
        
        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end

function Base.show(io::IO, route::AbstractHydroRoute)
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroRoute(")
        print(io, "inputs=[", join(get_input_names(route), ","), "]")
        print(io, ", states=[", join(get_state_names(route), ","), "]")
        print(io, ", outputs=[", join(get_output_names(route), ","), "]")
        print(io, ", params=[", join(get_param_names(route), ","), "]")
        print(io, ", nns=[", join(get_nn_names(route), ","), "]")
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "HydroBucket", color=:light_blue, bold=true)
        printstyled(io, "{$(get_name(route))}", color=:light_black)
        println(io)
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:  ", color=:light_green)
        println(io, "[", join(get_input_names(route), ", "), "]")
        
        # Print states with a different color
        print(io, "│ ")
        printstyled(io, "States:  ", color=:blue)
        println(io, "[", join(get_state_names(route), ", "), "]")
        
        # Print outputs with a different color
        print(io, "│ ")
        printstyled(io, "Outputs: ", color=:light_yellow)
        println(io, "[", join(get_output_names(route), ", "), "]")
        
        # Print parameters with a different color
        print(io, "│ ")
        printstyled(io, "Params:  ", color=:light_magenta)
        println(io, "[", join(get_param_names(route), ", "), "]")
        
        # Print neural networks with a different color
        print(io, "│ ")
        printstyled(io, "NNs:     ", color=:light_cyan)
        println(io, "[", join(get_nn_names(route), ", "), "]")

        # Print fluxes expressions
        println(io, "│")
        print(io, "│ ")
        printstyled(io, "Fluxes:", color=:yellow)
        println(io)
        for flux in route.rfluxes
            for (output, ex) in zip(get_outputs(flux), get_exprs(flux))
                print(io, "│   ")
                println(io, output ~ ex)
            end
        end

        println(io, "│")
        print(io, "│ ")
        printstyled(io, "State Fluxes:", color=:yellow)
        println(io)
        for dflux in route.dfluxes
            for (st, ex) in zip(get_states(dflux), get_exprs(dflux))
                print(io, "│   ")
                println(io, st ~ ex)
            end
        end
        
        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end

function Base.show(io::IO, model::AbstractModel)
    fluxes_in_model = filter(x -> x isa AbstractFlux, model.components)
    buckets_in_model = filter(x -> x isa AbstractBucket, model.components)
    routes_in_model = filter(x -> x isa AbstractRoute, model.components)
    uh_in_model = filter(x -> x isa AbstractHydrograph, model.components)
    @assert(length(routes_in_model) <= 1, "Only one route is allowed in a model")
    compact = get(io, :compact, false)
    if compact
        print(io, "HydroModel(")
        print(io, "name: ", get_name(model))
        print(io, ", components: ", length(model.components))
        print(io, ")")
    else
        printstyled(io, "┌ ", color=:light_blue, bold=true)
        printstyled(io, "HydroModel: ", color=:light_blue, bold=true)
        printstyled(io, get_name(model), color=:white, bold=true)
        println(io)
        
        # Print components with a different color
        print(io, "│ ")
        printstyled(io, "Components: ", color=:cyan)
        println(io, join(map(c -> get_name(c), model.components), ", "))
        
        # Print inputs with a different color
        print(io, "│ ")
        printstyled(io, "Inputs:     ", color=:light_green)
        println(io, "[", join(get_input_names(model), ", "), "]")
        
        # Print states with a different color
        print(io, "│ ")
        printstyled(io, "States:     ", color=:blue)
        println(io, "[", join(get_state_names(model), ", "), "]")
        
        # Print outputs with a different color
        print(io, "│ ")
        printstyled(io, "Outputs:    ", color=:light_yellow)
        println(io, "[", join(get_output_names(model), ", "), "]")
        
        # Print parameters with a different color
        print(io, "│ ")
        printstyled(io, "Parameters: ", color=:light_magenta)
        println(io, "[", join(get_param_names(model), ", "), "]")
        
        # Print component summary
        print(io, "│ ")
        printstyled(io, "Summary:", color=:white, bold=true)
        println(io)
        
        # Print fluxes count
        print(io, "│   ")
        printstyled(io, "Fluxes:  ", color=:light_cyan)
        println(io, length(fluxes_in_model), " flux", length(fluxes_in_model) == 1 ? "" : "es")
        
        # Print buckets count
        print(io, "│   ")
        printstyled(io, "Buckets: ", color=:light_cyan)
        println(io, length(buckets_in_model), " bucket", length(buckets_in_model) == 1 ? "" : "s")
        
        # Print routes count
        print(io, "│   ")
        printstyled(io, "Routes:  ", color=:light_cyan)
        println(io, length(routes_in_model), " route", length(routes_in_model) == 1 ? "" : "s")

        # Print uh count
        print(io, "│   ")
        printstyled(io, "UnitHydrographs:  ", color=:light_cyan)
        println(io, length(uh_in_model), " uh", length(uh_in_model) == 1 ? "" : "s")
        
        printstyled(io, "└─", color=:light_blue)
        println(io)
    end
end