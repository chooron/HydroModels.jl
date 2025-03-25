function build_flux_func(
    inputs::Vector{Num},
    outputs::Vector{Num},
    params::Vector{Num},
    exprs::Vector{Num},
)
    input_names, output_names = Symbolics.tosymbol.(inputs), Symbolics.tosymbol.(outputs)
    param_names = Symbolics.tosymbol.(params)
    flux_exprs = toexprv2.(unwrap.(exprs))
    input_assign_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    compute_calls = [:($o = $expr) for (o, expr) in zip(output_names, flux_exprs)]
    return_calls = :(return [$(output_names...)])
    flux_func_expr = :(function (inputs, pas)
        $(input_assign_calls...)
        $(params_assign_calls...)
        $(compute_calls...)
        $(return_calls)
    end)
    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    return generated_flux_func
end

function build_ele_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    meta::ComponentVector,
)
    input_names, output_names = tosymbol.(meta.inputs), tosymbol.(meta.outputs)
    state_names, param_names = tosymbol.(meta.states), tosymbol.(meta.params)
    nfluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)

    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(state_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($nn = pas.nns.$nn) for nn in [nflux.nninfos[:nns] for nflux in nfluxes]]
    define_calls = reduce(vcat, [input_define_calls, state_define_calls, params_assign_calls, nn_params_assign_calls])

    # varibles definitions expressions
    state_compute_calls, multi_state_compute_calls, flux_compute_calls, multi_flux_compute_calls = [], [], [], []
    for f in fluxes
        if f isa AbstractNeuralFlux
            append!(state_compute_calls, [:($(f.nninfos[:inputs]) = [$(get_input_names(f)...)]) for f in nfluxes])
            append!(multi_state_compute_calls, [:($(f.nninfos[:inputs]) = permutedims(reduce(hcat, [$(get_input_names(f)...)]))) for f in nfluxes])
            append!(flux_compute_calls, [:($(f.nninfos[:inputs]) = permutedims(reduce(hcat, [$(get_input_names(f)...)]))) for f in nfluxes])
            append!(multi_flux_compute_calls, [:($(f.nninfos[:inputs]) = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [$(get_input_names(f)...)]), (3, 1, 2))) for f in nfluxes])

            push!(state_compute_calls, :($(f.nninfos[:outputs]) = $(f.func)($(f.nninfos[:inputs]), $(f.nninfos[:nns]))))
            append!(state_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[$i]) for (i, nm) in enumerate(get_output_names(f))])

            push!(multi_state_compute_calls, :($(f.nninfos[:outputs]) = $(f.func)($(f.nninfos[:inputs]), $(f.nninfos[:nns]))))
            append!(multi_state_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[$i, :]) for (i, nm) in enumerate(get_output_names(f))])

            push!(flux_compute_calls, :($(f.nninfos[:outputs]) = $(f.func)($(f.nninfos[:inputs]), $(f.nninfos[:nns]))))
            append!(flux_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[$i, :]) for (i, nm) in enumerate(get_output_names(f))])

            push!(multi_flux_compute_calls, :($(f.nninfos[:outputs]) = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), $(f.func).(eachslice($(f.nninfos[:inputs]), dims=2), Ref($(f.nninfos[:nns])))), (1, 3, 2))))
            append!(multi_flux_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[$i, :, :]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($nm = $(toexpr(expr))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(multi_state_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(multi_flux_compute_calls, [:($nm = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    # Create return expressions with concrete values
    return_flux = :(return [$(output_names...)])
    return_state = :(return [$(map(expr -> :($(toexpr(expr))), reduce(vcat, get_exprs.(dfluxes)))...)])
    return_multi_state = :(return [$(map(expr -> :($(toexprv2(unwrap(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)])

    # Create fcuntion expression
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(flux_compute_calls...)
        $(return_flux)
    end)

    multi_flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(multi_flux_compute_calls...)
        $(return_flux)
    end)

    diff_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(state_compute_calls...)
        $(return_state)
    end)

    multi_diff_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(multi_state_compute_calls...)
        $(return_multi_state)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_multi_flux_func = @RuntimeGeneratedFunction(multi_flux_func_expr)
    generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
    generated_multi_diff_func = @RuntimeGeneratedFunction(multi_diff_func_expr)
    return generated_flux_func, generated_multi_flux_func, generated_diff_func, generated_multi_diff_func
end

function build_route_func(
    fluxes::AbstractVector{<:AbstractFlux},
    dfluxes::AbstractVector{<:AbstractStateFlux},
    rfunc::Function,
    meta::ComponentVector,
)
    input_names, output_names = tosymbol.(meta.inputs), tosymbol.(meta.outputs)
    state_names, param_names = tosymbol.(meta.states), tosymbol.(meta.params)
    nfluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)

    input_define_calls = [:($i = inputs[$idx]) for (idx, i) in enumerate(input_names)]
    state_define_calls = [:($s = states[$idx]) for (idx, s) in enumerate(state_names)]
    params_assign_calls = [:($p = pas.params.$p) for p in param_names]
    nn_params_assign_calls = [:($nn = pas.nns.$nn) for nn in [nflux.nninfos[:nns] for nflux in nfluxes]]
    define_calls = reduce(vcat, [input_define_calls, state_define_calls, params_assign_calls, nn_params_assign_calls])
    state_compute_calls, flux_compute_calls, = [], []
    for f in fluxes
        if f isa AbstractNeuralFlux
            append!(state_compute_calls, [:($(f.nninfos[:inputs]) = permutedims(reduce(hcat, [$(get_input_names(f)...)]))) for f in nfluxes])
            push!(state_compute_calls, :($(f.nninfos[:outputs]) = $(f.func)($(f.nninfos[:inputs]), $(f.nninfos[:nns]))))
            append!(state_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[$i, :]) for nm in get_output_names(f)])
            
            append!(flux_compute_calls, [:($(f.nninfos[:inputs]) = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), [$(get_input_names(f)...)]), (3, 1, 2))) for f in nfluxes])
            push!(flux_compute_calls, :($(f.nninfos[:outputs]) = permutedims(reduce((m1, m2) -> cat(m1, m2, dims=3), $(f.func).(eachslice($(f.nninfos[:inputs]), dims=2), Ref($(f.nninfos[:nns])))), (1, 3, 2))))
            append!(flux_compute_calls, [:($(nm) = $(f.nninfos[:outputs])[$i, :, :]) for (i, nm) in enumerate(get_output_names(f))])
        else
            append!(state_compute_calls, [:($(nm) = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
            append!(flux_compute_calls, [:($(nm) = $(toexprv2(unwrap(expr)))) for (nm, expr) in zip(get_output_names(f), f.exprs)])
        end
    end

    dfluxes_outflows = reduce(vcat, [tosymbol.(dflux.meta.outflows) for dflux in dfluxes])
    return_state = :(return [$(map(expr -> :($(toexprv2(unwrap(expr)))), reduce(vcat, get_exprs.(dfluxes)))...)], [$(dfluxes_outflows...)])
    # return_state = :(return [$(map((expr, out) -> :($(toexprv2(unwrap(expr))) .+ $(rfunc)($(out))), zip(dfluxes_exprs, dfluxes_outflows))...)])
    # return_state = :(return [$(map((expr, out) -> :($(expr) .+ rfunc($(out))), zip(dfluxes_exprs, dfluxes_outflows))...)])
    # return_state = :(return [$(map((expr, out) -> :($(expr) .+ rfunc($(out))), zip(dfluxes_exprs, dfluxes_outflows))...)])
    # Create function expression
    meta_exprs = [:(Base.@_inline_meta)]

    flux_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(flux_compute_calls...)
        $(:(return [$(output_names...)]))
    end)

    diff_func_expr = :(function (inputs, states, pas)
        $(meta_exprs...)
        $(define_calls...)
        $(state_compute_calls...)
        $(return_state)
    end)

    generated_flux_func = @RuntimeGeneratedFunction(flux_func_expr)
    generated_diff_func = @RuntimeGeneratedFunction(diff_func_expr)
    return generated_flux_func, generated_diff_func
end
