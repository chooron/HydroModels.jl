"""
SymbolicUtils-based builders.

These builders keep the existing runtime API, but generate the function body through
`SymbolicUtils.Code` combinators instead of assembling the whole anonymous function by hand.

Notes:
- `SymbolicUtils.Code.toexpr` can emit array construction code, but it does not automatically
  broadcast scalar intrinsics such as `sin` over vector or matrix inputs.
- To preserve the current `build.jl` behavior, `Dim1` and `Dim2` wrap each symbolic RHS in `@.`.
- Neural flux embedding is intentionally left out for now.
"""

@inline function _symbolic_codegen_rhs(expr, ::NoBroadcast)
    return expr
end

@inline function _symbolic_codegen_rhs(expr, ::Broadcast)
    return LiteralExpr(:(@. $(simplify_expr(toexpr(expr)))))
end

@inline _symbolic_codegen_rhs(expr, dim::AbstractDimConfig) =
    _symbolic_codegen_rhs(expr, broadcast_strategy(dim))

@inline function _symbolic_var_assignments(
    vars::AbstractVector{Symbol},
    target::Symbol,
    dim_config::AbstractDimConfig,
    config::BuildConfig = DEFAULT_BUILD_CONFIG;
    prefix::String = ""
)
    return map(enumerate(vars)) do (idx, var)
        rhs = make_index_expr(target, idx, dim_config)
        var_name = Symbol(prefix, var)
        if config.mode == Fast
            Assignment(var_name, :(@inbounds $(rhs)))
        else
            Assignment(var_name, rhs)
        end
    end
end

@inline _symbolic_var_assignments(; vars, target, dims=0, prefix="", config=DEFAULT_BUILD_CONFIG) =
    _symbolic_var_assignments(vars, target, to_dim_type(Val(dims)), config; prefix=prefix)

@inline function _symbolic_param_assignments(
    params::AbstractVector{Symbol};
    target::Symbol = :pas
)
    return [Assignment(p, :($(target).params.$p)) for p in params]
end

@inline function _symbolic_all_assignments(
    infos::HydroModelCore.HydroInfos,
    dim_config::AbstractDimConfig,
    config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    return vcat(
        _symbolic_var_assignments(infos.inputs, :inputs, dim_config, config),
        _symbolic_var_assignments(infos.states, :states, dim_config, config),
        _symbolic_param_assignments(infos.params)
    )
end

@inline function _symbolic_flux_assignments(
    flux::AbstractHydroFlux,
    dim_config::AbstractDimConfig
)
    rhs_mode = broadcast_strategy(dim_config)
    return [
        Assignment(name, _symbolic_codegen_rhs(expr, rhs_mode))
        for (name, expr) in zip(get_output_names(flux), get_exprs(flux))
    ]
end

@inline function _symbolic_flux_assignments(
    fluxes::AbstractVector,
    dim_config::AbstractDimConfig
)
    return reduce(vcat, (_symbolic_flux_assignments(flux, dim_config) for flux in fluxes); init=Assignment[])
end

@inline function _symbolic_make_output_array(names::AbstractVector{Symbol})
    return MakeArray(collect(names), Array)
end

function _symbolic_state_return_expr(
    dfluxes::AbstractVector,
    dim_config::AbstractDimConfig
)
    exprs = reduce(vcat, get_exprs.(dfluxes); init=[])
    isempty(exprs) && return :(nothing)

    if uses_broadcast(dim_config)
        pieces = [:(@. $(simplify_expr(toexpr(expr)))) for expr in exprs]
        return LiteralExpr(:(vcat($(pieces...))))
    end

    return MakeArray(collect(exprs), Array)
end

function _symbolic_route_return_expr(
    infos::HydroModelCore.HydroInfos,
    dfluxes::AbstractVector
)
    exprs = reduce(vcat, get_exprs.(dfluxes); init=[])
    if isempty(exprs)
        return _symbolic_make_output_array(infos.outputs)
    end

    state_terms = [:(@. $(simplify_expr(toexpr(expr)))) for expr in exprs]
    return LiteralExpr(:([$(infos.outputs...), vcat($(state_terms...))]))
end

function _symbolic_build_runtime_function(
    args::Vector,
    assignments::Vector;
    result,
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG,
    enable_cse::Bool = true
)
    pre = build_config.inline_hints ? Any[:(Base.@_inline_meta)] : Any[]
    body_ir = Let(assignments, result, false)
    if enable_cse
        body_ir = cse(body_ir)
    end
    func_ir = Func(args, Assignment[], body_ir, pre)
    func_expr = toexpr(func_ir)

    if build_config.debug
        println("=" ^ 80)
        println("Generated SymbolicUtils Function (Mode: $(build_config.mode)):")
        println("=" ^ 80)
        println(func_expr)
        println("=" ^ 80)
    end

    return @RuntimeGeneratedFunction(func_expr)
end

@inline function _throw_neural_flux_not_supported()
    throw(ArgumentError(
        "SymbolicUtils builders currently support hydro fluxes only; neural network embedding is skipped in symbol_build.jl"
    ))
end

@inline function _ensure_hydro_fluxes_only(fluxes::AbstractVector)
    any(f -> f isa AbstractNeuralFlux, fluxes) && _throw_neural_flux_not_supported()
    return nothing
end

function build_symbolic_flux_func(
    exprs::Vector{Num},
    infos::HydroModelCore.HydroInfos;
    dims::Int = 0,
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG,
    enable_cse::Bool = true
)
    dim_config = to_dim_type(Val(dims))
    assignments = vcat(
        _symbolic_var_assignments(infos.inputs, :inputs, dim_config, build_config),
        _symbolic_param_assignments(infos.params),
        [
            Assignment(name, _symbolic_codegen_rhs(expr, dim_config))
            for (name, expr) in zip(infos.outputs, exprs)
        ]
    )

    return _symbolic_build_runtime_function(
        Any[:inputs, :pas],
        assignments;
        result=_symbolic_make_output_array(infos.outputs),
        build_config=build_config,
        enable_cse=enable_cse
    )
end

function build_symbolic_bucket_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos,
    multiply::Bool;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG,
    enable_cse::Bool = true
)
    _ensure_hydro_fluxes_only(fluxes)

    flux_dim = multiply ? Dim2() : Dim1()
    diff_dim = multiply ? Dim1() : Dim0()

    flux_assignments = vcat(
        _symbolic_all_assignments(infos, flux_dim, build_config),
        _symbolic_flux_assignments(fluxes, flux_dim)
    )
    flux_func = _symbolic_build_runtime_function(
        Any[:inputs, :states, :pas],
        flux_assignments;
        result=_symbolic_make_output_array(infos.outputs),
        build_config=build_config,
        enable_cse=enable_cse
    )

    if isempty(infos.states)
        return flux_func, (_...) -> nothing
    end

    diff_assignments = vcat(
        _symbolic_all_assignments(infos, diff_dim, build_config),
        _symbolic_flux_assignments(fluxes, diff_dim)
    )
    diff_func = _symbolic_build_runtime_function(
        Any[:inputs, :states, :pas],
        diff_assignments;
        result=_symbolic_state_return_expr(dfluxes, diff_dim),
        build_config=build_config,
        enable_cse=enable_cse
    )

    return flux_func, diff_func
end

function build_symbolic_route_func(
    fluxes::Vector{<:AbstractHydroFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG,
    enable_cse::Bool = true
)
    _ensure_hydro_fluxes_only(fluxes)

    flux_assignments = vcat(
        _symbolic_all_assignments(infos, Dim2(), build_config),
        _symbolic_flux_assignments(fluxes, Dim1())
    )
    flux_func = _symbolic_build_runtime_function(
        Any[:inputs, :states, :pas],
        flux_assignments;
        result=_symbolic_make_output_array(infos.outputs),
        build_config=build_config,
        enable_cse=enable_cse
    )

    if isempty(infos.states)
        return flux_func, (_...) -> nothing
    end

    diff_assignments = vcat(
        _symbolic_all_assignments(infos, Dim1(), build_config),
        _symbolic_flux_assignments(fluxes, Dim1())
    )
    diff_func = _symbolic_build_runtime_function(
        Any[:inputs, :states, :pas],
        diff_assignments;
        result=_symbolic_route_return_expr(infos, dfluxes),
        build_config=build_config,
        enable_cse=enable_cse
    )

    return flux_func, diff_func
end

function build_symbolic_uh_func(
    uh_conds::AbstractVector{<:Pair},
    params::AbstractVector{Symbol},
    max_lag::Number;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG,
    enable_cse::Bool = true
)
    param_assignments = _symbolic_param_assignments(params)
    bounds = vcat([0], reverse(first.(uh_conds)))
    values = reverse(last.(uh_conds))

    condition_exprs = map(eachindex(values)) do i
        quote
            if $(toexpr(bounds[i])) <= t <= $(toexpr(bounds[i + 1]))
                return $(toexpr(values[i]))
            end
        end
    end

    uh_func = _symbolic_build_runtime_function(
        Any[:t, :pas],
        param_assignments;
        result=LiteralExpr(quote
            $(condition_exprs...)
            1.0
        end),
        build_config=build_config,
        enable_cse=enable_cse
    )

    max_lag_func = _symbolic_build_runtime_function(
        Any[:pas],
        param_assignments;
        result=LiteralExpr(:(ceil($(toexpr(max_lag))))),
        build_config=build_config,
        enable_cse=enable_cse
    )

    return uh_func, max_lag_func
end
