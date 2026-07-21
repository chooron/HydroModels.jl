"""
HydroModelCore Unified Build System
===================================

Merged build system combining the best of V1 and V2:
- Efficient code generation with broadcasting control
- Suitable for ForwardDiff and Mooncake-compatible generated code
- Type-stable dimension handling
- Optional performance optimizations

Key Features:
1. Three dimension modes: Dim0 (scalar), Dim1 (vector), Dim2 (matrix)
2. Automatic broadcasting based on dimension configuration
3. AD-friendly operations with predictable allocation behavior
4. Debug mode for code inspection
5. Performance modes: safe, fast, autodiff
"""

#==============================================================================
# Core Type System
==============================================================================#

"""Abstract dimension configuration for compile-time dispatch"""
abstract type AbstractDimConfig end

"""0D: scalar computation (no broadcasting)"""
struct Dim0 <: AbstractDimConfig end

"""1D: vector computation (broadcasting)"""
struct Dim1 <: AbstractDimConfig end

"""2D: matrix computation (broadcasting)"""
struct Dim2 <: AbstractDimConfig end

# Type conversions
@inline to_dim_type(::Val{0}) = Dim0()
@inline to_dim_type(::Val{1}) = Dim1()
@inline to_dim_type(::Val{2}) = Dim2()
@inline to_dim_val(::Type{Dim0}) = Val(0)
@inline to_dim_val(::Type{Dim1}) = Val(1)
@inline to_dim_val(::Type{Dim2}) = Val(2)

"""Get dimension rank for dispatch"""
@inline dim_rank(::Dim0) = 0
@inline dim_rank(::Dim1) = 1
@inline dim_rank(::Dim2) = 2

#==============================================================================
# Broadcasting Strategy System
==============================================================================#

"""Broadcasting strategy trait for fine-grained control"""
abstract type BroadcastStrategy end
struct NoBroadcast <: BroadcastStrategy end      # Scalar operations
struct Broadcast <: BroadcastStrategy end         # Vector/Matrix broadcasting

"""Map dimension type to broadcasting strategy"""
@inline broadcast_strategy(::Dim0) = NoBroadcast()
@inline broadcast_strategy(::Dim1) = Broadcast()
@inline broadcast_strategy(::Dim2) = Broadcast()

"""Check if dimension uses broadcasting"""
@inline uses_broadcast(dim::AbstractDimConfig) = broadcast_strategy(dim) isa Broadcast

#==============================================================================
# Performance Mode Configuration
==============================================================================#

"""
Performance mode for code generation.

Modes:
- `:safe` - Maximum safety (default)
- `:fast` - Performance-oriented indexing
- `:autodiff` - Optimized for automatic differentiation
"""
@enum PerformanceMode Safe Fast AutoDiff

"""
Configuration for code generation behavior.

Fields:
- `mode`: Performance mode (:safe, :fast, :autodiff)
- `debug`: Enable debug output
- `inline_hints`: Add @inline hints
"""
struct BuildConfig
    mode::PerformanceMode
    debug::Bool
    inline_hints::Bool
end

# Predefined configurations
const SAFE_CONFIG = BuildConfig(Safe, false, true)
const FAST_CONFIG = BuildConfig(Fast, false, true)
const AUTODIFF_CONFIG = BuildConfig(AutoDiff, false, true)
const DEBUG_CONFIG = BuildConfig(Safe, true, true)

# Default is safe mode
const DEFAULT_BUILD_CONFIG = SAFE_CONFIG

"""Check whether a build configuration preserves AD-friendly indexing."""
@inline is_ad_safe(config::BuildConfig) = config.mode != Fast

#==============================================================================
# Expression Utilities
==============================================================================#

"""
$(SIGNATURES)

Simplify expression for readability.
Converts `:((+)(x, y))` to `:(x + y)`.
"""
@inline simplify_expr(expr) = Meta.parse(string(expr))

"""
$(SIGNATURES)

Apply broadcasting to expression based on strategy.
Keeps the default and autodiff modes friendly to ForwardDiff and Mooncake.
"""
@inline function apply_broadcast(expr::Expr, strategy::BroadcastStrategy)
    if strategy isa Broadcast
        :(@. $(simplify_expr(expr)))
    else
        expr
    end
end

@inline apply_broadcast(expr::Expr, dim::AbstractDimConfig) =
    apply_broadcast(expr, broadcast_strategy(dim))

#==============================================================================
# Indexing Expression Generation
==============================================================================#

"""
$(SIGNATURES)

Generate array indexing expression for variable extraction.

# Examples
```julia
make_index_expr(:inputs, 1, Dim0())  # => :(inputs[1])
make_index_expr(:inputs, 1, Dim1())  # => :(inputs[1, :])
make_index_expr(:inputs, 1, Dim2())  # => :(inputs[1, :, :])
```
"""
@inline function make_index_expr(target::Symbol, idx::Int, ::Dim0)
    return :($(target)[$idx])
end

@inline function make_index_expr(target::Symbol, idx::Int, ::Dim1)
    return :($(target)[$idx, :])
end

@inline function make_index_expr(target::Symbol, idx::Int, ::Dim2)
    return :($(target)[$idx, :, :])
end

#==============================================================================
# Variable Assignment Generation
==============================================================================#

"""
$(SIGNATURES)

Generate variable assignment expressions.

# Arguments
- `vars`: Variable names to assign
- `target`: Source array symbol
- `dim_config`: Dimension configuration
- `config`: Build configuration (optional)
- `prefix`: Optional prefix for variable names

# AD behavior
Safe and autodiff modes avoid unnecessary indexing transformations.
In `:fast` mode, may use @inbounds for performance.

# Example
```julia
generate_var_assignments([:temp, :prcp], :inputs, Dim1())
# => [:(temp = inputs[1, :]), :(prcp = inputs[2, :])]
```
"""
@inline function generate_var_assignments(
    vars::AbstractVector{Symbol},
    target::Symbol,
    dim_config::AbstractDimConfig,
    config::BuildConfig = DEFAULT_BUILD_CONFIG;
    prefix::String = ""
)
    assignments = map(enumerate(vars)) do (idx, var)
        var_name = Symbol(prefix, var)
        index_expr = make_index_expr(target, idx, dim_config)

        # Fast mode uses @inbounds; safe modes keep ordinary indexing.
        if config.mode == Fast
            :(@inbounds $(var_name) = $(index_expr))
        else
            :($(var_name) = $(index_expr))
        end
    end

    return assignments
end

# Convenience methods
@inline generate_var_assignments(; vars, target, dims=0, prefix="", config=DEFAULT_BUILD_CONFIG) =
    generate_var_assignments(vars, target, to_dim_type(Val(dims)), config; prefix=prefix)

#==============================================================================
# Parameter and Neural Network Assignments
==============================================================================#

"""
$(SIGNATURES)

Generate parameter extraction from ComponentArray.
AD-safe in all modes.
"""
@inline function generate_param_assignments(
    params::AbstractVector{Symbol};
    target::Symbol = :pas
)
    return [:($p = $(target).params.$p) for p in params]
end

@inline generate_param_assignments(; params, target=:pas) =
    generate_param_assignments(params; target=target)

"""
$(SIGNATURES)

Generate neural network extraction expressions.
AD-safe in all modes.
"""
@inline function generate_nn_assignments(
    nn_fluxes::AbstractVector;
    target::Symbol = :pas
)
    if isempty(nn_fluxes)
        return Any[]
    end

    nn_names = [get_nn_names(f)[1] for f in nn_fluxes]
    return [:($nn = $(target).nns.$nn) for nn in nn_names]
end

@inline generate_nn_assignments(; nnfluxes, target=:pas) =
    generate_nn_assignments(nnfluxes; target=target)

#==============================================================================
# Unified Assignment Generation
==============================================================================#

"""
$(SIGNATURES)

Generate all variable, parameter, and neural network assignments.
This is the main entry point for creating variable bindings.
"""
@inline function generate_all_assignments(
    infos::HydroModelCore.HydroInfos,
    nn_fluxes::AbstractVector,
    dim_config::AbstractDimConfig,
    config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    return vcat(
        generate_var_assignments(infos.inputs, :inputs, dim_config, config),
        generate_var_assignments(infos.states, :states, dim_config, config),
        generate_param_assignments(infos.params),
        generate_nn_assignments(nn_fluxes)
    )
end

@inline generate_all_assignments(
    infos::HydroModelCore.HydroInfos,
    nn_fluxes::Vector,
    dims::Int,
    config::BuildConfig = DEFAULT_BUILD_CONFIG
) = generate_all_assignments(infos, nn_fluxes, to_dim_type(Val(dims)), config)

#==============================================================================
# Computation Expression Generation
==============================================================================#

"""
$(SIGNATURES)

Generate computation expressions for HydroFlux.
Automatically applies broadcasting based on dimension.
AD-safe in all modes.
"""
@inline function generate_flux_expression(
    flux::AbstractHydroFlux,
    dim_config::AbstractDimConfig
)
    strategy = broadcast_strategy(dim_config)
    output_names = get_output_names(flux)

    return [
        :($(nm) = $(apply_broadcast(toexpr(expr), strategy)))
        for (nm, expr) in zip(output_names, flux.exprs)
    ]
end

"""
$(SIGNATURES)

Generate computation expressions for NeuralFlux.
Handles input stacking and output extraction.
AD-safe in all modes.
"""
@inline function generate_flux_expression(
    flux::AbstractNeuralFlux,
    dim_config::AbstractDimConfig
)
    nn_name = get_nn_names(flux)[1]
    nn_input = Symbol(nn_name, "_input")
    nn_output = Symbol(nn_name, "_output")
    input_names = flux.infos.inputs
    output_names = get_output_names(flux)

    # Determine input stacking based on dimension
    input_expr = if dim_rank(dim_config) == 0
        :([$(input_names...)])
    else
        # Stack keeps the generated expression allocation-transparent to AD.
        :(stack([$(input_names...)], dims=1))
    end

    # Generate computation pipeline.
    computations = [
        # Normalize inputs
        :($(nn_input) = $(input_expr) |> $(flux.norm_func)),
        # Apply neural network
        :($(nn_output) = $(flux.chain_func)($(nn_input), $(nn_name))),
        # Extract outputs
        [
            :($(nm) = $(make_index_expr(nn_output, i, dim_config)))
            for (i, nm) in enumerate(output_names)
        ]...
    ]

    return computations
end

"""Fallback for generic AbstractFlux"""
@inline generate_flux_expression(flux::AbstractFlux, ::Val{D}) where D =
    generate_flux_expression(flux, to_dim_type(Val(D)))

"""
$(SIGNATURES)

Generate all computation expressions for a vector of fluxes.
"""
@inline function generate_compute_calls(
    fluxes::AbstractVector,
    dim_config::AbstractDimConfig
)
    return reduce(vcat,
        [generate_flux_expression(f, dim_config) for f in fluxes],
        init=Any[]
    )
end

@inline generate_compute_calls(; fluxes, dims=0) =
    generate_compute_calls(fluxes, to_dim_type(Val(dims)))

#==============================================================================
# State Differential Return Expression
==============================================================================#

"""
$(SIGNATURES)

Generate return expression for state differentials.
Automatically applies broadcasting based on dimension.
AD-safe in all modes.

# Arguments
- `dfluxes`: Vector of state flux components
- `broadcast`: Whether to use broadcasting
"""
@inline function generate_states_return_expr(
    dfluxes::AbstractVector,
    broadcast::Bool
)
    all_exprs = reduce(vcat, get_exprs.(dfluxes), init=[])

    if isempty(all_exprs)
        return :(return nothing)
    end

    # Apply broadcasting.
    if broadcast
        broadcast_exprs = [:(@. $(simplify_expr(toexpr(expr)))) for expr in all_exprs]
        # Use vcat to preserve the existing state layout.
        return :(return vcat($(broadcast_exprs...)))
    else
        plain_exprs = [toexpr(expr) for expr in all_exprs]
        # Use a vector literal for scalar state updates.
        return :(return [$(plain_exprs...)])
    end
end

# Convenience methods
@inline generate_states_return_expr(dfluxes::AbstractVector, dim::AbstractDimConfig) =
    generate_states_return_expr(dfluxes, uses_broadcast(dim))

@inline generate_states_expression(; dfluxes, broadcast=false) =
    generate_states_return_expr(dfluxes, broadcast)

#==============================================================================
# Function Builder Configuration
==============================================================================#

"""
Configuration for building a runtime function.

Fields:
- `signature`: Function signature expression
- `assignments`: Variable assignment expressions
- `computations`: Computation expressions
- `return_expr`: Return statement
"""
struct FunctionBuildConfig
    signature::Expr
    assignments::Vector{Any}
    computations::Vector{Any}
    return_expr::Expr
end

"""
$(SIGNATURES)

Build RuntimeGeneratedFunction from configuration.

# Arguments
- `config`: Function build configuration
- `build_config`: Build configuration for optimization

# AD behavior
Generated functions use ordinary indexing in `:safe` and `:autodiff` modes.
"""
function build_runtime_function(
    config::FunctionBuildConfig;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    # Build function body
    body_exprs = Any[]

    # Add inline hint if requested
    if build_config.inline_hints
        push!(body_exprs, :(Base.@_inline_meta))
    end

    # Add assignments
    append!(body_exprs, config.assignments)

    # Add computations
    append!(body_exprs, config.computations)

    # Add return
    push!(body_exprs, config.return_expr)

    func_body = Expr(:block, body_exprs...)
    func_expr = Expr(:function, config.signature, func_body)

    # Debug output if requested
    if build_config.debug
        println("=" ^ 80)
        println("Generated Function (Mode: $(build_config.mode)):")
        println("=" ^ 80)
        println(func_expr)
        println("=" ^ 80)
        if !is_ad_safe(build_config)
            println("⚠  Warning: fast indexing mode may be less AD-friendly")
        else
            println("✓ AD-friendly indexing")
        end
        println("=" ^ 80)
    end

    return @RuntimeGeneratedFunction(func_expr)
end

"""
$(SIGNATURES)

Build a pair of functions (flux and diff).

# Returns
- `(flux_func, diff_func)` if has states
- `(flux_func, dummy_func)` if no states
"""
function build_function_pair(
    flux_config::FunctionBuildConfig,
    diff_config::Union{FunctionBuildConfig, Nothing},
    has_states::Bool;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    flux_func = build_runtime_function(flux_config; build_config=build_config)

    if has_states && !isnothing(diff_config)
        diff_func = build_runtime_function(diff_config; build_config=build_config)
        return flux_func, diff_func
    else
        # Dummy function for stateless models.
        return flux_func, (_...) -> nothing
    end
end

#==============================================================================
# High-Level Function Builders
==============================================================================#

"""
$(TYPEDSIGNATURES)

Build simple flux function from symbolic expressions.

# Arguments
- `exprs`: Vector of symbolic expressions
- `infos`: HydroInfos containing variable metadata
- `build_config`: Build configuration (optional)

# Returns
RuntimeGeneratedFunction that computes flux values.

# AD behavior
Generated code is pure in all modes.
"""
function build_flux_func(
    exprs::Vector{Num},
    infos::HydroModelCore.HydroInfos;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    # Apply broadcasting to expressions.
    flux_exprs = [:(@. $(simplify_expr(toexpr(expr)))) for expr in exprs]

    config = FunctionBuildConfig(
        :((inputs, pas)),
        vcat(
            generate_var_assignments(infos.inputs, :inputs, Dim0(), build_config),
            generate_param_assignments(infos.params)
        ),
        [:($(o) = $(expr)) for (o, expr) in zip(infos.outputs, flux_exprs)],
        :(return [$(infos.outputs...)])
    )

    return build_runtime_function(config; build_config=build_config)
end

"""
$(TYPEDSIGNATURES)

Build bucket function for hydrological models.

# Arguments
- `fluxes`: Vector of flux components
- `dfluxes`: Vector of state flux components
- `infos`: HydroInfos containing metadata
- `multiply`: Whether to use multiple grid mode
- `build_config`: Build configuration (optional)

# Returns
`(flux_func, diff_func)` pair of runtime-generated functions.

# Broadcasting Strategy
- Single grid (`multiply=false`): flux uses Dim1, diff uses Dim0
- Multiple grids (`multiply=true`): flux uses Dim2, diff uses Dim1

# AD behavior
Safe and autodiff modes use ordinary indexing; fast mode uses @inbounds.
"""
function build_bucket_func(
    fluxes::Vector{<:AbstractFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos,
    multiply::Bool;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    nn_fluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)
    has_states = !isempty(infos.states)

    # Dimension selection based on multiply mode
    flux_dim = multiply ? Dim2() : Dim1()
    diff_dim = multiply ? Dim1() : Dim0()

    # Build flux function configuration
    flux_config = FunctionBuildConfig(
        :((inputs, states, pas)),
        generate_all_assignments(infos, nn_fluxes, flux_dim, build_config),
        generate_compute_calls(fluxes, flux_dim),
        :(return [$(infos.outputs...)])
    )

    # Build diff function configuration
    diff_config = if has_states
        FunctionBuildConfig(
            :((inputs, states, pas)),
            generate_all_assignments(infos, nn_fluxes, diff_dim, build_config),
            generate_compute_calls(fluxes, diff_dim),
            generate_states_return_expr(dfluxes, diff_dim)
        )
    else
        nothing
    end

    return build_function_pair(flux_config, diff_config, has_states; build_config=build_config)
end

"""
$(TYPEDSIGNATURES)

Build routing function for river network calculations.

# Key Features
- flux: 2D assignment + 1D computation
- diff: 1D assignment + 1D computation + special return format

# AD behavior
Safe and autodiff modes use ordinary indexing.
"""
function build_route_func(
    fluxes::Vector{<:AbstractHydroFlux},
    dfluxes::Vector{<:AbstractStateFlux},
    infos::HydroModelCore.HydroInfos;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    nn_fluxes = filter(f -> f isa AbstractNeuralFlux, fluxes)
    has_states = !isempty(infos.states)

    # Flux configuration
    flux_config = FunctionBuildConfig(
        :((inputs, states, pas)),
        generate_all_assignments(infos, nn_fluxes, Dim2(), build_config),
        generate_compute_calls(fluxes, Dim1()),
        :(return [$(infos.outputs...)])
    )

    # Diff configuration with special return format
    diff_config = if has_states
        state_diff_exprs = [:(@. $(simplify_expr(toexpr(expr))))
                           for expr in reduce(vcat, get_exprs.(dfluxes), init=[])]

        FunctionBuildConfig(
            :((inputs, states, pas)),
            generate_all_assignments(infos, nn_fluxes, Dim1(), build_config),
            generate_compute_calls(fluxes, Dim1()),
            # Use an array literal to preserve the result layout.
            :(return [$(infos.outputs...), vcat($(state_diff_exprs...))])
        )
    else
        nothing
    end

    return build_function_pair(flux_config, diff_config, has_states; build_config=build_config)
end

"""
$(TYPEDSIGNATURES)

Build Unit Hydrograph (UH) function pair.

# Generates
1. `uh_func(t, pas)`: Compute weight at time t
2. `max_lag_func(pas)`: Compute maximum lag time

# AD behavior
Uses pure functions.
"""
function build_uh_func(
    uh_conds::AbstractVector{<:Pair},
    params::AbstractVector{Symbol},
    max_lag::Number;
    build_config::BuildConfig = DEFAULT_BUILD_CONFIG
)
    # Reverse conditions and values for proper evaluation order
    conditions_rev = vcat([0], reverse(first.(uh_conds)))
    values_rev = reverse(last.(uh_conds))
    param_assigns = generate_param_assignments(params)

    # Generate piecewise conditions
    condition_checks = map(eachindex(values_rev)) do i
        :(
            if $(toexpr(conditions_rev[i])) ≤ t ≤ $(toexpr(conditions_rev[i+1]))
                return $(toexpr(values_rev[i]))
            end
        )
    end

    # UH weight function
    uh_func_expr = :(function (t, pas)
        $(param_assigns...)
        $(condition_checks...)
        return 1.0  # Default fallback
    end)

    # Max lag function
    max_lag_expr = :(function (pas)
        $(param_assigns...)
        return ceil($(toexpr(max_lag)))
    end)

    if build_config.debug
        println("=" ^ 80)
        println("UH Function:")
        println(uh_func_expr)
        println("\nMax Lag Function:")
        println(max_lag_expr)
        println("=" ^ 80)
    end

    return (
        @RuntimeGeneratedFunction(uh_func_expr),
        @RuntimeGeneratedFunction(max_lag_expr)
    )
end

#==============================================================================
# Debugging and Analysis Utilities
==============================================================================#

"""
$(SIGNATURES)

Preview generated function expression for debugging.
"""
function preview_function_expr(config::FunctionBuildConfig)
    func_body = Expr(:block,
        :(Base.@_inline_meta),
        config.assignments...,
        config.computations...,
        config.return_expr
    )

    func_expr = Expr(:function, config.signature, func_body)

    println("\n" * "=" ^ 80)
    println("Generated Function Expression")
    println("=" ^ 80)
    println(func_expr)
    println("=" ^ 80 * "\n")
end

"""
$(SIGNATURES)

Analyze generated code complexity.
"""
function analyze_generated_code(config::FunctionBuildConfig)
    stats = Dict{Symbol, Any}(
        :num_assignments => length(config.assignments),
        :num_computations => length(config.computations),
        :uses_broadcast => any(ex -> ex isa Expr && ex.head == :macrocall &&
                              ex.args[1] == Symbol("@."),
                              config.computations),
        :complexity => length(config.assignments) + length(config.computations),
        :ad_safe => true  # Generated code is AD-safe by default
    )

    println("\n" * "=" ^ 80)
    println("Code Analysis")
    println("=" ^ 80)
    for (key, value) in stats
        println("  $key: $value")
    end
    println("=" ^ 80 * "\n")

    return stats
end
