"""
YAML expression parser - converts formula strings to Symbolics expressions.

This module provides functionality to parse mathematical formulas from YAML
configuration files and convert them into Symbolics.jl expressions.
"""

"""
    parse_formula(formula_str::AbstractString, vars::Dict{Symbol,Num}, params::Dict{Symbol,Num})

Parse a formula string into a Symbolics equation.

# Arguments
- `formula_str`: Formula string in the form "lhs ~ rhs"
- `vars`: Dictionary mapping variable names to Symbolics variables
- `params`: Dictionary mapping parameter names to Symbolics parameters

# Returns
- Symbolics equation (lhs ~ rhs)

# Example
```julia
@variables Q P ET
@parameters k
vars = Dict(:Q => Q, :P => P, :ET => ET)
params = Dict(:k => k)
eq = parse_formula("Q ~ k * (P - ET)", vars, params)
```
"""
function parse_formula(formula_str::AbstractString, vars::Dict{Symbol,Num}, params::Dict{Symbol,Num})
    # Split on '~' to get LHS and RHS
    parts = split(formula_str, "~", limit=2)
    if length(parts) != 2
        error("Invalid formula syntax (missing '~'): $formula_str")
    end

    lhs_str = strip(parts[1])
    rhs_str = strip(parts[2])

    # Parse LHS (should be a single variable)
    lhs_sym = Symbol(lhs_str)
    if !haskey(vars, lhs_sym)
        error("Unknown variable in formula LHS: $lhs_sym")
    end
    lhs = vars[lhs_sym]

    # Parse RHS (expression)
    rhs = parse_expression_string(rhs_str, vars, params)

    # Create equation
    return lhs ~ rhs
end

"""
    parse_expression_string(expr_str::AbstractString, vars::Dict, params::Dict)

Parse an expression string into a Symbolics expression.

# Arguments
- `expr_str`: Expression string (e.g., "k * (P - ET)")
- `vars`: Dictionary of variables
- `params`: Dictionary of parameters

# Returns
- Symbolics expression
"""
function parse_expression_string(expr_str::AbstractString, vars::Dict{Symbol,Num}, params::Dict{Symbol,Num})
    # Normalize to String so SubString and other AbstractString types work.
    return parse_expression_string(String(expr_str), vars, params)
end

function parse_expression_string(expr_str::String, vars::Dict{Symbol,Num}, params::Dict{Symbol,Num})
    # Parse string to Julia expression
    try
        expr = Meta.parse(expr_str)
        # Replace symbols with Symbolics variables/parameters and evaluate
        return evaluate_symbolic_expr(expr, vars, params)
    catch e
        error("Failed to parse expression '$expr_str': $e")
    end
end

"""
    evaluate_symbolic_expr(expr, vars::Dict, params::Dict)

Evaluate a Julia expression with Symbolics variables/parameters.

# Arguments
- `expr`: Julia expression (from Meta.parse)
- `vars`: Dictionary of variables
- `params`: Dictionary of parameters

# Returns
- Symbolics expression
"""
function evaluate_symbolic_expr(expr, vars::Dict{Symbol,Num}, params::Dict{Symbol,Num})
    if expr isa Symbol
        # Check if it's a variable or parameter
        if haskey(vars, expr)
            return vars[expr]
        elseif haskey(params, expr)
            return params[expr]
        else
            # Could be a function name - return as is for now
            return expr
        end
    elseif expr isa Expr && expr.head == :call
        # Handle function calls
        func_sym = expr.args[1]
        args = [evaluate_symbolic_expr(arg, vars, params) for arg in expr.args[2:end]]

        # Get the actual function and apply it
        if func_sym isa Symbol
            if func_sym == :+
                return +(args...)
            elseif func_sym == :-
                return length(args) == 1 ? -(args[1]) : -(args...)
            elseif func_sym == :*
                return *(args...)
            elseif func_sym == :/
                return /(args...)
            elseif func_sym == :^
                return ^(args...)
            elseif func_sym == :max
                return max(args...)
            elseif func_sym == :min
                return min(args...)
            elseif func_sym == :exp
                return exp(args[1])
            elseif func_sym == :log
                return log(args[1])
            elseif func_sym == :sqrt
                return sqrt(args[1])
            elseif func_sym == :abs
                return abs(args[1])
            elseif func_sym == :tanh
                return tanh(args[1])
            elseif func_sym == :step_func
                # Call the step_func from HydroModels
                return HydroModels.step_func(args[1])
            else
                error("Unknown function: $func_sym")
            end
        else
            error("Invalid function call: $func_sym")
        end
    elseif expr isa Number
        # Literal number
        return expr
    else
        # Other types - try to handle them
        error("Unsupported expression type: $(typeof(expr)) - $expr")
    end
end

"""
    replace_symbols(expr, vars::Dict, params::Dict)

Recursively replace symbols in expression with Symbolics variables/parameters.

# Arguments
- `expr`: Julia expression (from Meta.parse)
- `vars`: Dictionary of variables
- `params`: Dictionary of parameters

# Returns
- Modified expression with symbols replaced
"""
function replace_symbols(expr, vars::Dict{Symbol,Num}, params::Dict{Symbol,Num})
    if expr isa Symbol
        # Check if it's a variable or parameter
        if haskey(vars, expr)
            return vars[expr]
        elseif haskey(params, expr)
            return params[expr]
        else
            # Could be a function name (min, max, exp, etc.) - leave as is
            return expr
        end
    elseif expr isa Expr
        # Recursively process expression arguments
        return Expr(expr.head, [replace_symbols(arg, vars, params) for arg in expr.args]...)
    else
        # Literal value (number, string, etc.)
        return expr
    end
end

"""
    extract_variables_from_formula(formula_str::AbstractString)

Extract variable names from a formula string.

# Arguments
- `formula_str`: Formula string

# Returns
- Set of variable symbols found in the formula
"""
function extract_variables_from_formula(formula_str::AbstractString)
    # Simple regex-based extraction
    # Matches identifiers (letters, numbers, underscores)
    matches = eachmatch(r"\b[a-zA-Z_][a-zA-Z0-9_]*\b", formula_str)
    vars = Set{Symbol}()

    # Known functions to exclude
    known_functions = Set([:min, :max, :exp, :log, :sqrt, :abs, :sin, :cos, :tan,
                           :tanh, :step_func, :ceil, :floor, :round])

    for m in matches
        sym = Symbol(m.match)
        if !(sym in known_functions)
            push!(vars, sym)
        end
    end

    return vars
end
