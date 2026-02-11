"""
Custom error types for HydroModels with descriptive messages.
"""

abstract type HydroModelsError <: Exception end

"""
    ConfigurationError

Error thrown when a component is configured incorrectly.

# Fields
- `component::Symbol`: Name of the component with the error
- `field::Symbol`: Field that has the configuration error
- `expected::String`: Description of what was expected
- `got::String`: Description of what was actually provided
"""
struct ConfigurationError <: HydroModelsError
    component::Symbol
    field::Symbol
    expected::String
    got::String
end

function Base.showerror(io::IO, e::ConfigurationError)
    print(io, "ConfigurationError in $(e.component): ")
    print(io, "Expected $(e.field) to be $(e.expected), but got $(e.got)")
end

"""
    DimensionMismatchError

Error thrown when array dimensions don't match expected dimensions.

# Fields
- `component::Symbol`: Name of the component with the error
- `expected_dims::Tuple`: Expected dimensions
- `got_dims::Tuple`: Actual dimensions received
- `context::String`: Additional context about the error
"""
struct DimensionMismatchError <: HydroModelsError
    component::Symbol
    expected_dims::Tuple
    got_dims::Tuple
    context::String
end

function Base.showerror(io::IO, e::DimensionMismatchError)
    print(io, "DimensionMismatchError in $(e.component): ")
    print(io, "Expected dimensions $(e.expected_dims), got $(e.got_dims)")
    if !isempty(e.context)
        print(io, "\nContext: $(e.context)")
    end
end

"""
    MacroSyntaxError

Error thrown when macro syntax is incorrect.

# Fields
- `macro_name::Symbol`: Name of the macro
- `expected::String`: Description of expected syntax
- `got_expr::String`: The expression that was provided
- `suggestion::String`: Suggestion for fixing the error
"""
struct MacroSyntaxError <: HydroModelsError
    macro_name::Symbol
    expected::String
    got_expr::String
    suggestion::String
end

function Base.showerror(io::IO, e::MacroSyntaxError)
    print(io, "MacroSyntaxError in @$(e.macro_name): ")
    print(io, "Expected $(e.expected)\n")
    print(io, "Got: $(e.got_expr)\n")
    if !isempty(e.suggestion)
        print(io, "Suggestion: $(e.suggestion)")
    end
end

"""
    ParameterError

Error thrown when there's an issue with model parameters.

# Fields
- `component::Symbol`: Name of the component
- `param_name::Symbol`: Name of the parameter with the issue
- `issue::String`: Description of the issue
"""
struct ParameterError <: HydroModelsError
    component::Symbol
    param_name::Symbol
    issue::String
end

function Base.showerror(io::IO, e::ParameterError)
    print(io, "ParameterError in $(e.component): ")
    print(io, "Parameter '$(e.param_name)' $(e.issue)")
end

"""
    VariableError

Error thrown when there's an issue with variables.

# Fields
- `component::Symbol`: Name of the component
- `var_name::Symbol`: Name of the variable with the issue
- `issue::String`: Description of the issue
"""
struct VariableError <: HydroModelsError
    component::Symbol
    var_name::Symbol
    issue::String
end

function Base.showerror(io::IO, e::VariableError)
    print(io, "VariableError in $(e.component): ")
    print(io, "Variable '$(e.var_name)' $(e.issue)")
end
