"""
    @hydroflux(name, eqs...)

Create a `HydroFlux` from a set of equations, where:
- Left sides of equations are output variables
- Right sides contain input variables and parameters
- Parameters are identified using ModelingToolkit's `isparameter`
- `name` is an optional name for the flux (provide as first argument)

# Examples
```julia
@variables x, y, z
@parameters a, b

# Create a named flux with one equation
flux1 = @hydroflux :my_flux z = a*x + b*y

# Create a flux with multiple equations
flux2 = @hydroflux begin
    z₁ = a*x + b*y
    z₂ = x^2 + y^2
end
```
"""
macro hydroflux(args...)
    # Process arguments
    if length(args) == 1
        # No name provided, just equations
        name = nothing
        eqs_expr = args[1]
    else
        # First argument is the name
        name = args[1] 
        eqs_expr = args[2]
    end
    
    # Extract equations based on whether it's a block or single equation
    if Meta.isexpr(eqs_expr, :block)
        # Extract all equations from the block
        eqs_list = []
        for arg in eqs_expr.args
            if !(arg isa LineNumberNode) && !Meta.isexpr(arg, :line)
                push!(eqs_list, arg)
            end
        end
        
        # Create a new expression with all equations
        eqs_array = Expr(:vect)
        for eq in eqs_list
            push!(eqs_array.args, eq)
        end
        
        # Use the array of equations
        equations_expr = eqs_array
    else
        # For single equations, just wrap it in an array
        equations_expr = Expr(:vect, eqs_expr)
    end
    
    # Return the processed equations with common code
    return esc(quote
        let
            equations = $equations_expr
            
            # Process the equations
            processed_eqs = []
            for eq in equations
                if eq isa Equation
                    # Already an Equation object
                    push!(processed_eqs, eq)
                elseif Meta.isexpr(eq, :(=)) || Meta.isexpr(eq, :(~))
                    # Convert to Equation if it's a Julia expression
                    lhs = eq.args[1]
                    rhs = eq.args[2]
                    push!(processed_eqs, Symbolics.Equation(lhs, rhs))
                else
                    error("Expected equation (using = or ~), got: $eq")
                end
            end

            # Extract left and right sides
            lhs_terms = [eq.lhs for eq in processed_eqs]
            rhs_terms = [eq.rhs for eq in processed_eqs]

            # Extract outputs from left sides
            outputs = Num.(lhs_terms)   

            # Extract all variables from right sides
            all_vars = Set{Num}()
            for expr in rhs_terms
                union!(all_vars, get_variables(expr))
            end

            # Separate inputs and parameters
            inputs = Num.(filter(x -> !ModelingToolkit.isparameter(x), collect(all_vars)))
            params = Num.(filter(x -> ModelingToolkit.isparameter(x), collect(all_vars)))

            # Remove any output variables from inputs (they might appear on both sides)
            inputs = setdiff(inputs, outputs)

            # Create the flux
            HydroFlux(inputs, outputs, params, exprs=Num.(rhs_terms), name=$(name))
        end
    end)
end

"""
    @neuralflux(eq::Expr)

Create a `NeuralFlux` using the syntax: `output ~ chain(inputs)`, where:
- `output` is the output variable or a vector of output variables
- `~` is the separator
- `chain` is a Lux neural network chain
- `inputs` is a vector of input variables

# Examples
```julia
@variables x, y, z
chain = Chain(Dense(2 => 10, relu), Dense(10 => 1), name=:my_net)

# Create a neural flux with a single output
flux1 = @neuralflux z ~ chain([x, y])

# Create a neural flux with multiple outputs
chain2 = Chain(Dense(2 => 16, relu), Dense(16 => 2), name=:multi_net)
flux2 = @neuralflux [z₁, z₂] ~ chain2([x, y])
```
"""
macro neuralflux(args...)

    if length(args) == 1
        # No name provided, just equations
        name = nothing
        eqs_expr = args[1]
    else
        # First argument is the name
        name = args[1]
        eqs_expr = args[2]
    end

    if eqs_expr.head != :call || eqs_expr.args[1] != :~
        error("Expected equation in the form: output ~ chain(inputs)")
    end

    # Get left and right sides of the tilde
    lhs = eqs_expr.args[2]  # Output variable(s)
    rhs = eqs_expr.args[3]  # Chain call expression

    # Ensure the right side is a function call
    if !Meta.isexpr(rhs, :call)
        error("Right side must be a function call: chain([inputs])")
    end

    # Extract the chain and inputs
    chain_expr = rhs.args[1]

    # Check if inputs are provided as an array
    if length(rhs.args) != 2 || !Meta.isexpr(rhs.args[2], :vect)
        error("Inputs must be provided as an array: chain([x, y, ...])")
    end

    inputs_expr = rhs.args[2].args

    # Create the expression for the macro expansion
    return quote
        # Get the chain
        chain_obj = $(chain_expr)

        # Extract outputs
        outputs = $(lhs) isa Vector ? Num.($(lhs)) : [Num($(lhs))]

        # Extract inputs
        inputs = Num.([$(inputs_expr...)])

        # Create the NeuralFlux
        NeuralFlux(inputs, outputs, chain_obj)
    end
end

"""
    @stateflux_build(name, eq)

Create a `StateFlux` from an equation, where:
- Left side of the equation is the state variable
- Right side contains the expression for the state variable's change
- `name` is an optional name for the flux (provide as first argument)

# Examples
```julia
@variables x, y, z
@parameters a, b

# Create a named state flux
flux1 = @stateflux_build :my_flux z = a*x + b*y

# Create a state flux without a name
flux2 = @stateflux_build z = a*x + b*y
```
"""
macro stateflux(args...)
    # Process arguments
    if length(args) == 1
        # No name provided, just equation
        name = nothing
        eq_expr = args[1]
    else
        # First argument is the name
        name = args[1]
        eq_expr = args[2]
    end

    # Escape the entire expression to be evaluated in the caller's context
    return esc(quote
        let
            # Ensure the expression is an equation
            eq = if $eq_expr isa Equation
                # Already an Equation object
                $eq_expr
            elseif Meta.isexpr($eq_expr, :(=)) || Meta.isexpr($eq_expr, :(~))
                # Convert to Equation if it's a Julia expression
                lhs = $eq_expr.args[1]
                rhs = $eq_expr.args[2]
                Symbolics.Equation(lhs, rhs)
            else
                error("Expected equation (using = or ~), got: $($eq_expr)")
            end

            # Extract state variable from left side
            state = Num(eq.lhs)

            # Extract all variables from right side
            rhs_vars = get_variables(eq.rhs)

            # Separate inputs and parameters
            inputs = Num.(filter(x -> !isparameter(x), rhs_vars))
            params = Num.(filter(x -> isparameter(x), rhs_vars))

            # Create the StateFlux
            StateFlux(
                inputs,
                state,
                params,
                expr=Num(eq.rhs),
                name=$(name)
            )
        end
    end)
end

"""
    @hydrobucket(name, block)

Create a `HydroBucket` with the specified fluxes and state fluxes.

# Syntax
```julia
# With name
@hydrobucket :bucket_name begin
    @fluxes begin
        flux1
        flux2
        # ...
    end
    
    @dfluxes begin
        dflux1
        dflux2
        # ...
    end
end

# Without name
@hydrobucket begin
    @fluxes begin
        flux1
        flux2
        # ...
    end
    
    # @dfluxes is optional
end
```

# Arguments
- `name`: Optional symbol for naming the bucket
- `@fluxes`: Required section defining the flux components
- `@dfluxes`: Optional section defining the state flux components

# Returns
A `HydroBucket` instance with the specified configuration.

# Examples
```julia
# Create a bucket with fluxes and state fluxes
@hydrobucket :my_bucket begin
    @fluxes begin
        precipitation_flux
        evaporation_flux
        runoff_flux
    end
    
    @dfluxes begin
        storage_flux
    end
end
```
"""
macro hydrobucket(args...)
    # Parse arguments
    if length(args) == 1
        # No name provided, just block
        name = nothing
        block = args[1]
    else
        # First argument is the name
        name = args[1]
        block = args[2]
    end

    # Ensure the block is a block expression
    if !Meta.isexpr(block, :block)
        error("Expected a block of sections, got: $block")
    end

    # Initialize variables with default values
    fluxes_expr = :AbstractHydroFlux[]
    dfluxes_expr = :AbstractStateFlux[]

    # Process each section in the block
    for expr in block.args
        if Meta.isexpr(expr, :macrocall)
            section_name = expr.args[1]

            if section_name == Symbol("@fluxes")
                # Process @fluxes section
                section_block = expr.args[3]
                if !Meta.isexpr(section_block, :block)
                    error("Expected a block after @fluxes, got: $section_block")
                end

                # Extract flux components from the block
                flux_components = filter(ex -> !Meta.isexpr(ex, :line), section_block.args)
                fluxes_expr = :([$([esc(flux) for flux in flux_components]...)])

            elseif section_name == Symbol("@dfluxes")
                # Process @dfluxes section
                section_block = expr.args[3]
                if !Meta.isexpr(section_block, :block)
                    error("Expected a block after @dfluxes, got: $section_block")
                end

                # Extract state flux components from the block
                dflux_components = filter(ex -> !Meta.isexpr(ex, :line), section_block.args)
                dfluxes_expr = :([$([esc(dflux) for dflux in dflux_components]...)])
            end
        end
    end

    # Create the expression for the macro expansion
    return quote
        HydroBucket(
            fluxes=$(fluxes_expr),
            dfluxes=$(dfluxes_expr),
            name=$(name),
            sort_fluxes=true
        )
    end
end

function (chain::LuxCore.AbstractLuxLayer)(var::Vector{Num})
    #* assert the chain has a name
    @assert chain.name isa Symbol "Neural network chain should have a name with Symbol type"
    #* Get the neural network name (neural flux param name) and object
    chain_name = chain.name
    #* Initialize model parameter type for model parameter dimension definition
    init_params = ComponentVector(Lux.initialparameters(StableRNG(42), chain))
    chain_params = first(@parameters $chain_name[1:length(init_params)] = Vector(init_params))
    lazy_params = Symbolics.array_term((x, axes) -> ComponentVector(x, axes), chain_params, getaxes(init_params), size=size(chain_params))

    nn_input_name = Symbol(chain_name, :_input)
    nn_input = first(@variables $(nn_input_name)[1:length(inputs)])

    #* Constructing a calculation expression based on a neural network
    flux_expr = LuxCore.stateless_apply(chain, nn_input, lazy_params)

    return (chain=chain, input=var, expr=flux_expr, params=chain_params)
end