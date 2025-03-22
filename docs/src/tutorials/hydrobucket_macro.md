# HydroBucket Macro

The `@hydrobucket` macro provides a convenient way to create `HydroBucket` instances using a structured, block-based syntax similar to ModelingToolkit's `@mtkmodel`. This approach allows for clear organization of flux components and state flux components.

## Syntax

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

## Arguments

- `name`: Optional symbol for naming the bucket
- `@fluxes`: Required section defining the flux components
- `@dfluxes`: Optional section defining the state flux components

## Return Value

The macro returns a fully configured `HydroBucket` instance with the specified fluxes, state fluxes, and name. By default, fluxes are sorted automatically to ensure proper calculation order.

## Examples

### Basic Usage

```julia
using HydroModels
using ModelingToolkit, Symbolics

# Define variables and parameters
@variables P, ET, Q, S
@parameters a, b, c

# Create flux components using macros
precipitation_flux = @hydroflux_build :precipitation P = P
evaporation_flux = @hydroflux_build :evaporation ET = a * S
runoff_flux = @hydroflux_build :runoff Q = b * S^c

# Create a state flux for storage
storage_flux = @stateflux_build :storage S = P - ET - Q

# Create a bucket with all components
bucket = @hydrobucket :my_bucket begin
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

### Without a Name

If you don't specify a name, one will be automatically generated based on the hash of the bucket's metadata:

```julia
# Create a bucket without specifying a name
bucket = @hydrobucket begin
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

### Without State Fluxes

You can create a bucket with only flux components:

```julia
# Create a bucket with only flux components
bucket = @hydrobucket :flux_only begin
    @fluxes begin
        precipitation_flux
        evaporation_flux
        runoff_flux
    end
end
```

## Comparison with Direct Constructor

The `@hydrobucket` macro provides a more readable and structured way to create bucket instances compared to using the constructor directly:

### Using the Macro

```julia
bucket = @hydrobucket :my_bucket begin
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

### Using the Constructor Directly

```julia
bucket = HydroBucket(
    fluxes = [precipitation_flux, evaporation_flux, runoff_flux],
    dfluxes = [storage_flux],
    name = :my_bucket,
    sort_fluxes = true
)
```

## Implementation Details

The macro performs the following steps:
1. Parses the input arguments to extract the optional name and the block of sections
2. Processes each section (`@fluxes`, `@dfluxes`) to extract the components
3. Constructs a call to the `HydroBucket` constructor with the extracted values
4. Sets `sort_fluxes = true` by default to ensure proper calculation order

## See Also

- [`@hydroflux_build`](hydroflux_build_macro.md): Macro for creating flux components
- [`@stateflux_build`](stateflux_build_macro.md): Macro for creating state flux components
- [`@neuralflux`](neuralflux_macro.md): Macro for creating neural network flux components
