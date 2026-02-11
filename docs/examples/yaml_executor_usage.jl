"""
YAML Executor Usage Example

This script demonstrates how to use the YAML executor to run hydrological models
with configuration and data specified in YAML files.
"""

using HydroModels
using YAML
using CSV
using DataFrames
using Plots

# ============================================
# Example 1: Create sample data and execute
# ============================================

function example_csv_execution()
    println("=" ^ 60)
    println("Example 1: CSV Data Execution")
    println("=" ^ 60)

    # Create sample forcing data
    n_days = 365
    forcing_df = DataFrame(
        date = 1:n_days,
        Precipitation = max.(0.0, 5.0 .+ 10.0 .* rand(n_days)),
        PotentialET = 3.0 .+ 2.0 .* sin.(range(0, 4π, length=n_days))
    )

    # Save to CSV
    CSV.write("forcing_data.csv", forcing_df)
    println("✓ Created forcing_data.csv with $n_days days")

    # Create sample parameters
    params_df = DataFrame(
        Smax = [150.0],
        Qmax = [12.0],
        f = [2.5],
        Df = [0.6]
    )
    CSV.write("parameters.csv", params_df)
    println("✓ Created parameters.csv")

    # Create YAML configuration
    yaml_content = """
version: "1.0"
schema: "hydromodels"

parameters:
  Smax:
    description: "Maximum soil moisture"
    units: "mm"
    default: 150.0
  Qmax:
    description: "Maximum baseflow"
    units: "mm/day"
    default: 12.0
  f:
    description: "Baseflow exponent"
    units: "-"
    default: 2.5
  Df:
    description: "Drainage factor"
    units: "-"
    default: 0.6

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "baseflow ~ Qmax * exp(-f * max(0.0, (Smax - soilwater) / Smax))"
      - formula: "evap ~ min(PET, soilwater * Df)"
    state_fluxes:
      - formula: "soilwater ~ Prcp - evap - baseflow"

model:
  type: HydroModel
  name: simple_bucket
  components: [soil]

config:
  solver: MutableSolver
  interpolator: DirectInterpolation
  min_value: 1.0e-6

data:
  type: csv
  path: "forcing_data.csv"
  variables:
    Prcp: Precipitation
    PET: PotentialET
  time_column: date

parameters_data:
  type: csv
  path: "parameters.csv"

initial_states:
  soilwater: 75.0
"""

    write("model_config.yaml", yaml_content)
    println("✓ Created model_config.yaml")

    # Execute model
    println("\nExecuting model...")
    output, model, config, data = execute_from_yaml("model_config.yaml", return_components=true)

    println("\n" * "=" ^ 60)
    println("Execution Results")
    println("=" ^ 60)
    println("Model: $(model.name)")
    println("Components: $(length(model.components))")
    println("Timesteps: $(size(output.soil, 1))")
    println("State variables: $(size(output.soil, 2))")

    # Extract results
    soilwater = output.soil[:, 1]
    baseflow = output.soil[:, 2]
    evap = output.soil[:, 3]

    println("\nSoil Water Statistics:")
    println("  Mean: $(round(mean(soilwater), digits=2)) mm")
    println("  Min:  $(round(minimum(soilwater), digits=2)) mm")
    println("  Max:  $(round(maximum(soilwater), digits=2)) mm")

    println("\nBaseflow Statistics:")
    println("  Mean: $(round(mean(baseflow), digits=2)) mm/day")
    println("  Total: $(round(sum(baseflow), digits=2)) mm")

    println("\nEvapotranspiration Statistics:")
    println("  Mean: $(round(mean(evap), digits=2)) mm/day")
    println("  Total: $(round(sum(evap), digits=2)) mm")

    # Plot results
    try
        p1 = plot(forcing_df.date, forcing_df.Precipitation,
                  label="Precipitation", ylabel="mm/day", title="Forcing Data",
                  legend=:topright)
        plot!(p1, forcing_df.date, forcing_df.PotentialET,
              label="PET")

        p2 = plot(1:n_days, soilwater,
                  label="Soil Water", ylabel="mm", title="Model States",
                  legend=:topright, color=:blue)

        p3 = plot(1:n_days, baseflow,
                  label="Baseflow", ylabel="mm/day", title="Model Fluxes",
                  legend=:topright, color=:green)
        plot!(p3, 1:n_days, evap,
              label="Evapotranspiration", color=:red)

        plot(p1, p2, p3, layout=(3, 1), size=(800, 800))
        savefig("yaml_executor_results.png")
        println("\n✓ Saved plot to yaml_executor_results.png")
    catch e
        println("\n⚠ Could not create plot: $e")
    end

    # Cleanup
    rm("forcing_data.csv")
    rm("parameters.csv")
    rm("model_config.yaml")

    return output
end

# ============================================
# Example 2: Direct parameter values
# ============================================

function example_direct_parameters()
    println("\n" * "=" ^ 60)
    println("Example 2: Direct Parameter Values")
    println("=" ^ 60)

    # Create forcing data
    n_days = 100
    forcing_df = DataFrame(
        date = 1:n_days,
        Prcp = max.(0.0, 8.0 .* rand(n_days)),
        ET = 2.5 .+ 0.5 .* rand(n_days)
    )
    CSV.write("forcing_simple.csv", forcing_df)

    # YAML with direct parameter values
    yaml_content = """
version: "1.0"

parameters:
  Smax:
    default: 100.0
  k:
    default: 0.5

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "runoff ~ k * max(0.0, soilwater - Smax)"
    state_fluxes:
      - formula: "soilwater ~ Prcp - ET - runoff"

model:
  type: HydroModel
  name: simple_model
  components: [soil]

config:
  solver: MutableSolver
  interpolator: DirectInterpolation

data:
  type: csv
  path: "forcing_simple.csv"
  variables:
    Prcp: Prcp
    ET: ET

parameters_data:
  type: values
  values:
    Smax: 120.0
    k: 0.7

initial_states:
  soilwater: 60.0
"""

    write("simple_config.yaml", yaml_content)
    println("✓ Created simple_config.yaml with direct parameter values")

    # Execute
    output = execute_from_yaml("simple_config.yaml")

    println("\nExecution completed successfully!")
    println("Output shape: $(size(output.soil))")

    # Cleanup
    rm("forcing_simple.csv")
    rm("simple_config.yaml")

    return output
end

# ============================================
# Run examples
# ============================================

if abspath(PROGRAM_FILE) == @__FILE__
    println("HydroModels YAML Executor Examples\n")

    # Run Example 1
    output1 = example_csv_execution()

    # Run Example 2
    output2 = example_direct_parameters()

    println("\n" * "=" ^ 60)
    println("All examples completed successfully!")
    println("=" ^ 60)
end
