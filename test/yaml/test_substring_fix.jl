"""
Test expression parser with SubString

This test verifies that the expression parser can handle SubString inputs.
"""

using HydroModels
using YAML
using Test

@testset "Expression Parser SubString Test" begin
    # Create a simple YAML file
    yaml_content = """
version: "1.0"

parameters:
  k:
    default: 0.5

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "Q ~ k * S"
    state_fluxes:
      - formula: "S ~ P - Q"

model:
  type: HydroModel
  name: test_model
  components: [soil]

config:
  solver: MutableSolver
  interpolator: DirectInterpolation
"""

    # Write to temporary file
    yaml_file = tempname() * ".yaml"
    write(yaml_file, yaml_content)

    try
        # Try to load the model
        model = load_model_from_yaml(yaml_file)

        @test model isa HydroModel
        @test model.name == :test_model
        @test length(model.components) == 1

        println("✓ Model loaded successfully with SubString handling")
    catch e
        println("✗ Error loading model: $e")
        rethrow(e)
    finally
        # Cleanup
        rm(yaml_file, force=true)
    end
end

println("\n✓ SubString test completed!")
