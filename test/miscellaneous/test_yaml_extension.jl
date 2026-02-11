# Test for YAML Extension
# This file tests that the HydroModelsYAMLExt extension loads correctly

using Test
using HydroModels

# Test that extension is not loaded without YAML
@testset "Extension not loaded without YAML" begin
    @test !isdefined(Main, :load_model_from_yaml)
end

# Load YAML to trigger extension
using YAML

# Test that extension functions are now available
@testset "Extension loaded with YAML" begin
    @test isdefined(HydroModelsYAMLExt, :load_model_from_yaml)
    @test isdefined(HydroModelsYAMLExt, :load_config_from_yaml)
    @test isdefined(HydroModelsYAMLExt, :load_parameters_from_yaml)
end

# Test loading a simple YAML model
@testset "Load YAML model" begin
    # Create a temporary YAML file
    yaml_content = """
    version: "1.0"
    schema: "hydromodels"

    parameters:
      k:
        description: "Runoff coefficient"
        units: "dimensionless"
        default: 0.5
        bounds: [0.1, 1.0]

    components:
      - type: HydroBucket
        name: simple_bucket
        fluxes:
          - formula: "Q ~ k * P"
        state_fluxes: []

    model:
      type: HydroModel
      name: test_model
      components: [simple_bucket]

    config:
      solver: MutableSolver
      interpolator: DirectInterpolation
      min_value: 1.0e-6
    """

    # Write to temporary file
    temp_file = tempname() * ".yaml"
    write(temp_file, yaml_content)

    try
        # Test loading model
        model = HydroModelsYAMLExt.load_model_from_yaml(temp_file)
        @test model isa HydroModels.HydroModel
        @test model.name == :test_model

        # Test loading config
        config = HydroModelsYAMLExt.load_config_from_yaml(temp_file)
        @test config isa HydroModels.HydroConfig
        @test config.solver == HydroModels.MutableSolver

        # Test loading parameters
        params_info = HydroModelsYAMLExt.load_parameters_from_yaml(temp_file)
        @test haskey(params_info, :k)
        @test params_info[:k][:description] == "Runoff coefficient"
        @test params_info[:k][:bounds] == (0.1, 1.0)

        println("✓ All YAML extension tests passed!")
    finally
        # Clean up
        rm(temp_file, force=true)
    end
end
