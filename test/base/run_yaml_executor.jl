"""
Test YAML Executor functionality

This test demonstrates the YAML executor that loads model configuration,
forcing data, and parameters from YAML files and executes the model.
"""

using HydroModels
using YAML
using CSV
using DataFrames
using ComponentArrays
using Test

@testset "YAML Executor Tests" begin
    # Create temporary directory for test files
    test_dir = mktempdir()

    @testset "CSV Data Executor" begin
        # Create test forcing data CSV
        forcing_csv = joinpath(test_dir, "forcing.csv")
        forcing_df = DataFrame(
            date = 1:100,
            Temp = 20.0 .+ 5.0 .* sin.(range(0, 2π, length=100)),
            Prcp = max.(0.0, 10.0 .* rand(100)),
            PET = 3.0 .+ 1.0 .* sin.(range(0, 2π, length=100))
        )
        CSV.write(forcing_csv, forcing_df)

        # Create test parameter CSV
        params_csv = joinpath(test_dir, "params.csv")
        params_df = DataFrame(
            Smax = [100.0],
            Qmax = [10.0],
            f = [2.0],
            Df = [0.5]
        )
        CSV.write(params_csv, params_df)

        # Create YAML configuration
        yaml_config = joinpath(test_dir, "model_config.yaml")
        yaml_content = """
version: "1.0"
schema: "hydromodels"

parameters:
  Smax:
    description: "Maximum soil moisture"
    units: "mm"
    default: 100.0
    bounds: [50.0, 500.0]
  Qmax:
    description: "Maximum baseflow"
    units: "mm/day"
    default: 10.0
    bounds: [1.0, 50.0]
  f:
    description: "Baseflow exponent"
    units: "-"
    default: 2.0
    bounds: [0.1, 5.0]
  Df:
    description: "Drainage factor"
    units: "-"
    default: 0.5
    bounds: [0.0, 1.0]

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "baseflow ~ Qmax * exp(-f * max(0.0, (Smax - soilwater) / Smax))"
      - formula: "evap ~ min(PET, soilwater)"
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
  path: "$forcing_csv"
  variables:
    Prcp: Prcp
    PET: PET
  time_column: date

parameters_data:
  type: csv
  path: "$params_csv"

initial_states:
  soilwater: 50.0
"""
        write(yaml_config, yaml_content)

        # Test execution
        @testset "Execute from YAML" begin
            output = execute_from_yaml(yaml_config)

            @test output isa NamedTuple
            @test haskey(output, :soil)
            @test size(output.soil, 1) == 100  # 100 timesteps
        end

        @testset "Execute with components return" begin
            output, model, config, data = execute_from_yaml(yaml_config, return_components=true)

            @test output isa NamedTuple
            @test model isa HydroModel
            @test config isa HydroConfig
            @test data isa Dict
            @test haskey(data, :Prcp)
            @test haskey(data, :PET)
        end
    end

    @testset "Direct Parameter Values" begin
        # Create YAML with direct parameter values
        yaml_config = joinpath(test_dir, "model_direct_params.yaml")

        forcing_csv = joinpath(test_dir, "forcing.csv")

        yaml_content = """
version: "1.0"
schema: "hydromodels"

parameters:
  Smax:
    description: "Maximum soil moisture"
    units: "mm"
    default: 100.0
  Qmax:
    description: "Maximum baseflow"
    units: "mm/day"
    default: 10.0

components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "baseflow ~ Qmax * (soilwater / Smax)"
    state_fluxes:
      - formula: "soilwater ~ Prcp - baseflow"

model:
  type: HydroModel
  name: simple_model
  components: [soil]

config:
  solver: MutableSolver
  interpolator: DirectInterpolation

data:
  type: csv
  path: "$forcing_csv"
  variables:
    Prcp: Prcp
    PET: PET

parameters_data:
  type: values
  values:
    Smax: 120.0
    Qmax: 8.0

initial_states:
  soilwater: 60.0
"""
        write(yaml_config, yaml_content)

        output = execute_from_yaml(yaml_config)
        @test output isa NamedTuple
        @test haskey(output, :soil)
    end

    @testset "Error Handling" begin
        # Test missing data section
        yaml_config = joinpath(test_dir, "model_no_data.yaml")
        yaml_content = """
version: "1.0"
parameters:
  Smax:
    default: 100.0
components:
  - type: HydroBucket
    name: soil
    fluxes:
      - formula: "baseflow ~ Qmax * soilwater"
    state_fluxes:
      - formula: "soilwater ~ Prcp - baseflow"
model:
  type: HydroModel
  name: test
  components: [soil]
"""
        write(yaml_config, yaml_content)

        @test_throws ErrorException execute_from_yaml(yaml_config)
    end

    # Cleanup
    rm(test_dir, recursive=true)
end

println("✓ YAML Executor tests completed")
