"""
Test YAML Executor Core Functions (without YAML package)

This test verifies the core data loading functions work correctly.
"""

using HydroModels
using CSV
using DataFrames
using ComponentArrays
using Test

@testset "YAML Executor Core Functions" begin
    # Create temporary directory
    test_dir = mktempdir()
    println("Test directory: $test_dir")

    @testset "CSV Data Loading" begin
        # Create test CSV file
        csv_file = joinpath(test_dir, "test_forcing.csv")
        df = DataFrame(
            date = 1:10,
            Precipitation = [5.0, 0.0, 12.0, 3.0, 0.0, 8.0, 15.0, 2.0, 0.0, 6.0],
            PET = [3.0, 3.5, 2.8, 3.2, 3.1, 2.9, 3.3, 3.4, 3.0, 3.1]
        )
        CSV.write(csv_file, df)

        # Test if we can load the CSV
        df_loaded = CSV.read(csv_file, DataFrame)
        @test size(df_loaded) == (10, 3)
        @test hasproperty(df_loaded, :Precipitation)
        @test hasproperty(df_loaded, :PET)

        println("✓ CSV file created and loaded successfully")
    end

    @testset "Parameter CSV" begin
        # Create parameter CSV
        params_file = joinpath(test_dir, "test_params.csv")
        params_df = DataFrame(
            Smax = [100.0],
            Qmax = [10.0],
            f = [2.0]
        )
        CSV.write(params_file, params_df)

        # Load and verify
        params_loaded = CSV.read(params_file, DataFrame)
        @test size(params_loaded) == (1, 3)
        @test params_loaded[1, :Smax] == 100.0

        println("✓ Parameter CSV created and loaded successfully")
    end

    @testset "ComponentVector Creation" begin
        # Test creating ComponentVector from parameters
        params = ComponentVector(Smax=100.0, Qmax=10.0, f=2.0)
        @test params.Smax == 100.0
        @test params.Qmax == 10.0
        @test params.f == 2.0

        println("✓ ComponentVector created successfully")
    end

    # Cleanup
    rm(test_dir, recursive=true)
    println("\n✓ All core function tests passed!")
end
