# Test building model from YAML and loading parameters

using YAML

yaml_file = joinpath(dirname(dirname(@__DIR__)), "examples", "exphydro.yaml")

@testset "Load model from YAML" begin
    model = load_model_from_yaml(yaml_file)
    @test model isa HydroModels.HydroModel
    @test model.name == :exphydro

    # Verify model structure
    @test length(model.components) == 2
    @test :temp in HydroModels.get_input_names(model)
    @test :prcp in HydroModels.get_input_names(model)
    @test :lday in HydroModels.get_input_names(model)
    @test :snowpack in HydroModels.get_state_names(model)
    @test :soilwater in HydroModels.get_state_names(model)
end

@testset "Load config from YAML" begin
    config = load_config_from_yaml(yaml_file)
    @test config isa HydroModels.HydroConfig
    @test config.solver == HydroModels.MutableSolver
end

@testset "Load parameters from YAML" begin
    params_info = load_parameters_from_yaml(yaml_file)

    @test haskey(params_info, :Tmin)
    @test haskey(params_info, :Tmax)
    @test haskey(params_info, :Df)
    @test haskey(params_info, :Smax)
    @test haskey(params_info, :Qmax)
    @test haskey(params_info, :f)

    # Check metadata
    @test params_info[:Smax][:description] == "Maximum soil moisture"
    @test params_info[:Smax][:default] == 1709.0
    @test params_info[:Smax][:bounds] == (50.0, 5000.0)
end

@testset "Run YAML-loaded model" begin
    model = load_model_from_yaml(yaml_file)
    params_info = load_parameters_from_yaml(yaml_file)

    # Build parameter vector from YAML defaults
    param_names = HydroModels.get_param_names(model)
    param_values = NamedTuple{Tuple(param_names)}(
        params_info[p][:default] for p in param_names
    )
    pas = ComponentVector(params=param_values)

    # Load test data and run
    ts = collect(1:100)
    input_ntp, input_mat, df = load_test_data(:exphydro, ts)

    config = create_test_config(solver=MutableSolver, timeidx=ts)
    initstates = ComponentVector(snowpack=0.0, soilwater=1303.0)
    result = model(input_mat, pas, config; initstates=initstates)

    # Basic shape and sanity checks
    expected_outputs = length(HydroModels.get_state_names(model)) +
                       length(HydroModels.get_output_names(model))
    @test size(result, 1) == expected_outputs
    @test size(result, 2) == length(ts)
    @test all(isfinite, result)
end
