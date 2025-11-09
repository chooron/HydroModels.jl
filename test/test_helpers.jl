# Test Helper Functions
# Provides common utilities to reduce code duplication across test files

"""
    load_test_data(dataset::Symbol, timesteps::AbstractVector{Int})

Load test data from a specified dataset.

# Arguments
- `dataset`: Dataset identifier (:exphydro, :gr4j, :m50)
- `timesteps`: Vector of time indices to load

# Returns
- `input_ntp`: NamedTuple of input variables
- `input_mat`: Matrix of input data
- `df`: Original DataFrame
"""
function load_test_data(dataset::Symbol, timesteps::AbstractVector{Int})
    if dataset == :exphydro
        df = DataFrame(CSV.File("../data/exphydro/01013500.csv"))
        input_ntp = (
            lday = df[timesteps, "dayl(day)"],
            temp = df[timesteps, "tmean(C)"],
            prcp = df[timesteps, "prcp(mm/day)"]
        )
        input_mat = Matrix(reduce(hcat, collect(input_ntp[[:temp, :lday, :prcp]]))')
        return input_ntp, input_mat, df
        
    elseif dataset == :gr4j
        df = DataFrame(CSV.File("../data/gr4j/sample.csv"))
        for col in names(df)[3:end]
            df[ismissing.(df[:, col]), col] .= 0.0
        end
        input_ntp = (prcp = df[!, "prec"], ep = df[!, "pet"])
        input_mat = Matrix(reduce(hcat, collect(input_ntp[[:prcp, :ep]]))')
        return input_ntp, input_mat, df
        
    elseif dataset == :m50
        df = DataFrame(CSV.File("../data/m50/01013500.csv"))
        prcp_vec = df[timesteps, "Prcp"]
        temp_vec = df[timesteps, "Temp"]
        dayl_vec = df[timesteps, "Lday"]
        input_ntp = (prcp = prcp_vec, temp = temp_vec, lday = dayl_vec)
        input_mat = Matrix(reduce(hcat, collect(input_ntp[[:prcp, :temp, :lday]]))')
        return input_ntp, input_mat, df
        
    else
        error("Unknown dataset: $dataset")
    end
end

"""
    create_multinode_input(input_mat::Matrix, num_nodes::Int)

Convert single-node input matrix to multi-node format.

# Arguments
- `input_mat`: Single-node input matrix (variables × time)
- `num_nodes`: Number of nodes

# Returns
- Multi-node input array (variables × nodes × time)
"""
function create_multinode_input(input_mat::Matrix, num_nodes::Int)
    return repeat(
        reshape(input_mat, size(input_mat)[1], 1, size(input_mat)[2]),
        1, num_nodes, 1
    )
end

"""
    create_multinode_params(params::ComponentVector, num_nodes::Int)

Replicate parameters for multiple nodes.

# Arguments
- `params`: Single-node parameters (ComponentVector with :params field)
- `num_nodes`: Number of nodes

# Returns
- Multi-node parameters (ComponentVector with expanded params)
"""
function create_multinode_params(params::ComponentVector, num_nodes::Int)
    param_nt = NamedTuple(params[:params])
    node_params = NamedTuple{keys(param_nt)}(
        fill(v, num_nodes) for (k, v) in pairs(param_nt)
    )
    return ComponentVector(params = node_params)
end

"""
    create_multinode_states(states::ComponentVector, num_nodes::Int)

Replicate initial states for multiple nodes.

# Arguments
- `states`: Single-node initial states
- `num_nodes`: Number of nodes

# Returns
- Multi-node initial states (ComponentVector with expanded states)
"""
function create_multinode_states(states::ComponentVector, num_nodes::Int)
    state_nt = NamedTuple(states)
    node_states = NamedTuple{keys(state_nt)}(
        fill(v, num_nodes) for (k, v) in pairs(state_nt)
    )
    return ComponentVector(node_states)
end

"""
    create_test_config(;
        solver = MutableSolver,
        interpolator = Val(DirectInterpolation),
        timeidx = Int[],
        min_value = 1e-6,
        parallel = false
    )

Create a HydroConfig for testing.
"""
function create_test_config(;
    solver = MutableSolver,
    interpolator = Val(DirectInterpolation),
    timeidx = Int[],
    min_value = 1e-6,
    parallel = false
)
    return HydroConfig(
        solver = solver,
        interpolator = interpolator,
        timeidx = timeidx,
        min_value = min_value,
        parallel = parallel
    )
end

"""
    test_component_interface(component)

Test that a component has the expected interface methods.
"""
function test_component_interface(component)
    @test hasmethod(HydroModels.get_input_names, (typeof(component),))
    @test hasmethod(HydroModels.get_output_names, (typeof(component),))
    @test hasmethod(HydroModels.get_param_names, (typeof(component),))
end

"""
Standard ExpHydro parameters for testing.
"""
const EXPHYDRO_PARAMS = (
    Df = 2.674548848,
    Tmax = 0.175739196,
    Tmin = -2.092959084,
    Smax = 1709.46,
    f = 0.0167,
    Qmax = 18.47
)

"""
Standard ExpHydro initial states for testing.
"""
const EXPHYDRO_STATES = (
    snowpack = 0.0,
    soilwater = 1303.0
)

"""
Standard GR4J parameters for testing.
"""
const GR4J_PARAMS = (
    x1 = 320.11,
    x2 = 2.42,
    x3 = 69.63,
    x4 = 1.39
)

"""
Standard GR4J initial states for testing.
"""
const GR4J_STATES = (
    soilwater = 235.97,
    routingstore = 45.47
)

