using HydroModels
using HydroModelTools
using DataInterpolations
using ComponentArrays  # For parameter handling
using CSV, DataFrames
using Plots

# Define the step function used in hydrological processes
step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Define variables and parameters (the macro @variables and @parameters are from ModelingToolkit.jl)
@variables temp lday prcp pet snowfall rainfall snowpack melt
@variables soilwater evap baseflow surfaceflow flow
@parameters Tmin Tmax Df Smax Qmax f

# Define the snow component
snow_bucket = @hydrobucket :snow begin
    fluxes = begin
        @hydroflux pet ~ 29.8 * lday * 24 * 0.611 * exp((17.3 * temp) / (temp + 237.3)) / (temp + 273.2)
        @hydroflux snowfall ~ step_func(Tmin - temp) * prcp
        @hydroflux rainfall ~ step_func(temp - Tmin) * prcp
        @hydroflux melt ~ step_func(temp - Tmax) * step_func(snowpack) * min(snowpack, Df * (temp - Tmax))
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

# Define the soil water component
soil_bucket = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux evap ~ step_func(soilwater) * pet * min(1.0, soilwater / Smax)
        @hydroflux baseflow ~ step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))
        @hydroflux surfaceflow ~ max(0.0, soilwater - Smax)
        @hydroflux flow ~ baseflow + surfaceflow
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + melt) - (evap + flow)
    end
end

# Combine components into the complete ExpHydro model
exphydro_model = @hydromodel :exphydro begin
    snow_bucket
    soil_bucket
end

# Load forcing data
file_path = "data/exphydro/01013500.csv"
df = DataFrame(CSV.File(file_path))

# Select a time period for simulation
ts = collect(1:10000)  # One year of daily data

# Extract required input variables
input_data = (
    lday = df[ts, "dayl(day)"],   # Day length
    temp = df[ts, "tmean(C)"],    # Mean temperature
    prcp = df[ts, "prcp(mm/day)"] # Precipitation
)

q_data = df[ts, "flow(mm)"]

# Convert to matrix format required by HydroModels.jl
# We need to make the input matrix in the sort of the input names of the model
input_matrix = reduce(hcat, input_data[HydroModels.get_input_names(exphydro_model)]) |> permutedims

# Define parameter values (calibrated for basin 01013500)
# Create parameter ComponentVector
params = ComponentVector(
    f = 0.0167 , # Baseflow recession parameter
    Smax = 1709.46 , # Maximum soil water storage (mm)
    Qmax = 18.47 , # Maximum baseflow rate (mm/day)
    Df = 2.674 , # Degree-day factor for snowmelt (mm/°C/day)
    Tmax = 0.17 , # Temperature threshold for snowmelt (°C)
    Tmin = -2.09 # Temperature threshold for snow/rain partitioning (°C)
)

# Set initial states
init_states = ComponentVector(
    snowpack = 0.0,       # Initial snow water equivalent (mm)
    soilwater = 1303.00   # Initial soil moisture (mm)
)

# Combine into parameter array structure
pas = ComponentVector(params = params)

# Define simulation configuration
config = (
    # Time indices for simulation
    timeidx = ts,
    # Solver settings
    solver = HydroModels.ManualSolver(mutable = true),
    # Interpolation method for inputs
    interp = LinearInterpolation
)

# Execute the model simulation
results = exphydro_model(
    input_matrix,   # Input data matrix
    pas,            # Parameters
    initstates = init_states,  # Initial states
    config = config, # Configuration
)

# Extract results into named variables
output_names = vcat(
    HydroModels.get_state_names(exphydro_model),
    HydroModels.get_output_names(exphydro_model)
)
output_data = NamedTuple{Tuple(output_names)}(eachslice(results, dims=1)) |> DataFrame

using BenchmarkTools

# Benchmark model execution
@btime exphydro_model(
    $input_matrix,
    $pas,
    initstates = $init_states,
    config = $config
)

# Plot the results
using Plots

# Create a multi-panel plot
p1 = plot(ts, output_data.snowpack, label="Snowpack (mm)", ylabel="SWE (mm)")
p2 = plot(ts, output_data.soilwater, label="Soil Moisture (mm)", ylabel="SM (mm)")
p3 = plot(ts, output_data.flow, label="Simulated Flow", ylabel="Flow (mm/day)")

# Add observed flow if available
if "flow(mm)" in names(df)
    plot!(p3, ts, q_data, label="Observed Flow")
end

# Combine plots
plot(p1, p2, p3, layout=(3,1), size=(800, 600), dpi=800)

savefig("E:\\JlCode\\HydroModels\\docs\\src\\image\\get_start_en\\exphydro_results.png")

observed_flow = df[ts, "flow(mm)"]
simulated_flow = output_data.flow

# Nash-Sutcliffe Efficiency
nse = 1 - sum((observed_flow - simulated_flow).^2) / sum((observed_flow .- mean(observed_flow)).^2)

# Percent Bias
pbias = 100 * (sum(simulated_flow) - sum(observed_flow)) / sum(observed_flow)

println("Nash-Sutcliffe Efficiency: ", nse)
println("Percent Bias: ", pbias, "%")