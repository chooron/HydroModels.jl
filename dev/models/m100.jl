using ModelingToolkit
using Lux
using StableRNGs

step_func(x) = (tanh(5.0 * x) + 1.0) * 0.5

# Model parameters
# Physical process parameters
@parameters Tmin Tmax Df Smax f Qmax
# Normalization parameters
@parameters snowpack_std snowpack_mean
@parameters soilwater_std soilwater_mean
@parameters prcp_std prcp_mean
@parameters temp_std temp_mean

# Model variables
# Input variables
@variables prcp temp lday
# State variables
@variables snowpack soilwater
# Process variables
@variables pet rainfall snowfall melt
# Neural network variables
@variables log_evap_div_lday log_flow asinh_melt asinh_ps asinh_pr
@variables norm_snw norm_slw norm_temp norm_prcp

# Neural network definition
m100_nn = Lux.Chain(
    Lux.Dense(4, 32, tanh),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 32, leakyrelu),
    Lux.Dense(32, 5),
    name=:m100nn
)

# Neural network flux
m100_nn_flux = NeuralFlux(
    [norm_snw, norm_slw, norm_prcp, norm_temp] =>
        [log_evap_div_lday, log_flow, asinh_melt, asinh_ps, asinh_pr],
    m100_nn
)

# Main model bucket
m100_bucket = @hydrobucket :m100_bucket begin
    fluxes = [
        @hydroflux begin
            norm_snw ~ (snowpack - snowpack_mean) / snowpack_std
            norm_slw ~ (soilwater - soilwater_mean) / soilwater_std
            norm_prcp ~ (prcp - prcp_mean) / prcp_std
            norm_temp ~ (temp - temp_mean) / temp_std
        end,
        @neuralflux [log_evap_div_lday, log_flow, asinh_melt, asinh_ps, asinh_pr] ~ m100_nn([norm_snw, norm_slw, norm_prcp, norm_temp]),
        @hydroflux melt ~ relu(sinh(asinh_melt) * step_func(snowpack))
    ]
    dfluxes = [
        @stateflux snowpack ~ relu(sinh(asinh_ps)) * step_func(-temp) - melt,
        @stateflux soilwater ~ relu(sinh(asinh_pr)) + melt - step_func(soilwater) * lday * exp(log_evap_div_lday) - step_func(soilwater) * exp(log_flow)
    ]
end

# Complete model
m100_model = @hydromodel :m100 begin
    m100_bucket
end

export m100_model, m100_bucket