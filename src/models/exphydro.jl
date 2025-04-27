module exphydro
using ..HydroModels
using ..HydroModels: step_func
# define variables and parameters
@variables T Ep P snowfall rainfall snowpack melt
@variables soilwater Ep evap baseflow surfaceflow flow rainfall

@parameters f [description="Minimum Terature for snowfall", bounds=(0.0, 0.1), unit="celsius"]
@parameters Smax [description="Maximum Terature for snowmelt", bounds=(100.0, 2000.0), unit="celsius"]
@parameters Qmax [description="Degree-day factor", bounds=(10.0, 50.0), unit="mm/degree celsius/d"]
@parameters Df [description="Maximum soil moisture storage", bounds=(0.0, 5.0), unit="mm"]
@parameters Tmax [description="Maximum runoff", bounds=(0.0, 3.0), unit="mm/d"]
@parameters Tmin [description="Runoff coefficient", bounds=(-3.0, 0.0), unit="d-1"]

# define model components
bucket_1 = @hydrobucket :surface begin
    fluxes = begin
        @hydroflux begin
            snowfall ~ step_func(Tmin - T) * P
            rainfall ~ step_func(T - Tmin) * P
        end
        @hydroflux melt ~ step_func(T - Tmax) * step_func(snowpack) * min(snowpack, Df * (T - Tmax))
    end
    dfluxes = begin
        @stateflux snowpack ~ snowfall - melt
    end
end

bucket_2 = @hydrobucket :soil begin
    fluxes = begin
        @hydroflux evap ~ step_func(soilwater) * Ep * min(1.0, soilwater / Smax)
        @hydroflux baseflow ~ step_func(soilwater) * Qmax * exp(-f * (max(0.0, Smax - soilwater)))
        @hydroflux surfaceflow ~ max(0.0, soilwater - Smax)
        @hydroflux flow ~ baseflow + surfaceflow
    end
    dfluxes = begin
        @stateflux soilwater ~ (rainfall + melt) - (evap + flow)
    end
end

model = @hydromodel :exphydro begin
    bucket_1
    bucket_2
end

end