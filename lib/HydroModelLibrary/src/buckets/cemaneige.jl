@reexport module Cemaneige

using ..HydroModels
ifelse_func(x) = ifelse(x > 0, 1, 0)

@variables snowfall rainfall prcp mean_temp max_temp min_temp prcp_ mean_temp_ max_temp_ min_temp_
@variables melt solid_frac infiltration d_thermal snowwater new_snowwater thermal new_thermal
@parameters height altitude zthresh snwthresh CTG Kf

funcs = [
    HydroFlux([prcp] => [prcp_], [zthresh, altitude, height],
        exprs=[ifelse(zthresh > altitude, prcp * exp((altitude - height) * 0.0004), prcp * ifelse(height <= zthresh, exp(zthresh - height) * 0.0004, 1))]
    ),
    HydroFlux([mean_temp, min_temp, max_temp] => [mean_temp_, min_temp_, max_temp_], [altitude, height],
        exprs=[(altitude - height) * (-0.0065) + mean_temp, (altitude - height) * (-0.0065) + min_temp, (altitude - height) * (-0.0065) + max_temp]
    ),
    HydroFlux([prcp, mean_temp, max_temp, min_temp] => [solid_frac], [altitude, zthresh],
        exprs=[ifelse_func(zthresh - altitude) * (ifelse_func(max_temp) * ifelse_func(-min_temp) * (1.0 - max_temp / (max_temp - min_temp)) + ifelse_func(-max_temp)) +
               ifelse_func(altitude - zthresh) * (ifelse_func(-mean_temp) + ifelse_func(3 - mean_temp) * ifelse_func(mean_temp) * (1 - (mean_temp + 1) / 4.0))]
    ),
    HydroFlux([prcp, solid_frac] => [snowfall, rainfall], exprs=[prcp * solid_frac, prcp * (1 - solid_frac)]),
    HydroFlux([thermal, mean_temp] => [new_thermal], [CTG], exprs=[min(0.0, CTG * thermal + (1 - CTG) * mean_temp)]),
    HydroFlux([snowfall, snowwater, new_thermal, mean_temp] => [melt], [Kf, CTG, snwthresh],
        exprs=[ifelse((new_thermal == 0) & (mean_temp > 0), min(Kf * mean_temp, snowwater + snowfall), 0.0) * (0.9 * min(1.0, (snowwater + snowfall) / snwthresh) + 0.1)]
    ),
    HydroFlux([snowfall, melt, snowwater] => [new_snowwater], exprs=[snowwater + snowfall - melt]),
    HydroFlux([rainfall, melt] => [infiltration], exprs=[rainfall + melt]),
]

dfuncs = [StateFlux(new_snowwater => snowwater), StateFlux(new_thermal => thermal)]

cemaneige_bukcet = HydroBucket(name=:cemaneige, funcs=funcs, dfuncs=dfuncs,)
export cemaneige_bukcet
end