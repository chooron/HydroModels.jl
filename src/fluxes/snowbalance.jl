module SnowBalance
using ..HydroModels

"""
    Simple melt 
"""
function SNOBAL_SIMPLE_MELT(;
    snowstorage::Number=first(@variables snowstorage),
    snowmelt::Number=first(@variables snowmelt),
    cumulmelt::Number=first(@variables cumulmelt),
    overflow::Number=first(@variables overflow),
    potential_melt::Number=first(@variables potential_melt),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    temp::Number=first(@variables temp),
    ma::Number=first(@variables ma),
    mamax::Number=first(@parameters mamax),
    mamin::Number=first(@parameters mamin),
    malpha::Number=first(@parameters malpha),
    tbm::Number=first(@parameters tbm),
    name::Union{Nothing,Symbol}=nothing
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux ma ~ min(mamax, mamin * (1 + malpha * cumulmelt))
            @hydroflux potential_melt ~ ma * max(0.0, temp - tbm)
            @hydroflux snowmelt ~ min(snowstorage, potential_melt)
            @hydroflux overflow ~ snowmelt + rainfall
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt
            @stateflux cumulmelt ~ ifelse(snowstorage > 0, snowmelt, -cumulmelt)
        end
    end
end

"""
    HBV snow balance
"""
function SNOBAL_HBV(;
    snowstorage::Number=first(@variables snowstorage),
    liquidstorage::Number=first(@variables liquidstorage),
    temp::Number=first(@variables temp [description = "the temperature"]),
    potential_melt::Number=first(@variables potential_melt),
    cumulmelt::Number=first(@variables cumulmelt),
    refreeze::Number=first(@variables refreeze),
    overflow::Number=first(@variables overflow),
    snowmelt::Number=first(@variables snowmelt),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    ma::Number=first(@variables ma),
    ka::Number=first(@parameters ka [description = "the refreeze coefficient"]),
    mamax::Number=first(@parameters mamax),
    mamin::Number=first(@parameters mamin),
    malpha::Number=first(@parameters malpha),
    tbm::Number=first(@parameters tbm),
    tbf::Number=first(@parameters tbf),
    const_swi::Number=first(@parameters const_swi [description = "the maximum possible liquid water storage of the snowpack"]),
    name::Union{Nothing,Symbol}=nothing
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux ma ~ min(mamax, mamin * (1 + malpha * cumulmelt))
            @hydroflux potential_melt ~ ma * max(0.0, temp - tbm)
            @hydroflux snowmelt ~ min(snowstorage, potential_melt)
            @hydroflux refreeze ~ min(max(0.0, liquidstorage), ka * max(tbf - temp, 0.0))
            @hydroflux overflow ~ max(0.0, liquidstorage + rainfall + snowmelt - const_swi * snowstorage)
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt + refreeze
            @stateflux liquidstorage ~ snowmelt + rainfall - refreeze - overflow
            @stateflux cumulmelt ~ ifelse(snowstorage > 0, snowmelt, -cumulmelt)
        end
    end
end

"""
    HMETS snow balance
"""
function SNOBAL_HMETS(;
    snowstorage::Number=first(@variables snowstorage),
    liquidstorage::Number=first(@variables liquidstorage),
    cumulmelt::Number=first(@variables cumulmelt),
    tmean::Number=first(@variables tmean),
    potential_melt::Number=first(@variables potential_melt),
    swi::Number=first(@variables swi),
    refreeze::Number=first(@variables refreeze),
    snowmelt::Number=first(@variables snowmelt),
    snowfall::Number=first(@variables snowfall),
    rainfall::Number=first(@variables rainfall),
    overflow::Number=first(@variables overflow),
    ma::Number=first(@variables ma),
    alpha::Number=first(@parameters alpha),
    swimax::Number=first(@parameters swimax),
    swimin::Number=first(@parameters swimin),
    kf::Number=first(@parameters kf),
    tbm::Number=first(@parameters tbm),
    tbf::Number=first(@parameters tbf),
    f::Number=first(@parameters f),
    mamax::Number=first(@parameters mamax),
    mamin::Number=first(@parameters mamin),
    malpha::Number=first(@parameters malpha),
    name::Union{Nothing,Symbol}=nothing
)
    return HydroModels.@hydrobucket name begin
        fluxes = begin
            @hydroflux ma ~ min(mamax, mamin * (1 + malpha * cumulmelt))
            @hydroflux potential_melt ~ ma * max(0.0, tmean - tbm)
            @hydroflux snowmelt ~ min(snowstorage, potential_melt)
            @hydroflux refreeze ~ min(max(0.0, liquidstorage), kf * max(0.0, tbf - tmean)^f)
            @hydroflux swi ~ max(swimin, swimax * (1 - alpha * cumulmelt))
            @hydroflux overflow ~ max(0.0, liquidstorage + rainfall + snowmelt - swi * snowstorage)
        end
        dfluxes = begin
            @stateflux snowstorage ~ snowfall - snowmelt + refreeze
            @stateflux liquidstorage ~ snowmelt + rainfall - refreeze - overflow
            @stateflux cumulmelt ~ ifelse(snowstorage > 0, snowmelt, -cumulmelt)
        end
    end
end

end