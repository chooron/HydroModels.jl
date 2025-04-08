# standard hydro bucket
struct SnowpackBucket{F} <: AbstractHydroBucket where {F<:AbstractHydroFlux}
    snowfall_fluxes::AbstractVector{F}
    refreeze_fluxes::AbstractVector{F}
    melt_fluxes::AbstractVector{F}

    function SnowpackBucket{F}(;
        snowfall_fluxes::AbstractVector{F},
        refreeze_fluxes::AbstractVector{F},
        melt_fluxes::AbstractVector{F}
    ) where {F<:AbstractHydroFlux}
        return new{F}(snowfall_fluxes, refreeze_fluxes, melt_fluxes)
    end
end

struct SoilwaterBucket{F} <: AbstractHydroBucket where {F<:AbstractHydroFlux}
    infiltration_fluxes::AbstractVector{F}
    evaporation_fluxes::AbstractVector{F}
end

Struct 