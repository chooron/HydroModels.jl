module HydroModelsMooncakeExt

using HydroModels
using Mooncake

"""
HydroModels model descriptions are static during a differentiated solve.

The numerical values that must receive derivatives (parameters, initial states,
and forcing data) are passed separately to the model and are therefore not
discarded by these structural tangent rules.
"""
Mooncake.tangent_type(::Type{<:HydroModels.AbstractComponent}) = Mooncake.NoTangent
Mooncake.tangent_type(::Type{<:HydroModels.HydroConfig}) = Mooncake.NoTangent

# Components/configurations are static descriptions during differentiation.
# Marking them as raw friendly tangents also prevents Mooncake from recursively
# traversing generated-function closures and metadata containers.
Mooncake.friendly_tangent_cache(::HydroModels.AbstractComponent) =
    Mooncake.FriendlyTangentCache{Mooncake.AsRaw}(nothing)
Mooncake.friendly_tangent_cache(::HydroModels.HydroConfig) =
    Mooncake.FriendlyTangentCache{Mooncake.AsRaw}(nothing)

end # module
