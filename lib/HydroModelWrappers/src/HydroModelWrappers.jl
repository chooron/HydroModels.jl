module HydroModelWrappers

using Chain

abstract type AbstractHydroWrapper end
abstract type AbstractDataPreprocessor <: AbstractHydroWrapper end
abstract type AbstractDataPostprocessor <: AbstractHydroWrapper end
abstract type AbstractComponentDecorator <: AbstractHydroWrapper end

include("namedtuple_wrapper.jl")
# include("estimate_params.jl")
# include("neural_wrapper.jl")
# include("record_states.jl")
include("stats_outlet.jl")

export NamedTuplePreprocessor, NamedTuplePostprocessor, SelectComponentOutlet

end # module HydroModelWrappers
