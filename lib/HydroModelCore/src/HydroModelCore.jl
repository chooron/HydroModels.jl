module HydroModelCore

using ModelingToolkit: Num
using Symbolics: tosymbol, unwrap
## Abstract Component Types
abstract type AbstractComponent end

abstract type AbstractFlux <: AbstractComponent end
abstract type AbstractHydroFlux <: AbstractFlux end
abstract type AbstractNeuralFlux <: AbstractHydroFlux end
abstract type AbstractStateFlux <: AbstractFlux end

abstract type AbstractElement <: AbstractComponent end
abstract type AbstractBucket <: AbstractElement end
abstract type AbstractHydrograph <: AbstractElement end
abstract type AbstractRoute <: AbstractElement end
abstract type AbstractHydroRoute <: AbstractRoute end
abstract type AbstractModel <: AbstractComponent end

abstract type AbstractNNLayer <: AbstractComponent end
abstract type AbstractNNModel <: AbstractComponent end

include("attributes.jl")
include("display.jl")

end # module HydroModelCore
