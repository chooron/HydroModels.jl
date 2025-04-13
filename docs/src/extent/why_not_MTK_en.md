# HydroModels.jl and ModelingToolkit.jl

## Why Not Use [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) Directly

The design philosophy of HydroModels.jl is actually completely aligned with ModelingToolkit.jl. Both frameworks use Symbolics.jl for symbolic computation. ModelingToolkit.jl is a more general symbolic system that supports the construction of various problems, and [ModelingToolkitStandardLibrary.jl](https://github.com/SciML/ModelingToolkitStandardLibrary.jl) provides support for basic components in specific domains.

However, I encountered several issues when using ModelingToolkit.jl:

- Although ModelingToolkit.jl can express the ordinary differential equations of hydrological models through symbolic programming, its support for calculations like unit hydrograph and routing processes in hydrological models is not particularly ideal
- Within the same module, such as evaporation calculation formulas and soil calculation modules, the modules needed for different regions may differ due to different formulas, requiring separate construction of calculation formulas to support model building
- Hydrological models typically include multiple ordinary differential equations. Using a single ODESystem becomes relatively chaotic, while using multiple ODESystems becomes relatively complex
- Additionally, ODESystem may not directly support multi-node input and spatial routing process calculations
- Although ModelingToolkitNeuralNets.jl exists for neural network model embedding, the integration is not as straightforward as imagined
- I also encountered issues with automatic differentiation support based on ModelingToolkit.jl, especially with Zygote.jl

Therefore, I decided to build my own model library, referencing the symbolic programming model construction approach, to better support the modeling needs of hydrological models.

## Future Work

It's undeniable that ModelingToolkit.jl is an excellent framework. After meeting the requirements of hydrological models, I will try to make the modules in HydroModels.jl become similar modules in ModelingToolkitStandardLibrary.jl, and inherit the support of ModelingToolkit.jl's AbstractSystem class to provide consistent functional support, thereby maintaining compatibility with more computational capabilities in the SciML ecosystem.
