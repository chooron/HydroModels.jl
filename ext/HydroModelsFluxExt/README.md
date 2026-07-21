# HydroModelsFluxExt

This optional extension activates after `using Flux`. It provides forward
integration for stateless `Flux.Dense` and `Flux.Chain` networks and does not
make Flux a HydroModels core dependency.

Flux parameters come from the model's current `Flux.destructure` vector and
are stored under the HydroModels neural-network parameter tree. The extension
rebuilds a runtime model for each forward call so supplied calibration
parameters are always honored. This is intentionally correct-first and can
allocate in ODE RHS loops; see `benchmark/nn_backend_forward_benchmark.jl`.
