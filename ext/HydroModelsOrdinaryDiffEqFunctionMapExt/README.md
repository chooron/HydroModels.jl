# HydroModelsOrdinaryDiffEqFunctionMapExt

`DiscreteSolver` support through `FunctionMap` is activated only after the
caller loads `OrdinaryDiffEqFunctionMap`:

```julia
using OrdinaryDiffEqFunctionMap
using HydroModels
```

Loading `OrdinaryDiffEq` or `DifferentialEquations` alone does not activate
this extension. The `OrdinaryDiffEqFunctionMap` package is intentionally an
optional dependency; it is not promoted to a HydroModels core dependency.
