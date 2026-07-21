# Neural network backends

HydroModels keeps neural-network integrations optional. The existing Lux
backend remains the full-featured backend; Flux and SimpleChains are loaded
only when their packages are loaded.

| Backend | Forward | Explicit parameters | Explicit state | Mooncake | Recommended use |
| --- | :---: | :---: | :---: | :---: | --- |
| Lux | yes | yes | yes | project-level validation | stateful and general networks |
| Flux | yes | `Flux.destructure` | no in this extension | future candidate | stateless Dense/Chain models |
| SimpleChains | yes | native flat vector | no | no | lightweight CPU forward evaluation |

## Flux

Load Flux before constructing a Flux-backed component:

```julia
using HydroModels
using Flux

chain = Flux.Chain(Flux.Dense(2 => 8, tanh), Flux.Dense(8 => 1))
flux = NeuralFlux(inputs, outputs, chain; chain_name=:runoff_net)
```

The extension calls `Flux.destructure` once when the component is created. The
current model parameters are stored under `nns.<chain_name>.params` as a
`ComponentVector`. Each `nn_func(x, p)` call rebuilds a runtime model from the
provided parameter vector, so parameter changes are never hidden by a stale
cache. This is correct for calibration but can allocate in a high-frequency
ODE RHS.

The supported first-stage subset is stateless `Flux.Dense` and `Flux.Chain`
compositions using ordinary pure forward activations such as `identity`,
`relu`, `tanh`, and `sigmoid`. BatchNorm, Dropout, recurrent layers, GPU
execution, and custom mutable layers are not part of this compatibility
promise.

## SimpleChains

```julia
using HydroModels
using SimpleChains

chain = SimpleChain(
    SimpleChains.static(2),
    TurboDense(tanh, 8),
    TurboDense(identity, 1),
)
flux = NeuralFlux(inputs, outputs, chain; chain_name=:runoff_net)
```

SimpleChains is natively structure/parameter separated. HydroModels calls
`SimpleChains.init_params(chain, Float32; rng=...)` during parameter
initialization and stores the resulting flat vector under `params`. The
backend supports the feed-forward `SimpleChain`/`TurboDense` subset used by
the HydroModels factories and is CPU forward-only.

SimpleChains backend is forward-only in HydroModels. Direct Mooncake
differentiation is unsupported because the tested SimpleChains execution path
uses llvmcall-based kernels that Mooncake cannot currently translate. This is
a tested NO-GO for end-to-end Mooncake training, not an untested feature.

Factories select the backend explicitly:

```julia
create_neural_bucket(Val(:flux); ...)
create_simple_neural_bucket(Val(:flux); ...)
create_neural_bucket(Val(:simplechains); ...)
create_simple_neural_bucket(Val(:simplechains); ...)
```

The default no-argument factory continues to select Lux when Lux is loaded.
No Flux or SimpleChains package is loaded by `using HydroModels`.
