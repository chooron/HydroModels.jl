# Migrating from v0.6 to v0.7

HydroModels v0.7 removes AD-framework-era API constraints and makes the
callable and parameter-tree contracts explicit. This is a breaking release:
update configuration values and any code that relied on the old `Val(Type)`
interpolation convention.

## Summary of breaking changes

| v0.6 | v0.7 | Reason |
| --- | --- | --- |
| `interpolator = Val(ConstantInterpolation)` | `interpolator = ConstantInterpolation` | Types are already dispatch values; wrapping them in `Val` is unnecessary public API. |
| only `Val(Type)` interpolation specifications | a type or a pre-built callable | Reuse an external interpolator without reconstructing it. |
| function-only component constructors | any Julia callable | Functor structs retain concrete type and compose with the ecosystem. |
| extension-local `update_ca` merge | `merge_componentvectors` | A recursive, structural overlay has defined nested-field semantics. |

## Interpolator specifications

Pass a type when HydroModels should construct the interpolation object from the
model input and `timeidx`:

```julia
config = HydroConfig(
    solver=ODESolver,
    interpolator=LinearInterpolation,
    timeidx=1:365,
)
```

Pass an already-constructed callable when its data, time coordinates, or
options are owned by the caller. It must implement `itp(t)` and must be
consistent with the forcing data supplied to the model:

```julia
itp = DataInterpolations.CubicSpline(forcing, times)
config = HydroConfig(solver=ODESolver, interpolator=itp, timeidx=times)
```

`Val(LinearInterpolation)`, `Val(ConstantInterpolation)`, and
`Val(DataInterpolations.LinearInterpolation)` are no longer valid public
specifications. Replace them with the corresponding type.

For ODEs, a built-in `ConstantInterpolation` type or instance causes forcing
knot times to be passed as `tstops`. For an external discontinuous callable,
configure its discontinuities through the SciML solve options; HydroModels
cannot infer them safely.

## Callable component constructors

The functional constructors for `HydroFlux`, `HydroBucket`, `NeuralFlux`,
`HydroRoute`, `RouteIRF`, and `UnitHydrograph` now accept any callable object,
not only values whose type is a subtype of `Function`.

```julia
struct RunoffLaw
    scale::Float64
end

(law::RunoffLaw)(inputs, params) = (law.scale .* params.k .* inputs[1],)

flux = HydroFlux(RunoffLaw(1.0);
    inputs=[:precipitation], outputs=[:runoff], params=[:k])
```

Existing anonymous functions continue to work unchanged.

## ComponentVector overlays and optimisation

`merge_componentvectors(base, overrides...; strict=false)` recursively merges
the structural `NamedTuple` representation and creates one final
`ComponentVector`. A later leaf overrides an earlier leaf, while unrelated
fields remain in the result. `strict=true` rejects fields that do not exist in
the base tree.

```julia
base = ComponentVector(params=(k=1.0, bias=2.0), nns=(net=(weight=3.0,),))
fixed = ComponentVector(params=(bias=5.0,))
full = merge_componentvectors(base, fixed; strict=true)
# (params = (k = 1.0, bias = 5.0), nns = (net = (weight = 3.0,),))
```

The `OptimizationProblem` extension uses this strict overlay for
`fixed_params`. A fixed field must already be present in `initial_params` (or
the component's default parameter tree); misspellings now raise an
`ArgumentError` instead of being silently ignored. The parameter tree must be
fixed for the lifetime of an optimization problem. Do not add, remove, or
rename `ComponentVector` fields inside the differentiated objective.

Mooncake regression coverage verifies that an overridden field receives its
gradient only through the retained value. For a merge of `b` with
`b_override`, the gradient of the discarded `b` is zero and the gradient of
`b_override` is retained.

## Mooncake support baseline

v0.7 tests the ComponentVector overlay path with Mooncake. The supported
baseline is ComponentArrays 0.15.34 or newer and Mooncake 0.5.25 or newer.
The values and the parameter-tree shape have different roles: values may be
active, but field names and nesting are static metadata.

## Other compatibility notes

- `HydroConfig.timeidx` preserves the caller's element type; it is no longer
  coerced to `Vector{Int}`. This permits non-integer numeric time coordinates.
- `DirectInterpolation` remains an alias for `ConstantInterpolation` in this
  release, but new code should use `ConstantInterpolation`.
- The internal solver enum still dispatches with `Val`; this is an
  implementation detail. The public `solver=MutableSolver` style is unchanged.

## Migration checklist

1. Replace `Val(SomeInterpolator)` with `SomeInterpolator`.
2. If you already own the forcing interpolation object, pass that callable as
   `interpolator` and keep its data/time axis consistent with the model input.
3. Replace custom recursive ComponentArray merge code with
   `merge_componentvectors`; use `strict=true` for parameter overrides.
4. Run your gradient tests with a scalar loss and `AutoMooncake()` after
   changing a parameter-tree schema.
