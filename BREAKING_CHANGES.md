# Breaking Changes in v0.6.0

## Interpolation Method Names

| Old Name | New Name |
|----------|----------|
| `DirectInterpolation` | `ConstantInterpolation` |
| `EnzymeInterpolation` | `LinearInterpolation` |
| `EnzymeCompatibleInterpolation` | `LinearInterpolation` |

## Multi-Node Parameter Name

`hru_types` → `htypes` (in component constructors)

**Action required:** Update `hru_types` to `htypes` in your code.

## Migration

See [Migration Guide](docs/src/migration_guide_v06.md) for detailed instructions.
