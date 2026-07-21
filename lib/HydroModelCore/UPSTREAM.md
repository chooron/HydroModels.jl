# Vendored HydroModelCore

`lib/HydroModelCore` is the authoritative copy used by HydroModels.  The
standalone `E:\JlCode\HydroModelCore` checkout is a development mirror and can
be refreshed with:

```text
julia scripts/sync_hydromodelcore.jl --export
```

The vendored source originated from the HydroModelCore repository at commit
`cfcf519c670b18355677a76fd26e4c0a0b60b146` and remains MIT licensed.  Changes to
the standalone checkout should be reviewed and then imported into `lib` as a
deliberate change; the synchronization script only exports the vendored copy.
