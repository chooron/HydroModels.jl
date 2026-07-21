"""Check or export the vendored HydroModelCore mirror.

The copy under `lib/HydroModelCore` is authoritative for HydroModels.  The
optional export mode updates the sibling standalone checkout for development;
it never changes the vendored copy.
"""

const root = normpath(joinpath(@__DIR__, ".."))
const vendored = joinpath(root, "lib", "HydroModelCore")
const mirror = normpath(joinpath(root, "..", "HydroModelCore"))
const tracked = [
    "Project.toml",
    "LICENSE",
    joinpath("src", "HydroModelCore.jl"),
    joinpath("src", "build.jl"),
    joinpath("src", "symbol_build.jl"),
    joinpath("test", "runtests.jl"),
    joinpath("test", "test_build.jl"),
    joinpath("test", "test_parameters.jl"),
    joinpath("test", "test_symbol_build.jl"),
]

function same_file(a, b)
    isfile(a) && isfile(b) && read(a) == read(b)
end

function check_mirror()
    missing = String[]
    changed = String[]
    for rel in tracked
        src = joinpath(vendored, rel)
        dst = joinpath(mirror, rel)
        !isfile(dst) ? push!(missing, rel) : !same_file(src, dst) && push!(changed, rel)
    end
    isempty(missing) && isempty(changed) || begin
        !isempty(missing) && println("Missing in sibling mirror: ", join(missing, ", "))
        !isempty(changed) && println("Different from vendored copy: ", join(changed, ", "))
        return false
    end
    println("HydroModelCore mirror is synchronized.")
    true
end

function export_mirror()
    isdir(vendored) || error("Vendored Core does not exist: $vendored")
    for rel in tracked
        src = joinpath(vendored, rel)
        dst = joinpath(mirror, rel)
        isfile(src) || error("Tracked vendored file is missing: $src")
        mkpath(dirname(dst))
        cp(src, dst; force=true)
    end
    println("Exported vendored HydroModelCore to $mirror")
end

if "--export" in ARGS
    export_mirror()
elseif "--check" in ARGS || isempty(ARGS)
    exit(check_mirror() ? 0 : 1)
else
    println("Usage: julia scripts/sync_hydromodelcore.jl [--check|--export]")
    exit(2)
end
