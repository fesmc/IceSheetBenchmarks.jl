# ----------------------------------------------------------------------
# Shared NetCDF restart writer.
#
# Every analytical `write_fixture!` serialises the same restart schema:
# the (xc, yc, zeta, zeta_ac, zeta_rock, zeta_rock_ac) coordinate system
# plus a set of 2D (xc, yc) state variables and global provenance
# attributes. `_write_restart!` owns that boilerplate so the individual
# benchmarks only declare their variable list, zeta axes, and attributes.
# ----------------------------------------------------------------------

# Default sigma axes shared by the benchmarks that use uniform layers
# (Bueler, HOM-C, MISMIP3D). Benchmarks needing a different vertical
# discretisation (EISMINT-1, CalvingMIP) build their own `zeta_ac` and
# reuse `_DEFAULT_ZETA_ROCK_AC` for the bedrock column.
const _DEFAULT_ZETA_AC      = collect(range(0.0, 1.0; length=11))
const _DEFAULT_ZETA_ROCK_AC = collect(range(0.0, 1.0; length=5))

# Common floating-point element type for a benchmark built from mixed
# numeric inputs. `float` guarantees an `AbstractFloat` result (so Int
# kwargs promote to Float64), while genuine Float32 / Dual inputs are
# preserved. Used by the smart constructors to pick each benchmark's
# `T` parameter.
@inline _bench_eltype(xs::Real...) = float(promote_type(map(typeof, xs)...))

# _write_restart!(path, xc, yc, zeta_ac, zeta_rock_ac, vars, attrs) -> [path]
#
# Write a NetCDF restart at `path` with the standard benchmark schema.
#
#   - xc, yc                : cell-centre axes in metres (stored as km).
#   - zeta_ac, zeta_rock_ac : cell-edge sigma axes for ice and bedrock;
#                             cell-centre zeta / zeta_rock are the edge
#                             midpoints.
#   - vars  : iterable of (name, data, units, long_name) for the 2D
#             (xc, yc) state variables, written in the given order.
#   - attrs : iterable of name => value global attributes.
#
# Returns a 1-element Vector{String} containing `path`, matching the
# write_fixture! contract.
function _write_restart!(path::AbstractString,
                         xc::AbstractVector{<:Real},
                         yc::AbstractVector{<:Real},
                         zeta_ac::AbstractVector{<:Real},
                         zeta_rock_ac::AbstractVector{<:Real},
                         vars,
                         attrs)
    mkpath(dirname(path))
    isfile(path) && rm(path)

    Nx = length(xc); Ny = length(yc)
    Nz_ac      = length(zeta_ac)
    Nz_rock_ac = length(zeta_rock_ac)

    NCDataset(path, "c") do ds
        defDim(ds, "xc",           Nx)
        defDim(ds, "yc",           Ny)
        defDim(ds, "zeta",         Nz_ac - 1)
        defDim(ds, "zeta_ac",      Nz_ac)
        defDim(ds, "zeta_rock",    Nz_rock_ac - 1)
        defDim(ds, "zeta_rock_ac", Nz_rock_ac)

        xv = defVar(ds, "xc", Float64, ("xc",))
        xv[:] = xc ./ 1e3; xv.attrib["units"] = "km"
        yv = defVar(ds, "yc", Float64, ("yc",))
        yv[:] = yc ./ 1e3; yv.attrib["units"] = "km"

        zc = defVar(ds, "zeta", Float64, ("zeta",))
        zc[:] = 0.5 .* (zeta_ac[1:end-1] .+ zeta_ac[2:end]); zc.attrib["units"] = "1"
        zac = defVar(ds, "zeta_ac", Float64, ("zeta_ac",))
        zac[:] = zeta_ac; zac.attrib["units"] = "1"
        zrc = defVar(ds, "zeta_rock", Float64, ("zeta_rock",))
        zrc[:] = 0.5 .* (zeta_rock_ac[1:end-1] .+ zeta_rock_ac[2:end])
        zrc.attrib["units"] = "1"
        zrac = defVar(ds, "zeta_rock_ac", Float64, ("zeta_rock_ac",))
        zrac[:] = zeta_rock_ac; zrac.attrib["units"] = "1"

        for (name, data, units, longname) in vars
            v = defVar(ds, name, Float64, ("xc", "yc"))
            v[:, :] = data
            v.attrib["units"]     = units
            v.attrib["long_name"] = longname
        end

        for (k, val) in attrs
            ds.attrib[k] = val
        end
    end

    return [path]
end
