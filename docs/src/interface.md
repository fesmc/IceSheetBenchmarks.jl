# The `AbstractBenchmark` interface

Every benchmark in this package is a concrete subtype of
[`AbstractBenchmark`](@ref). The contract is intentionally small —
three methods plus one optional method — so that hosts can plug in
new benchmarks without touching the core package.

## Methods every benchmark implements

### `state(b, t) -> NamedTuple`

Returns the analytical state at time `t` (years) as a `NamedTuple`
keyed by schema names. The keys are a subset of:

| Key          | Shape        | Units        | Meaning                              |
| ------------ | ------------ | ------------ | ------------------------------------ |
| `xc`, `yc`   | `Vector`     | m            | Cell-centre grid axes                |
| `H_ice`      | `(Nx, Ny)`   | m            | Ice thickness                        |
| `z_bed`      | `(Nx, Ny)`   | m            | Bedrock elevation                    |
| `z_sl`       | `(Nx, Ny)`   | m            | Sea level                            |
| `smb_ref`    | `(Nx, Ny)`   | m/yr (i.e.q.) | Surface mass balance                |
| `T_srf`      | `(Nx, Ny)`   | K            | Surface temperature                  |
| `Q_geo`      | `(Nx, Ny)`   | mW m⁻²       | Geothermal heat flux                 |
| `bmb_shlf`   | `(Nx, Ny)`   | m/yr         | Basal mass balance under shelves     |
| `T_shlf`     | `(Nx, Ny)`   | K            | Sub-shelf ocean temperature          |
| `H_sed`      | `(Nx, Ny)`   | m            | Sediment thickness                   |
| `mask_ice`   | `(Nx, Ny)`   | —            | Ice mask (e.g. 0/1/2)                |
| `calv_mask`  | `(Nx, Ny)`   | —            | Calving mask                         |
| `lsf`        | `(Nx, Ny)`   | —            | Level-set function (ocean/ice)       |
| `ice_allowed`| `(Nx, Ny)`   | —            | Domain mask for allowed ice          |

Only the keys appropriate to the benchmark are populated. Benchmarks
without a closed-form solution may return only `(xc, yc, …)` at `t = 0`
and error for `t > 0` — see each benchmark page.

### `write_fixture!(b, path; times)`

Serialise `state(b, t)` to a NetCDF restart file at `path`. The
default convention is a single time. Returns a `Vector{String}` of
the paths written.

### `analytical_velocity(b, t)` (optional)

Closed-form depth-averaged ice-velocity field `(ux_bar, uy_bar)`,
returned on the face-staggered shapes `(Nx+1, Ny)` and `(Nx, Ny+1)`.
Only implemented when an analytical solution exists (currently the
Bueler / Halfar dome). The fallback throws an informative error.

## Calving-law hooks

Two model-agnostic skeletons are declared in the core module for the
CalvingMIP rate laws:

```julia
calvmip_exp1!(cr_x, cr_y, u_bar, v_bar, H_ice, f_ice, lsf, time;
              xc, yc, r_lim = 750e3)
calvmip_exp2!(cr_x, cr_y, u_bar, v_bar, H_ice, f_ice, lsf, time;
              xc, yc)
```

Hosts can extend these with array-only overloads or, via a package
extension, with overloads keyed on their own field types.

See the [API reference](api.md) for the full docstrings of these
symbols.
