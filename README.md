# IceSheetBenchmarks.jl

Portable, model-agnostic Julia package providing the standard
ice-sheet benchmark specifications: geometry, forcing, analytical or
empirical reference solutions. Benchmarks are returned as plain Julia
structs and as NamedTuples of fields keyed by ice-sheet schema names,
so any ice-sheet model implementation can consume them.

## Interface

Every benchmark `b <: AbstractBenchmark` implements:

- `state(b, t)` — analytical state at time `t` as a NamedTuple keyed
  by schema names (`:xc`, `:yc`, `:H_ice`, `:z_bed`, `:smb_ref`,
  `:T_srf`, `:Q_geo`, `:z_sl`, `:mask_ice`, …).
- `write_fixture!(b, path; times)` — serialise the analytical state
  to a NetCDF restart at `path`. Single-time only.
- `analytical_velocity(b, t)` (optional) — closed-form depth-averaged
  velocity `(ux_bar, uy_bar)`, when one exists.

Per-benchmark math helpers (`bueler_test_BC!`, `eismint_moving_smb`,
`calvmip_bed_circular`, `_trough_f17_zbed`, `_hom_c_beta`, …) are
exposed for hosts that need the formulae directly.

## Included benchmarks

| Benchmark                       | Reference                            | Analytical state | Analytical velocity |
| ------------------------------- | ------------------------------------ | :--------------: | :-----------------: |
| `BuelerBenchmark(:B / :C)`      | Bueler et al. 2005 — Halfar dome     | yes              | yes                 |
| `HOMCBenchmark(:C)`             | Pattyn et al. 2008 — ISMIP-HOM Exp C | IC only          | no                  |
| `TroughBenchmark(:F17)`         | Feldmann & Levermann 2017            | host-provided    | no                  |
| `EISMINT1MovingBenchmark`       | Huybrechts et al. 1996               | IC only          | no                  |
| `MISMIP3DBenchmark(:Stnd)`      | Pattyn et al. 2013                   | IC only          | no                  |
| `CalvingMIPBenchmark(:exp1/2)`  | Cornford et al. (CalvingMIP)         | IC only          | no                  |
| `InitMIPBenchmark`              | Goelzer et al. 2018 — InitMIP        | stub             | no                  |

For non-analytical benchmarks the host is expected to provide
`state(b, t > 0)` (e.g. by reading a host-produced reference fixture)
and `write_fixture!(b, …)` (e.g. by driving the host's own model).

## Calving-law helpers

Pure-array implementations of the CalvingMIP velocity-equilibrium /
oscillating-front calving rates are provided as

```julia
calvmip_exp1!(cr_x, cr_y, u_bar, v_bar, H_ice, f_ice, lsf, time;
              xc, yc, r_lim=750e3)
calvmip_exp2!(cr_x, cr_y, u_bar, v_bar, H_ice, f_ice, lsf, time;
              xc, yc)
```

Arguments are plain `Array`s in the natural staggered shapes
(`cr_x` / `u_bar` on x-faces `(Nx+1, Ny)`, `cr_y` / `v_bar` on
y-faces `(Nx, Ny+1)`, the rest on cell centres `(Nx, Ny)`). Hosts
with their own Field types wrap these in their `interior(...)` /
halo-fill conventions.

## Yelmo glue

A package extension `YelmoBenchmarks` activates when both
`IceSheetBenchmarks` and `Yelmo` are loaded. It provides:

- `Yelmo.YelmoModel(b::AbstractBenchmark, t; …)` — build a YelmoModel
  in memory directly from `state(b, t)` with no NetCDF round-trip.
- Field-typed wrappers around `calvmip_exp1!` / `calvmip_exp2!` that
  call `fill_halo_regions!` then delegate to the core array methods.

## Status

Pre-1.0 — interface may evolve. Not yet registered in the General
registry.
