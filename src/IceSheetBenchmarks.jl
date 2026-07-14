module IceSheetBenchmarks

using DocStringExtensions
using NCDatasets

# ----------------------------------------------------------------------
# AbstractBenchmark — interface contract.
#
# Concrete subtypes carry the parameters of a benchmark (grid axes,
# forcing parameters, Glen-flow parameters, …) and must implement:
#
#   - `state(b, t)` — analytical state at time `t`, as a NamedTuple
#     keyed by ice-sheet field names (`:xc`, `:yc`, `:H_ice`, `:z_bed`,
#     `:smb_ref`, `:T_srf`, `:Q_geo`, `:z_sl`, ...).
#   - `write_fixture!(b, path; times)` — serialize the analytical state
#     at one or more times to a NetCDF restart file.
#   - `analytical_velocity(b, t)` (optional) — closed-form
#     depth-averaged ice-velocity field `(ux_bar, uy_bar)`, when one
#     exists (e.g. Halfar dome). Benchmarks without an analytical
#     velocity solution leave this unimplemented.
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Supertype of every benchmark spec in this package. Concrete subtypes
carry the parameters of a benchmark (grid axes, forcing parameters,
Glen-flow parameters, …) and implement the small interface
documented under [`state`](@ref), [`write_fixture!`](@ref), and
optionally [`analytical_velocity`](@ref).

# Available subtypes:
 - [`HalfarDomeBenchmark`](@ref) — Halfar dome (Bueler et al. 2005).
 - [`HOMCBenchmark`](@ref) — HOM-C (ISMIP-HOM, Hindmarsh et al. 2008).
 - [`TroughBenchmark`](@ref) — Feldmann-Levermann (2017) trough.
 - [`EISMINT1MovingBenchmark`](@ref) — EISMINT1 moving-margin (Payne et al. 2000).
 - [`MISMIP3DBenchmark`](@ref) — MISMIP3D Stnd (Pattyn et al. 2013).
 - [`CalvingMIPBenchmark`](@ref) — CalvingMIP (Asay-Davis et al. 2021).
 - [`InitMIPBenchmark`](@ref) — InitMIP (Asay-Davis et al. 2021).
"""
abstract type AbstractBenchmark end

"""
$(TYPEDSIGNATURES)

Analytical state of benchmark `b` at time `t` (years). Returns a
`NamedTuple` keyed by ice-sheet schema names — see the
[Interface](@ref "The `AbstractBenchmark` interface") page for the
catalogue of keys. Benchmarks without a closed-form solution may
restrict to `t = 0` (initial condition) and error otherwise.
"""
function state end

"""
$(TYPEDSIGNATURES)

Serialise `state(b, t)` at one or more times to a NetCDF restart at
`path`. Returns the vector of file paths written. Some benchmarks
support only a single time.
"""
function write_fixture! end

"""
$(TYPEDSIGNATURES)

Closed-form depth-averaged ice-velocity field. `ux_bar` is shape
`(Nx+1, Ny)` and `uy_bar` is `(Nx, Ny+1)` (face-staggered). Only
implemented for benchmarks with a known analytical solution
(currently [`HalfarDomeBenchmark`](@ref)). The fallback throws an
informative error.
"""
function analytical_velocity end

analytical_velocity(b::AbstractBenchmark, t::Real) = error(
    "analytical_velocity not implemented for $(typeof(b)). " *
    "Use a concrete benchmark subtype with a closed-form velocity solution.")

"""
$(TYPEDSIGNATURES)

In-place velocity-equilibrium calving-rate law for CalvingMIP Exp1/Exp3.
Skeleton declared here so hosts (or the `YelmoBenchmarks` extension)
can extend it. See the [CalvingMIP page](@ref "CalvingMIP") for the
arguments and staggered shapes.
"""
function calvmip_exp1! end

"""
$(TYPEDSIGNATURES)

In-place oscillating-front calving-rate law for CalvingMIP Exp2/Exp4.
Skeleton declared here so hosts can extend it.
"""
function calvmip_exp2! end

export AbstractBenchmark
export state, write_fixture!, analytical_velocity
export calvmip_exp1!, calvmip_exp2!

include("fixture_io.jl")
include("bueler.jl")
include("hom_c.jl")
include("trough.jl")
include("eismint_moving.jl")
include("mismip3d.jl")
include("calvingmip.jl")
include("initmip.jl")

export HalfarDomeBenchmark, bueler_gamma, bueler_test_BC!
export HOMCBenchmark
export TroughBenchmark
export EISMINT1MovingBenchmark, eismint_moving_smb
export MISMIP3DBenchmark
export CalvingMIPBenchmark, calvmip_bed_circular, calvmip_bed_thule
export InitMIPBenchmark

end # module
