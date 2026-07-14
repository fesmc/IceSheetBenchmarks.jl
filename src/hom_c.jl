# ----------------------------------------------------------------------
# HOMCBenchmark — ISMIP-HOM Experiment C ("ice stream over a sloping bed
# with spatially-varying basal friction").
#
# Reference: Pattyn et al. (2008), "Benchmark experiments for higher-
# order and full-Stokes ice sheet models (ISMIP-HOM)", The Cryosphere
# 2, 95-108. Yelmo Fortran implementation: `yelmo/tests/yelmo_ismiphom.f90`
# (the `EXPC` branch, lines 131-146) and `yelmo/par/yelmo_ISMIPHOM.nml`.
#
# Geometry (Yelmo Fortran convention):
#
#   alpha = 0.1° = 0.1 π / 180 rad
#   omega = 2π / L                 # β-perturbation wavenumber
#
#   z_srf = -x · tan(alpha)
#   z_bed =  z_srf - 1000 m         (uniform 1000 m thick slab over a
#                                    sloping bed)
#   H_ice = 1000 m
#
# Basal friction (Pa·yr·m^-1):
#
#   β(x, y) = β₀ + β_amp · sin(omega · x) · sin(omega · y)
#
# with β₀ = 1000 and β_amp = 1000. Note: this is the **Yelmo Fortran**
# convention (perturbation amplitude 1.0). The published Pattyn 2008
# formula uses amplitude 0.9 → β oscillates over [100, 1900]; the Yelmo
# Fortran version oscillates over [0, 2000]. We mirror the Fortran
# reference (amplitude exposed as a struct field for future flexibility).
#
# Boundary conditions: fully periodic in x and y. The Fortran reference
# pads the domain with `f_extend = 0.5` (half a period on each side) to
# minimise edge effects under clamped BC. Yelmo.jl uses fully-periodic
# BC directly so f_extend is unnecessary.
#
# Validation strategy: HOM-C has no closed-form velocity solution; the
# Pattyn 2008 paper publishes inter-model centre-line velocities for
# the HOM intercomparison. The Yelmo.jl regression target is **180°
# rotational anti-symmetry** of the SSA solution — under fully-periodic
# BC the (x, y) → (L-x, L-y) rotation flips the slope-driven driving
# stress sign while preserving β / H / ATT / forcing, so the velocity
# must satisfy ux(x, y) = -ux(L-x, L-y) and uy(x, y) = -uy(L-x, L-y).
#
# `analytical_velocity` is intentionally NOT implemented (falls through
# to the default error stub).
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

ISMIP-HOM Experiment C benchmark struct. Carries domain axes, slope,
material parameters, and basal-friction perturbation parameters.

Only `variant = :C` is supported in this milestone.

`dx_km` defaults to `0.25 · (L_km / 10)` (the Yelmo Fortran convention,
`yelmo_ismiphom.f90:60`). For the canonical L = 80 km case this gives
`dx_km = 2.0` and Nx = Ny = 40.

The benchmark deliberately uses an unpadded `[0, L_km] × [0, L_km]`
domain with fully-periodic BC; the Fortran `f_extend = 0.5` padding is
unnecessary under periodic BC.
"""
struct HOMCBenchmark <: AbstractBenchmark
    variant::Symbol
    xc::Vector{Float64}
    yc::Vector{Float64}
    L_km::Float64
    dx_km::Float64
    H::Float64
    alpha_rad::Float64
    A_glen::Float64
    n_glen::Float64
    beta0::Float64
    beta_amp::Float64
    rho_ice::Float64
    g::Float64
end

function HOMCBenchmark(variant::Symbol = :C;
                       L_km::Real      = 80.0,
                       dx_km::Real     = 0.25 * (Float64(L_km) / 10.0),
                       alpha_deg::Real = 0.1,
                       H::Real         = 1000.0,
                       A_glen::Real    = 1e-16,
                       n_glen::Real    = 3.0,
                       beta0::Real     = 1000.0,
                       beta_amp::Real  = 1000.0,
                       rho_ice::Real   = 910.0,
                       g::Real         = 9.81)
    variant === :C || error(
        "HOMCBenchmark: only variant :C is supported (got $(variant)).")

    L_km_f  = Float64(L_km)
    dx_km_f = Float64(dx_km)
    Nx_f = L_km_f / dx_km_f
    isinteger(Nx_f) || error(
        "HOMCBenchmark: L_km/dx_km must be integer (got $Nx_f).")
    Nx = Int(Nx_f)

    # Cell-centre axes on [0, L] × [0, L]. Under fully-periodic BC the
    # cell at i = 1 has centre dx/2 and the cell at i = Nx has centre
    # L - dx/2 (the "missing" cell at L is the periodic copy of cell 1).
    dx_m = dx_km_f * 1e3
    xc = collect(range(0.5 * dx_m, (Nx - 0.5) * dx_m; length=Nx))
    yc = copy(xc)

    return HOMCBenchmark(variant,
                         xc, yc,
                         L_km_f, dx_km_f,
                         Float64(H),
                         Float64(alpha_deg) * π / 180.0,
                         Float64(A_glen),
                         Float64(n_glen),
                         Float64(beta0),
                         Float64(beta_amp),
                         Float64(rho_ice),
                         Float64(g))
end

"""
$(TYPEDSIGNATURES)

Analytical IC at time `t`. Returns a NamedTuple with:
  - `xc`, `yc`     — grid axes (metres).
  - `H_ice`        — 1000 m uniform.
  - `z_bed`        — `-x · tan α - H` (linear in x).
  - `z_sl`         — −10 000 m (matches `yelmo_ismiphom.f90:180`,
                     forcing all ice to be grounded so the bed slope
                     drives the surface slope and the driving stress
                     is non-trivial).
  - `smb_ref`      — zero (HOM-C has no mass-balance forcing).
  - `T_srf`        — 263.15 K (arbitrary; HOM-C is isothermal).
  - `Q_geo`        — 50.0 mW/m² (arbitrary).

The basal-friction perturbation β is **not** returned in `state` — it
lives on `dyn.beta` / `dyn.beta_acx` / `dyn.beta_acy` (which are not
part of the schema-routed Center-aligned state). Use
`_setup_hom_c_beta!` after constructing the YelmoModel to fill the
β fields from the analytical formula.
"""
function state(b::HOMCBenchmark, t::Real)
    Nx, Ny = length(b.xc), length(b.yc)
    H_ice = fill(b.H, Nx, Ny)
    z_srf = [-b.xc[i] * tan(b.alpha_rad) for i in 1:Nx, j in 1:Ny]
    z_bed = z_srf .- b.H
    z_sl  = fill(-10_000.0, Nx, Ny)   # force grounded ice (Fortran convention)
    smb   = zeros(Nx, Ny)
    Tsrf  = fill(263.15, Nx, Ny)
    Qgeo  = fill(50.0,   Nx, Ny)
    return (xc = b.xc, yc = b.yc,
            H_ice = H_ice, z_bed = z_bed, z_sl = z_sl,
            smb_ref = smb, T_srf = Tsrf, Q_geo = Qgeo)
end

"""
$(TYPEDSIGNATURES)

Serialise the analytical HOM-C IC at time `t = first(times)` to a
NetCDF restart at `path`. Single-time only (HOM-C has no time
evolution at the IC level — β is steady, geometry is steady). Uses the
default uniform 11-point ice / 5-point rock sigma axes.

Returns a 1-element `Vector{String}` containing `path`.
"""
function write_fixture!(b::HOMCBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    length(times) == 1 ||
        error("write_fixture!(HOMCBenchmark, …): multi-time fixtures " *
              "not supported (got $(length(times)) times).")
    t = Float64(first(times))
    s = state(b, t)

    vars = (
        ("H_ice",   s.H_ice,   "m",       "Ice thickness (HOM-C uniform slab)"),
        ("z_bed",   s.z_bed,   "m",       "Bedrock elevation (sloping bed, alpha=0.1°)"),
        ("z_sl",    s.z_sl,    "m",       "Sea level (very negative — keeps ice grounded)"),
        ("smb_ref", s.smb_ref, "m/yr",    "Surface mass balance (zero for HOM-C)"),
        ("T_srf",   s.T_srf,   "K",       "Surface temperature (HOM-C is isothermal)"),
        ("Q_geo",   s.Q_geo,   "mW m^-2", "Geothermal flux (arbitrary; isothermal test)"),
    )
    attrs = (
        "benchmark"       => "ISMIPHOM-$(string(b.variant))",
        "solution_type"   => "analytical-IC",
        "time_yr"         => t,
        "L_km"            => b.L_km,
        "dx_km"           => b.dx_km,
        "alpha_deg"       => b.alpha_rad * 180.0 / π,
        "H_m"             => b.H,
        "A_glen_Pa-3yr-1" => b.A_glen,
        "n_glen"          => b.n_glen,
        "beta0"           => b.beta0,
        "beta_amp"        => b.beta_amp,
    )
    return _write_restart!(path, b.xc, b.yc,
                           _DEFAULT_ZETA_AC, _DEFAULT_ZETA_ROCK_AC,
                           vars, attrs)
end

# Closed-form β at an (x_m, y_m) point (metres):
#   β = β₀ + β_amp · sin(2π x / L) · sin(2π y / L)
#
# Used by the Yelmo-side `_setup_hom_c_beta!` to fill dyn.beta fields
# on a constructed YelmoModel. Pure math, exposed here so any host
# can reuse the formula.
@inline function _hom_c_beta(b::HOMCBenchmark, x_m::Real, y_m::Real)
    omega = 2π / (b.L_km * 1e3)
    return b.beta0 + b.beta_amp * sin(omega * Float64(x_m)) * sin(omega * Float64(y_m))
end
