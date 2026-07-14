# ----------------------------------------------------------------------
# MISMIP3DBenchmark — Marine Ice-Sheet Model Intercomparison Project
# Phase 3D, experiment "Stnd" (steady-state buildup phase).
#
# Reference:
#   Pattyn, F., et al. (2013). "Grounding-line migration in plan-view
#   marine ice-sheet models", J. Glaciol. 59(215), 410-422.
#
# Geometry (Stnd):
#   - Domain: x ∈ [0, 800] km (Bounded), y ∈ [-50, +50] km (Periodic).
#   - Bed:    z_bed(x, y) = -100 - x_km   (m, y-invariant; slope -1/1000).
#   - IC:     H_ice = 10 m where z_bed ≥ -500 m, else 0
#             (i.e. ice-bearing region is x_km ≤ 400).
#
# Boundaries (Stnd, all constant in time):
#   - smb_ref = 0.5 m/yr
#   - T_srf   = 273.15 K
#   - Q_geo   = 42 mW/m²
#   - calv_mask(nx, :) = .TRUE. (only the eastern column allows calving).
#     Approximated here via `ice_allowed[Nx, :] = 0`.
#
# This is the model-agnostic spec. The literal Fortran 10 m all-floating
# IC is preserved here for fidelity. A consumer (e.g. Yelmo.jl run.jl)
# may swap in the thicker grounded variant `H_ice = max(0, 1000 - 0.9·z_bed)`
# documented in `mismip3D.f90:62-64` to keep its SSA solver well-posed
# from step 1 — that's a model concern, not a benchmark concern.
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

MISMIP3D Stnd benchmark spec. Carries domain axes, bed geometry, IC
slab thickness, surface forcing, and Glen-flow / friction parameters.

`dx_km = 16.0` gives Nx = 51, Ny = 7 — odd Ny yields a centerline cell
at j = 4 for stronger y-symmetry tests. `dx_km = 20.0` (Fortran default)
gives Nx = 41, Ny = 6 (even Ny, no centerline).

Only `variant = :Stnd` is supported.
"""
struct MISMIP3DBenchmark <: AbstractBenchmark
    variant::Symbol
    xc::Vector{Float64}      # cell-centre x [m]
    yc::Vector{Float64}      # cell-centre y [m]
    xmax_km::Float64
    Ly_km::Float64
    dx_km::Float64
    H0::Float64
    z_bed_floor::Float64
    bed_intercept::Float64
    bed_slope::Float64
    smb_const::Float64
    T_srf_const::Float64
    Q_geo_const::Float64
    n_glen::Float64
    A_glen::Float64
    cf_ref::Float64
    N_eff_const::Float64
end

function MISMIP3DBenchmark(variant::Symbol = :Stnd;
                            dx_km::Real         = 16.0,
                            xmax_km::Real       = 800.0,
                            Ly_km::Real         = 100.0,
                            H0::Real            = 10.0,
                            z_bed_floor::Real   = -500.0,
                            bed_intercept::Real = -100.0,
                            bed_slope::Real     = 1.0,
                            smb_const::Real     = 0.5,
                            T_srf_const::Real   = 273.15,
                            Q_geo_const::Real   = 42.0,
                            n_glen::Real        = 3.0,
                            A_glen::Real        = 3.1536e-18,
                            cf_ref::Real        = 3.165176e4,
                            N_eff_const::Real   = 1.0)
    variant === :Stnd || error(
        "MISMIP3DBenchmark: only variant :Stnd is supported (got $(variant)).")

    dx_km_f   = Float64(dx_km)
    xmax_km_f = Float64(xmax_km)
    Ly_km_f   = Float64(Ly_km)

    # Nx / Ny follow the Fortran node-count convention `int(extent/dx)+1`.
    # With dx_km=16: Nx = 51, Ny = 7 (odd → centerline cell at j=4).
    Nx = Int(floor(xmax_km_f / dx_km_f)) + 1
    Ny = Int(floor(Ly_km_f / dx_km_f)) + 1

    dx_m = dx_km_f * 1e3
    xc = collect(range(0.5 * dx_m, (Nx - 0.5) * dx_m; length=Nx))
    yc = collect(range(-0.5 * Ny * dx_m + 0.5 * dx_m,
                        0.5 * Ny * dx_m - 0.5 * dx_m; length=Ny))

    return MISMIP3DBenchmark(variant, xc, yc,
                              xmax_km_f, Ly_km_f, dx_km_f,
                              Float64(H0),
                              Float64(z_bed_floor),
                              Float64(bed_intercept),
                              Float64(bed_slope),
                              Float64(smb_const),
                              Float64(T_srf_const),
                              Float64(Q_geo_const),
                              Float64(n_glen),
                              Float64(A_glen),
                              Float64(cf_ref),
                              Float64(N_eff_const))
end

"""
$(TYPEDSIGNATURES)

Analytical Stnd IC at `t = 0`. Returns a NamedTuple keyed by
`:xc, :yc, :H_ice, :z_bed, :z_sl, :smb_ref, :T_srf, :Q_geo,
:bmb_shlf, :T_shlf, :H_sed, :ice_allowed`.

Only `t = 0` is supported here; non-zero times require a forward
simulation and are out of scope for this spec package.
"""
function state(b::MISMIP3DBenchmark, t::Real)
    Float64(t) == 0.0 || error(
        "MISMIP3DBenchmark.state: only t = 0 is supported " *
        "(got t = $t). Run a forward simulation to obtain non-zero times.")
    return _mismip3d_analytical_state(b)
end

function _mismip3d_analytical_state(b::MISMIP3DBenchmark)
    Nx, Ny = length(b.xc), length(b.yc)

    # z_bed = bed_intercept - bed_slope·x_km (y-invariant).
    z_bed = [b.bed_intercept - b.bed_slope * (b.xc[i] / 1e3) for i in 1:Nx, j in 1:Ny]

    # H_ice = H0 where z_bed ≥ z_bed_floor, else 0 (literal Fortran IC).
    H_ice = [z_bed[i, j] >= b.z_bed_floor ? b.H0 : 0.0 for i in 1:Nx, j in 1:Ny]

    # Yelmo `mask_ice` convention (0 = zero, 1 = fixed, 2 = dynamic).
    # All cells dynamic except the eastern (calving-boundary) column,
    # which is forced to zero ice.
    mask_ice = fill(2.0, Nx, Ny)
    mask_ice[Nx, :] .= 0.0

    smb_ref  = fill(b.smb_const,   Nx, Ny)
    T_srf    = fill(b.T_srf_const, Nx, Ny)
    Q_geo    = fill(b.Q_geo_const, Nx, Ny)
    z_sl     = zeros(Nx, Ny)
    bmb_shlf = zeros(Nx, Ny)
    T_shlf   = fill(b.T_srf_const, Nx, Ny)
    H_sed    = zeros(Nx, Ny)

    return (xc = b.xc, yc = b.yc,
            H_ice = H_ice, z_bed = z_bed, z_sl = z_sl,
            smb_ref = smb_ref, T_srf = T_srf, Q_geo = Q_geo,
            bmb_shlf = bmb_shlf, T_shlf = T_shlf, H_sed = H_sed,
            mask_ice = mask_ice)
end

"""
$(TYPEDSIGNATURES)

Write the analytical Stnd IC at `t = 0` to `path` as a NetCDF restart.
Only `t = 0` is supported.
"""
function write_fixture!(b::MISMIP3DBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    length(times) == 1 ||
        error("write_fixture!(MISMIP3DBenchmark, …): pass `times` " *
              "with exactly one entry per call (got $(length(times))).")
    t = Float64(first(times))
    t == 0.0 ||
        error("MISMIP3DBenchmark.write_fixture!: only t = 0 is supported " *
              "(got t = $t).")
    return _write_mismip3d_analytical_fixture!(b, path, t)
end

# Analytical fixture writer (t = 0). Factored out so a host whose
# `write_fixture!` override handles both analytical (t = 0) and
# host-driven (t > 0) branches can delegate the analytical path here.
function _write_mismip3d_analytical_fixture!(b::MISMIP3DBenchmark,
                                              path::AbstractString,
                                              t::Float64)
    s = _mismip3d_analytical_state(b)

    vars = (
        ("H_ice",    s.H_ice,    "m",       "Ice thickness (Stnd 10 m slab)"),
        ("z_bed",    s.z_bed,    "m",       "Bedrock elevation (-100 - x_km)"),
        ("z_sl",     s.z_sl,     "m",       "Sea level"),
        ("smb_ref",  s.smb_ref,  "m/yr",    "Surface mass balance (Stnd: 0.5 m/yr)"),
        ("T_srf",    s.T_srf,    "K",       "Surface temperature"),
        ("Q_geo",    s.Q_geo,    "mW m^-2", "Geothermal flux"),
        ("bmb_shlf", s.bmb_shlf, "m/yr",    "Shelf bmb (zero)"),
        ("T_shlf",   s.T_shlf,   "K",       "Shelf base temperature"),
        ("H_sed",    s.H_sed,    "m",       "Sediment thickness (zero)"),
        ("mask_ice", s.mask_ice, "1",       "Ice mask (0=none, 1=fixed, 2=dynamic; eastern column = 0)"),
    )
    attrs = (
        "benchmark"       => "MISMIP3D-$(string(b.variant))",
        "solution_type"   => "analytical-IC",
        "time_yr"         => t,
        "xmax_km"         => b.xmax_km,
        "Ly_km"           => b.Ly_km,
        "dx_km"           => b.dx_km,
        "H0_m"            => b.H0,
        "bed_intercept"   => b.bed_intercept,
        "bed_slope"       => b.bed_slope,
        "z_bed_floor_m"   => b.z_bed_floor,
        "smb_const_m_yr"  => b.smb_const,
        "T_srf_K"         => b.T_srf_const,
        "Q_geo_mWm2"      => b.Q_geo_const,
        "n_glen"          => b.n_glen,
        "A_glen_Pa-3yr-1" => b.A_glen,
        "cf_ref"          => b.cf_ref,
        "N_eff_Pa"        => b.N_eff_const,
    )
    return _write_restart!(path, b.xc, b.yc,
                           _DEFAULT_ZETA_AC, _DEFAULT_ZETA_ROCK_AC,
                           vars, attrs)
end
