# ----------------------------------------------------------------------
# EISMINT1MovingBenchmark — EISMINT-1 moving-margin experiment.
#
# Reference:
#   Huybrechts, P., Payne, A., et al. (1996).
#   "The EISMINT benchmarks for testing ice-sheet models",
#   Annals of Glaciology 23, 1-12.
#
# Geometry (constant forcing):
#   - Domain: 1500 × 1500 km nominal. Cell-centred grid with
#     `Nx = Ny = floor(L/dx) + 1`. With `dx_km = 50` (EISMINT-1
#     default): `Nx = Ny = 31`, extent 1550 km, summit at i=j=16
#     → (775, 775) km.
#   - Bed: flat (`z_bed = 0` everywhere).
#   - IC: zero ice (`H_ice = 0`); the dome ramps up from radial smb.
#
# Forcing (constant in time):
#   - `T_srf = 270 - 0.01·H_ice` [K] (cosmetic when therm is disabled).
#   - `smb = min(smb_max, s·(R_el − dist_km))` [m/yr] with
#     `R_el = 450 km`, `s = 0.01 m/yr/km`, `smb_max = 0.5 m/yr`.
#   - `Q_geo = 42 mW/m²`.
#
# This is the model-agnostic spec — it returns a NamedTuple of fields
# from `state(b, t)` and writes the same fields to NetCDF via
# `write_fixture!`. The Yelmo-specific YelmoModel constructor lives
# in the `YelmoBenchmarks` package extension.
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

EISMINT-1 moving-margin benchmark [huybrechts_eismint_1996](@citep).

# Fields
 - `xc`, `yc` : cell-centre x/y coordinates (m)
 - `L_km`     : domain extent (km). Default 1500 km.
 - `dx_km`    : cell size (km). Default 50 km.
 - `x_summit`, `y_summit` : summit coordinates (m)
 - `R_el_km`  : radial extent of positive smb (km). Default 450 km.
 - `smb_max`  : maximum smb (m/yr). Default 0.5 m/yr.
 - `smb_grad` : radial smb gradient (m/yr/km). Default 0.01 m/yr/km.
 - `T_srf_const` : surface temperature (K). Default 270 K.
 - `Q_geo_const` : geothermal flux (mW/m²). Default 42 mW/m².
 - `n_glen` : Glen exponent. Default 3.
 - `A_glen` : Glen flow factor (Pa⁻³ yr⁻¹). Default 1e-16 Pa⁻³ yr⁻¹.
"""
struct EISMINT1MovingBenchmark <: AbstractBenchmark
    xc::Vector{Float64}      # cell-centre x [m]
    yc::Vector{Float64}      # cell-centre y [m]
    L_km::Float64
    dx_km::Float64
    x_summit::Float64        # [m]
    y_summit::Float64        # [m]
    R_el_km::Float64
    smb_max::Float64
    smb_grad::Float64        # [m/yr / km]
    T_srf_const::Float64
    Q_geo_const::Float64
    n_glen::Float64
    A_glen::Float64
end

function EISMINT1MovingBenchmark(; dx_km::Real        = 50.0,
                                   L_km::Real         = 1500.0,
                                   R_el_km::Real      = 450.0,
                                   smb_max::Real      = 0.5,
                                   smb_grad::Real     = 0.01,
                                   T_srf_const::Real  = 270.0,
                                   Q_geo_const::Real  = 42.0,
                                   n_glen::Real       = 3.0,
                                   A_glen::Real       = 1e-16)
    dx_f = Float64(dx_km); L_f = Float64(L_km)
    Nx = Int(floor(L_f / dx_f)) + 1
    Ny = Nx

    dx_m = dx_f * 1e3
    xc = collect(range(0.5 * dx_m, (Nx - 0.5) * dx_m; length=Nx))
    yc = collect(range(0.5 * dx_m, (Ny - 0.5) * dx_m; length=Ny))

    x_summit = 0.5 * (xc[1] + xc[end])
    y_summit = 0.5 * (yc[1] + yc[end])

    return EISMINT1MovingBenchmark(xc, yc, L_f, dx_f,
                                    x_summit, y_summit,
                                    Float64(R_el_km),
                                    Float64(smb_max),
                                    Float64(smb_grad),
                                    Float64(T_srf_const),
                                    Float64(Q_geo_const),
                                    Float64(n_glen),
                                    Float64(A_glen))
end

"""
$(TYPEDSIGNATURES)

Surface mass balance at point `(x, y)` in metres:

```julia
    smb = min(smb_max, smb_grad · (R_el_km − dist_km))
```

where `dist_km` is the radial distance from the summit in kilometres.
"""
@inline function eismint_moving_smb(b::EISMINT1MovingBenchmark, x_m::Real, y_m::Real)
    dist_km = sqrt((x_m - b.x_summit)^2 + (y_m - b.y_summit)^2) / 1e3
    return min(b.smb_max, b.smb_grad * (b.R_el_km - dist_km))
end

"""
$(TYPEDSIGNATURES)

Analytical zero-ice IC at `t = 0`. Returns a NamedTuple of model-agnostic
fields keyed by `:xc, :yc, :H_ice, :z_bed, :z_sl, :smb_ref, :T_srf,
:Q_geo, :bmb_shlf, :T_shlf, :H_sed, :ice_allowed, :calv_mask`.

Only `t = 0` is supported here; non-zero times require a model-driven
trajectory and are out of scope for this spec package.
"""
function state(b::EISMINT1MovingBenchmark, t::Real)
    Float64(t) == 0.0 || error(
        "EISMINT1MovingBenchmark.state: only t = 0 is supported " *
        "(got t = $t). Run a forward simulation to obtain non-zero times.")
    return _eismint_moving_analytical_state(b)
end

function _eismint_moving_analytical_state(b::EISMINT1MovingBenchmark)
    Nx, Ny = length(b.xc), length(b.yc)

    H_ice = zeros(Nx, Ny)
    z_bed = zeros(Nx, Ny)
    # Sea-level dropped to -1000 m so every cell is grounded with
    # H_grnd > 0 from the start (cf. mass_balance.jl positive-smb gate
    # on (f_grnd == 0 AND H_ice == 0) cells).
    z_sl  = fill(-1000.0, Nx, Ny)

    smb_ref = [eismint_moving_smb(b, b.xc[i], b.yc[j]) for i in 1:Nx, j in 1:Ny]

    T_srf  = fill(b.T_srf_const, Nx, Ny)
    Q_geo  = fill(b.Q_geo_const, Nx, Ny)
    T_shlf = fill(b.T_srf_const, Nx, Ny)

    bmb_shlf  = zeros(Nx, Ny)
    H_sed     = zeros(Nx, Ny)
    # Yelmo `mask_ice` convention (0 = zero, 1 = fixed, 2 = dynamic).
    # All cells dynamic; EISMINT border handling is via the
    # `boundaries` string, not the ice mask.
    mask_ice  = fill(2.0, Nx, Ny)
    calv_mask = zeros(Nx, Ny)

    return (xc = b.xc, yc = b.yc,
            H_ice = H_ice, z_bed = z_bed, z_sl = z_sl,
            smb_ref = smb_ref, T_srf = T_srf, Q_geo = Q_geo,
            bmb_shlf = bmb_shlf, T_shlf = T_shlf, H_sed = H_sed,
            mask_ice = mask_ice, calv_mask = calv_mask)
end

# Default zeta axes for the analytical fixture writer (31 aa layers).
const _EISMINT_MOVING_NZ_AA = 31
function _eismint_moving_zeta_ac()
    nz_aa = _EISMINT_MOVING_NZ_AA
    zeta_aa = collect(range(0.5/nz_aa, 1.0 - 0.5/nz_aa; length=nz_aa))
    zeta_ac = vcat(0.0, 0.5*(zeta_aa[1:end-1].+zeta_aa[2:end]), 1.0)
    return zeta_aa, zeta_ac
end

"""
$(TYPEDSIGNATURES)

Write the analytical zero-ice IC at `t = 0` to `path` as a NetCDF
restart. Multi-time fixtures (`t > 0`) require a forward simulation
and are not supported here.
"""
function write_fixture!(b::EISMINT1MovingBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    length(times) == 1 ||
        error("write_fixture!(EISMINT1MovingBenchmark, …): pass `times` " *
              "with exactly one entry per call (got $(length(times))).")
    t = Float64(first(times))
    t == 0.0 ||
        error("EISMINT1MovingBenchmark.write_fixture!: only t = 0 supported " *
              "(got t = $t).")
    return _write_eismint_moving_analytical_fixture!(b, path, t)
end

# Analytical fixture writer (t = 0). Factored out so a host whose
# `write_fixture!` override handles both analytical (t = 0) and
# host-driven (t > 0) branches can delegate the analytical path here
# without an `invoke` dance.
function _write_eismint_moving_analytical_fixture!(b::EISMINT1MovingBenchmark,
                                                     path::AbstractString,
                                                     t::Float64)
    s = _eismint_moving_analytical_state(b)
    _, zeta_ac = _eismint_moving_zeta_ac()

    vars = (
        ("H_ice",     s.H_ice,     "m",       "Ice thickness (zero IC)"),
        ("z_bed",     s.z_bed,     "m",       "Bedrock elevation (flat)"),
        ("z_sl",      s.z_sl,      "m",       "Sea level"),
        ("smb_ref",   s.smb_ref,   "m/yr",    "Surface mass balance (radial moving-margin pattern)"),
        ("T_srf",     s.T_srf,     "K",       "Surface temperature"),
        ("Q_geo",     s.Q_geo,     "mW m^-2", "Geothermal flux"),
        ("bmb_shlf",  s.bmb_shlf,  "m/yr",    "Shelf bmb (zero)"),
        ("T_shlf",    s.T_shlf,    "K",       "Shelf base temperature"),
        ("H_sed",     s.H_sed,     "m",       "Sediment thickness (zero)"),
        ("mask_ice",  s.mask_ice,  "1",       "Ice mask (0=none, 1=fixed, 2=dynamic)"),
        ("calv_mask", s.calv_mask, "1",       "Calving mask (zero)"),
    )
    attrs = (
        "benchmark"       => "EISMINT1-moving",
        "solution_type"   => "analytical-IC",
        "time_yr"         => t,
        "L_km"            => b.L_km,
        "dx_km"           => b.dx_km,
        "R_el_km"         => b.R_el_km,
        "smb_max"         => b.smb_max,
        "smb_grad"        => b.smb_grad,
        "T_srf_K"         => b.T_srf_const,
        "Q_geo_mWm2"      => b.Q_geo_const,
        "n_glen"          => b.n_glen,
        "A_glen_Pa-3yr-1" => b.A_glen,
    )
    return _write_restart!(path, b.xc, b.yc,
                           zeta_ac, _DEFAULT_ZETA_ROCK_AC,
                           vars, attrs)
end
