# ----------------------------------------------------------------------
# TroughBenchmark — Feldmann-Levermann (2017) trough domain.
#
# Reference: Feldmann, J. & Levermann, A. (2017). "From cyclic ice
# streaming to Heinrich-like events: the grow-and-surge instability
# in the Parallel Ice Sheet Model", The Cryosphere 11, 1913-1932.
#
# Fortran reference: `yelmo/tests/yelmo_trough.f90` (`trough_f17_topo_init`
# at line 357) for the closed-form bed geometry and
# `yelmo/par/yelmo_TROUGH-F17.nml` for the domain-shape parameters
# (lx, ly, fc, dc, wc, x_cf).
#
# The closed-form geometry (`_trough_f17_zbed`) and the spec struct
# live here. The benchmark has **no** closed-form transient solution;
# host-side code is expected to provide `state(b, t)` (e.g. by reading
# a host-produced reference fixture) and `write_fixture!(b, …)` (host-
# specific fixture generation).
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Feldmann-Levermann (2017) "TROUGH-F17" trough benchmark spec.

`variant`:
  - `:F17` — the standard Feldmann-Levermann (2017, TC) setup. The
    only variant supported at present.

`dx_km` is the grid resolution in km. Default `4.0` matches the
namelist value. The grid has `Nx = lx/dx + 1` points spanning
`[0, lx]` in x and `Ny = ly/dx + 1` points spanning `[-ly/2, +ly/2]`
in y, matching the Fortran `yelmo_init_grid` call at
`yelmo_trough.f90:100-102`.

The trough geometry parameters (`fc, dc, wc, x_cf`) and forcing
(`Tsrf, smb, Q_geo`) default to the values in
`yelmo_TROUGH-F17.nml`.

This benchmark has no closed-form transient solution; the
`AbstractBenchmark` `state(b, t)` and `write_fixture!(b, …)`
methods are expected to be implemented by the host (e.g. Yelmo.jl
reads a reference fixture produced by its Fortran backend).
"""
struct TroughBenchmark <: AbstractBenchmark
    variant::Symbol
    xc::Vector{Float64}
    yc::Vector{Float64}
    dx_km::Float64

    # Trough geometry parameters (from yelmo_TROUGH-F17.nml /
    # yelmo_trough.f90:357 trough_f17_topo_init).
    lx_km::Float64
    ly_km::Float64
    fc_km::Float64        # characteristic side-wall width
    dc_m::Float64         # depth of bed trough below side walls
    wc_km::Float64        # half-width of bed trough
    x_cf_km::Float64      # calving-front x-position

    # Forcing (uniform).
    Tsrf_const::Float64   # [degC] surface T
    smb_const::Float64    # [m/yr] surface mass balance
    Qgeo_const::Float64   # [mW/m²] geothermal flux

    # Physical constants (echoed from yelmo_phys_const TROUGH section).
    rho_ice::Float64
    g::Float64
end

function TroughBenchmark(variant::Symbol;
                         dx_km::Real      = 4.0,
                         lx_km::Real      = 700.0,
                         ly_km::Real      = 160.0,
                         fc_km::Real      = 16.0,
                         dc_m::Real       = 500.0,
                         wc_km::Real      = 24.0,
                         x_cf_km::Real    = 640.0,
                         Tsrf_const::Real = -20.0,
                         smb_const::Real  = 0.3,
                         Qgeo_const::Real = 70.0,
                         rho_ice::Real    = 910.0,
                         g::Real          = 9.81)
    variant === :F17 ||
        error("TroughBenchmark: only :F17 variant is implemented (got $variant).")

    # Mirror the Fortran grid construction in yelmo_trough.f90:96-102:
    #   xmax =  lx                             # km
    #   ymax =  ly/2,  ymin = -ly/2            # km
    #   nx = int(xmax/dx)+1                    # cell-centre count
    #   ny = int((ymax-ymin)/dx)+1
    #   x[i] = 0   + (i-1)*dx                  # i = 1..nx, in km
    #   y[j] = ymin + (j-1)*dx                 # j = 1..ny, in km
    nx = Int(floor(lx_km / dx_km)) + 1
    ymax_km =  ly_km / 2.0
    ymin_km = -ly_km / 2.0
    ny = Int(floor((ymax_km - ymin_km) / dx_km)) + 1

    xc_m = collect(range(0.0, (nx - 1) * dx_km * 1e3; length=nx))
    yc_m = collect(range(ymin_km * 1e3, ymax_km * 1e3; length=ny))

    return TroughBenchmark(variant, xc_m, yc_m, Float64(dx_km),
                           Float64(lx_km), Float64(ly_km),
                           Float64(fc_km), Float64(dc_m),
                           Float64(wc_km), Float64(x_cf_km),
                           Float64(Tsrf_const), Float64(smb_const),
                           Float64(Qgeo_const),
                           Float64(rho_ice), Float64(g))
end

# F17 closed-form bedrock elevation at one (x_km, y_km) point.
# Ports `yelmo_trough.f90:382-396`. Pure math; used by host-side IC
# setters to fill `z_bed` from the analytical formula.
@inline function _trough_f17_zbed(x_km::Real, y_km::Real,
                                  fc_km::Real, dc_m::Real, wc_km::Real;
                                  zb_deep::Real = -720.0)
    zb_x = -150.0 - 0.84 * abs(Float64(x_km))                   # [m]
    e1 = -2.0 * (Float64(y_km) - Float64(wc_km)) / Float64(fc_km)
    e2 =  2.0 * (Float64(y_km) + Float64(wc_km)) / Float64(fc_km)
    zb_y = (Float64(dc_m) / (1.0 + exp(e1))) +
           (Float64(dc_m) / (1.0 + exp(e2)))                    # [m]
    return max(zb_x + zb_y, Float64(zb_deep))
end

function _trough_analytical_state(b::TroughBenchmark)
    Nx, Ny = length(b.xc), length(b.yc)
    z_bed = [_trough_f17_zbed(b.xc[i] / 1e3, b.yc[j] / 1e3,
                              b.fc_km, b.dc_m, b.wc_km)
             for i in 1:Nx, j in 1:Ny]
    H_ice   = zeros(Nx, Ny)                          # grows from smb
    z_sl    = zeros(Nx, Ny)
    smb_ref = fill(b.smb_const, Nx, Ny)
    T_srf   = fill(b.Tsrf_const + 273.15, Nx, Ny)    # struct field is °C
    Q_geo   = fill(b.Qgeo_const, Nx, Ny)
    return (xc = b.xc, yc = b.yc,
            H_ice = H_ice, z_bed = z_bed, z_sl = z_sl,
            smb_ref = smb_ref, T_srf = T_srf, Q_geo = Q_geo)
end

"""
$(TYPEDSIGNATURES)

Analytical F17 IC at `t = 0`: the closed-form trough bed
([`_trough_f17_zbed`](@ref)), zero ice, sea level at 0, and the
uniform forcing carried by the struct (`smb_const`, `Tsrf_const`
converted to K, `Qgeo_const`).

The benchmark has no closed-form transient solution, so only `t = 0`
is supported; non-zero times require a forward simulation and should
be provided by the host (e.g. by reading a host-produced reference
fixture).
"""
function state(b::TroughBenchmark, t::Real)
    Float64(t) == 0.0 || error(
        "TroughBenchmark.state: only t = 0 is supported (got t = $t). " *
        "Run a forward simulation for non-zero times.")
    return _trough_analytical_state(b)
end

const _TROUGH_DEFAULT_ZETA_AC      = _DEFAULT_ZETA_AC
const _TROUGH_DEFAULT_ZETA_ROCK_AC = _DEFAULT_ZETA_ROCK_AC

"""
$(TYPEDSIGNATURES)

Write the analytical F17 IC at `t = 0` to `path` as a NetCDF restart.
Only `t = 0` is supported; host-driven trajectories (`t > 0`) require a
forward simulation.
"""
function write_fixture!(b::TroughBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    length(times) == 1 ||
        error("write_fixture!(TroughBenchmark, …): pass `times` with " *
              "exactly one entry per call (got $(length(times))).")
    t = Float64(first(times))
    t == 0.0 ||
        error("TroughBenchmark.write_fixture!: only t = 0 is supported " *
              "(got t = $t).")

    s = _trough_analytical_state(b)
    vars = (
        ("H_ice",   s.H_ice,   "m",       "Ice thickness (zero IC)"),
        ("z_bed",   s.z_bed,   "m",       "Bedrock elevation (F17 trough)"),
        ("z_sl",    s.z_sl,    "m",       "Sea level"),
        ("smb_ref", s.smb_ref, "m/yr",    "Surface mass balance (constant)"),
        ("T_srf",   s.T_srf,   "K",       "Surface temperature (constant)"),
        ("Q_geo",   s.Q_geo,   "mW m^-2", "Geothermal flux (constant)"),
    )
    attrs = (
        "benchmark"     => "TROUGH-$(string(b.variant))",
        "solution_type" => "analytical-IC",
        "time_yr"       => t,
        "dx_km"         => b.dx_km,
        "lx_km"         => b.lx_km,
        "ly_km"         => b.ly_km,
        "fc_km"         => b.fc_km,
        "dc_m"          => b.dc_m,
        "wc_km"         => b.wc_km,
        "x_cf_km"       => b.x_cf_km,
        "smb_const"     => b.smb_const,
        "Tsrf_degC"     => b.Tsrf_const,
        "Qgeo_mWm2"     => b.Qgeo_const,
    )
    return _write_restart!(path, b.xc, b.yc,
                           _TROUGH_DEFAULT_ZETA_AC, _TROUGH_DEFAULT_ZETA_ROCK_AC,
                           vars, attrs)
end
