# ----------------------------------------------------------------------
# CalvingMIPBenchmark — circular- and Thule-domain benchmarks from the
# CalvingMIP intercomparison (https://github.com/JRowanJordan/CalvingMIP).
#
# References:
#   - CalvingMIP wiki: https://github.com/JRowanJordan/CalvingMIP/wiki
#   - Fortran driver:  yelmo/tests/yelmo_calving.f90
#   - Fortran geometry: yelmo/tests/calving_benchmarks.f90
#   - Calving laws:    yelmo/src/physics/calving/calving_ac.f90
#
# Domain (circular, exp1/2):
#   x, y ∈ [-800, 800] km  (1600 × 1600 km square)
#   Bed:  z_bed(r) = Bc − (Bc − Bl) · r² / R0²
#         R0 = 800 km, Bc = 900 m, Bl = −2000 m  (parabolic bowl)
#   Land cells (above SL) are pinned by the calving driver to lsf=−1.
#   Initial H_ice = 0 everywhere; ice grows from constant SMB.
#
# Boundary forcing (all experiments, constant in time):
#   smb       = 0.3 m/yr
#   T_srf     = 223.15 K (−50 °C)
#   Q_geo     = 42 mW/m²
#   z_sl      = 0 m
#   bmb_shlf  = 0 m/yr
#
# This file is the model-agnostic spec — `state(b, t)` returns a
# NamedTuple of fields and `write_fixture!` writes the same fields to
# NetCDF. Calving-law hooks (which require ice-model field types) live
# in the `YelmoBenchmarks` package extension.
# ----------------------------------------------------------------------

export CalvingMIPBenchmark
export calvmip_bed_circular, calvmip_bed_thule

# -----------------------------------------------------------------------
# Bed geometry (port of calving_benchmarks.f90)
# -----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

CalvingMIP circular-domain bed elevation (m) at point (x, y) in metres.
Parabolic bowl: z_bed = Bc − (Bc − Bl) · r² / R0².

Parameters match the CalvingMIP spec (calving_benchmarks.f90:63-87):
  R0 = 800 km, Bc = 900 m, Bl = −2000 m.
"""
function calvmip_bed_circular(x::Real, y::Real)
    R0 = 800e3
    Bc = 900.0
    Bl = -2000.0
    r = sqrt(x^2 + y^2)
    return Bc - (Bc - Bl) * r^2 / R0^2
end

"""
$(TYPEDSIGNATURES)

CalvingMIP Thule-domain bed elevation (m) at point (x, y) in metres.
Parabolic bowl with cosine undulations:
  a(r) = Bc − (Bc − Bl) · r² / R0²
  l(θ) = R0 − cos(2θ) · R0/2
  z_bed = Ba · cos(3π r / l(θ)) + a(r)

Parameters: R0 = 800 km, Bc = 900 m, Bl = −2000 m, Ba = 1100 m.
"""
function calvmip_bed_thule(x::Real, y::Real)
    R0 = 800e3
    Bc = 900.0
    Bl = -2000.0
    Ba = 1100.0
    r = sqrt(x^2 + y^2)
    θ = atan(y, x)
    l = R0 - cos(2θ) * R0 / 2
    a = Bc - (Bc - Bl) * r^2 / R0^2
    return Ba * cos(3π * r / l) + a
end

# -----------------------------------------------------------------------
# CalvingMIPBenchmark struct
# -----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

CalvingMIP benchmark on the circular (exp1/2) or Thule (exp3/4/5) domain.

Currently supported:
  - `:exp1` — equilibrium calving pinning the front at r = 750 km (circular).
  - `:exp2` — oscillating calving front, chained from an Exp1 state (circular).
  - `:exp3` — equilibrium calving pinning the front at r = 750 km (Thule).
  - `:exp4` — oscillating calving front, chained from an Exp3 state (Thule).

The calving law itself is supplied by the host model via a hook (e.g.
`YelmoHooks.calv_flt`); this struct only carries the model-agnostic
geometry and forcing.
"""
Base.@kwdef struct CalvingMIPBenchmark{T<:AbstractFloat} <: AbstractBenchmark
    exp::Symbol
    domain::Symbol
    xc::Vector{T}
    yc::Vector{T}
    dx_km::T
    smb_const::T = 0.3
    T_srf_const::T = 223.15
    Q_geo_const::T = 42.0
end

function CalvingMIPBenchmark(exp::Symbol;
    dx_km::Real=25.0,
    smb_const::Real=0.3,
    T_srf_const::Real=223.15,
    Q_geo_const::Real=42.0)
    exp in (:exp1, :exp2, :exp3, :exp4) || error(
        "CalvingMIPBenchmark: unsupported exp = $exp. Supported: :exp1, :exp2, :exp3, :exp4.")

    domain = exp in (:exp1, :exp2) ? :circular : :thule

    T = _bench_eltype(dx_km, smb_const, T_srf_const, Q_geo_const)
    dx_m = T(dx_km) * 1000
    # Domain: x, y ∈ [−800, 800] km  → 1600 km extent → 64 cells at 25 km.
    extent_m = T(1_600_000)
    N = Int(round(extent_m / dx_m))
    axis = collect(range(-extent_m/2 + dx_m/2, extent_m/2 - dx_m/2; length=N))

    return CalvingMIPBenchmark{T}(exp, domain, axis, copy(axis), T(dx_km),
        T(smb_const), T(T_srf_const), T(Q_geo_const))
end

# -----------------------------------------------------------------------
# Analytical t = 0 IC
# -----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Analytical CalvingMIP IC at `t = 0`: H_ice = 0 everywhere, z_bed from
the circular bowl formula, lsf = +1 everywhere.

Non-zero times require a model-driven trajectory (e.g. via a YelmoMirror
fixture or an Exp1 restart) and are out of scope for this spec package.
"""
function state(b::CalvingMIPBenchmark, t::Real)
    Float64(t) == 0.0 || error(
        "CalvingMIPBenchmark.state: only t = 0 is supported (got t = $t). " *
        "Run a forward simulation or load a restart for non-zero times.")
    return _calvingmip_analytical_state(b)
end

function _calvingmip_analytical_state(b::CalvingMIPBenchmark)
    Nx = length(b.xc);
    Ny = length(b.yc)

    bed_fn = b.domain == :thule ? calvmip_bed_thule : calvmip_bed_circular
    z_bed = [bed_fn(b.xc[i], b.yc[j]) for i in 1:Nx, j in 1:Ny]
    H_ice = zeros(Nx, Ny)
    # lsf = +1 (all ocean) — ice will grow from SMB; the calving step's
    # above-SL pin will force lsf = −1 over land each step.
    lsf = ones(Nx, Ny)

    smb_ref = fill(b.smb_const, Nx, Ny)
    T_srf = fill(b.T_srf_const, Nx, Ny)
    Q_geo = fill(b.Q_geo_const, Nx, Ny)
    z_sl = zeros(Nx, Ny)
    bmb_shlf = zeros(Nx, Ny)
    T_shlf = fill(b.T_srf_const, Nx, Ny)
    H_sed = zeros(Nx, Ny)
    ice_allowed = ones(Nx, Ny)
    calv_mask = zeros(Nx, Ny)

    return (xc=b.xc, yc=b.yc,
        H_ice=H_ice, z_bed=z_bed, z_sl=z_sl,
        smb_ref=smb_ref, T_srf=T_srf, Q_geo=Q_geo,
        bmb_shlf=bmb_shlf, T_shlf=T_shlf, H_sed=H_sed,
        ice_allowed=ice_allowed, calv_mask=calv_mask,
        lsf=lsf)
end

# Default zeta axes for the analytical fixture writer.
const _CALVINGMIP_NZ_AA = 11
function _calvingmip_zeta_ac()
    nz_aa = _CALVINGMIP_NZ_AA
    zeta_aa = collect(range(0.5/nz_aa, 1.0 - 0.5/nz_aa; length=nz_aa))
    zeta_ac = vcat(0.0, 0.5*(zeta_aa[1:(end-1)] .+ zeta_aa[2:end]), 1.0)
    return zeta_aa, zeta_ac
end

"""
$(TYPEDSIGNATURES)

Write the analytical zero-ice IC at `t = 0` to `path` as a NetCDF
restart. Multi-time fixtures (`t > 0`) require a forward simulation;
those live next to the host-model integration (e.g. the YelmoMirror
fixture writer in `test/benchmarks/calvingmip.jl`).
"""
function write_fixture!(b::CalvingMIPBenchmark, path::AbstractString;
    times::AbstractVector{<:Real}=[0.0])
    length(times) == 1 ||
        error("write_fixture!(CalvingMIPBenchmark, …): pass `times` " *
              "with exactly one entry per call (got $(length(times))).")
    t = Float64(first(times))
    t == 0.0 ||
        error("CalvingMIPBenchmark.write_fixture!: only t = 0 supported " *
              "(got t = $t). Use a model-side fixture writer for t > 0.")

    s = _calvingmip_analytical_state(b)
    bed_longname = b.domain == :thule ? "Bedrock elevation (Thule bowl + undulations)" :
                   "Bedrock elevation (parabolic bowl)"
    _, zeta_ac = _calvingmip_zeta_ac()

    vars = (
        ("H_ice", s.H_ice, "m", "Ice thickness (zero IC)"),
        ("z_bed", s.z_bed, "m", bed_longname),
        ("z_sl", s.z_sl, "m", "Sea level"),
        ("smb_ref", s.smb_ref, "m/yr", "Surface mass balance (constant)"),
        ("T_srf", s.T_srf, "K", "Surface temperature"),
        ("Q_geo", s.Q_geo, "mW m^-2", "Geothermal flux"),
        ("bmb_shlf", s.bmb_shlf, "m/yr", "Shelf bmb (zero)"),
        ("T_shlf", s.T_shlf, "K", "Shelf base temperature"),
        ("H_sed", s.H_sed, "m", "Sediment thickness (zero)"),
        ("ice_allowed", s.ice_allowed, "1", "Ice-allowed mask"),
        ("calv_mask", s.calv_mask, "1", "Calving mask (zero)"),
        ("lsf", s.lsf, "1", "Level-set function (+1 ocean, −1 ice)"),
    )
    attrs = (
        "benchmark" => "CalvingMIP-$(b.exp)",
        "solution_type" => "analytical-IC",
        "time_yr" => t,
        "dx_km" => b.dx_km,
        "smb_const" => b.smb_const,
        "T_srf_K" => b.T_srf_const,
        "Q_geo_mWm2" => b.Q_geo_const,
    )
    return _write_restart!(path, b.xc, b.yc,
        zeta_ac, _DEFAULT_ZETA_ROCK_AC,
        vars, attrs)
end

# ----------------------------------------------------------------------
# CalvingMIP calving-law math — plain-array core.
#
# Both `calvmip_exp1!` and `calvmip_exp2!` are documented in the
# package's `function … end` declarations in `IceSheetBenchmarks.jl`;
# this file owns the array bodies. Hosts with their own Field /
# halo-fill conventions wrap these — see the YelmoBenchmarks
# extension in `ext/YelmoBenchmarks.jl` for the Yelmo wrapper.
#
# Shape convention (matches the natural staggered grid layout):
#
#   cr_x, u_bar   :: AbstractMatrix  size (Nx+1, Ny)      [x-faces]
#   cr_y, v_bar   :: AbstractMatrix  size (Nx,   Ny+1)    [y-faces]
#   H_ice, f_ice,
#   lsf           :: AbstractMatrix  size (Nx,   Ny)      [aa-centres]
#   xc, yc        :: AbstractVector  length Nx, Ny        [aa-centres, metres]
#
# `H_ice`, `f_ice`, `lsf` are presently unused by the formulae but
# accepted so the array and Field signatures match.
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

CalvingMIP Exp1/3 velocity-equilibrium calving with the front pinned
at radius `r_lim` (default 750 km). Port of `calvmip_exp1` in
`yelmo/src/physics/calving/calving_ac.f90:395-482`.

Algorithm:
  1. `cr = −u` everywhere (velocity-equilibrium calving).
  2. For aa-centres inside `r < r_lim`, zero the calving rate on faces
     whose other neighbour is also inside `r_lim`. Faces straddling the
     `r_lim` boundary keep `cr = −u`. Result: front pinned at the
     `r_lim` circle.
"""
function calvmip_exp1!(cr_x::AbstractMatrix, cr_y::AbstractMatrix,
    u_bar::AbstractMatrix, v_bar::AbstractMatrix,
    H_ice, f_ice, lsf, time::Real;
    xc::AbstractVector{<:Real},
    yc::AbstractVector{<:Real},
    r_lim::Real=750e3)
    Nx = length(xc);
    Ny = length(yc)

    # Step 1: cr = −u (velocity equilibrium).
    @inbounds for j in axes(cr_x, 2), i in axes(cr_x, 1)
        cr_x[i, j] = -u_bar[i, j]
    end
    @inbounds for j in axes(cr_y, 2), i in axes(cr_y, 1)
        cr_y[i, j] = -v_bar[i, j]
    end

    # Step 2: zero faces that lie strictly inside `r_lim`.
    # x-face `i` sits between centres i−1 and i; the face to the RIGHT
    # of centre `i` is `cr_x[i+1, j]`. Same convention for y.
    @inbounds for j in 1:Ny, i in 1:Nx
        r_ij = sqrt(xc[i]^2 + yc[j]^2)
        r_ij >= r_lim && continue

        if i < Nx
            r_right = sqrt(xc[i+1]^2 + yc[j]^2)
            if r_right < r_lim
                cr_x[i+1, j] = 0.0
            end
        end
        if i > 1
            r_left = sqrt(xc[i-1]^2 + yc[j]^2)
            if r_left < r_lim
                cr_x[i, j] = 0.0
            end
        end
        if j < Ny
            r_top = sqrt(xc[i]^2 + yc[j+1]^2)
            if r_top < r_lim
                cr_y[i, j+1] = 0.0
            end
        end
        if j > 1
            r_bot = sqrt(xc[i]^2 + yc[j-1]^2)
            if r_bot < r_lim
                cr_y[i, j] = 0.0
            end
        end
    end
    return cr_x, cr_y
end

"""
$(TYPEDSIGNATURES)

CalvingMIP Exp2 oscillating-front calving rate. Port of
`calvmip_exp2` in `yelmo/src/physics/calving/calving_ac.f90:484-533`.

Net front velocity in the radial direction:
  w  = u + cr = (u/|u|) · wv,
  wv = −300 sin(2π t / 1000)   [m/yr],

so the front oscillates between ±300 m/yr radially with period
1000 yr. The face-local speed `|u|` uses the cross-staggered partner
velocity (4-point average of the orthogonal face values), regularised
with `max(|u|, 1e-8)`. The simpler "face-normal magnitude" form gives
a √2 speed bias at 45° and is unstable when the normal velocity is
near zero.
"""
function calvmip_exp2!(cr_x::AbstractMatrix, cr_y::AbstractMatrix,
    u_bar::AbstractMatrix, v_bar::AbstractMatrix,
    H_ice, f_ice, lsf, time::Real;
    xc::AbstractVector{<:Real},
    yc::AbstractVector{<:Real})
    wv = -300.0 * sinpi(2.0 * Float64(time) / 1000.0)

    Nxu, Nyu = size(u_bar, 1), size(u_bar, 2)   # x-faces: Nx+1, Ny
    Nxv, Nyv = size(v_bar, 1), size(v_bar, 2)   # y-faces: Nx,   Ny+1

    # x-faces: cross-stagger v from the 4 surrounding y-faces.
    @inbounds for j in 1:Nyu, i in 1:Nxu
        u = u_bar[i, j]
        i1 = max(1, i - 1);
        i2 = min(Nxv, i)
        jp1 = min(Nyv, j + 1)
        vcrs = 0.25 * (v_bar[i1, j] + v_bar[i1, jp1] +
                       v_bar[i2, j] + v_bar[i2, jp1])
        uxy = max(1e-8, sqrt(u*u + vcrs*vcrs))
        cr_x[i, j] = -u + (u / uxy) * wv
    end

    # y-faces: cross-stagger u from the 4 surrounding x-faces.
    @inbounds for j in 1:Nyv, i in 1:Nxv
        v = v_bar[i, j]
        ip1 = min(Nxu, i + 1)
        jm1 = max(1, j - 1);
        jj = min(Nyu, j)
        ucrs = 0.25 * (u_bar[i, jm1] + u_bar[ip1, jm1] +
                       u_bar[i, jj] + u_bar[ip1, jj])
        uxy = max(1e-8, sqrt(v*v + ucrs*ucrs))
        cr_y[i, j] = -v + (v / uxy) * wv
    end
    return cr_x, cr_y
end
