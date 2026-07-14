# ----------------------------------------------------------------------
# Bueler analytical ice-flow solutions, ported from Yelmo Fortran
# (`yelmo/tests/ice_benchmarks.f90`), and the corresponding
# `HalfarDomeBenchmark <: AbstractBenchmark` implementation of the
# AbstractBenchmark interface.
#
# Math layer:
#   - `bueler_gamma`    : SIA prefactor `γ = 2 A (ρ_i g)^n / (n + 2)`.
#   - `bueler_test_BC!` : Halfar (1981) similarity solution for the
#       time-dependent isothermal-SIA dome (Bueler et al. 2005,
#       Eqs. 10–11). Sets `H_ice` and `mbal` in-place at every
#       (xc, yc, time) point. `lambda = 0` reproduces pure decay
#       (BUELER-B); `lambda > 0` gives the BUELER-C variant with an
#       analytical mass balance.
#
# Benchmark layer:
#   - `HalfarDomeBenchmark`  : carries grid axes + Halfar parameters,
#       implements `state` / `write_fixture!` / `analytical_velocity`.
#
# Units: `xc`/`yc` in metres, `H` in metres, `mbal` in m/yr,
# `time` in years, `R0` in km, `H0` in metres, `A` in Pa^-3 yr^-1
# (the Yelmo Fortran convention; the Halfar formula's `time` in the
# similarity exponents is in the same units as `A`'s time, hence yr).
# ----------------------------------------------------------------------

# ----------------------------------------------------------------------
# Math layer (ported verbatim from the previous bueler.jl).
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

SIA flow prefactor `γ = 2 A (ρ_i g)^n / (n + 2)`. With Yelmo defaults
(`A = 1e-16` Pa⁻³ yr⁻¹, `n = 3`, `ρ_i = 910`, `g = 9.81`) this gives
`γ ≈ 2.84 × 10⁻²` m⁻³ yr⁻¹ (Bueler et al. 2005 default).
"""
@inline bueler_gamma(A, n, rho_ice, g) = 2.0 * A * (rho_ice * g)^n / (n + 2.0)

"""
$(TYPEDSIGNATURES)

Halfar similarity solution for the isothermal-SIA dome. Mutates
`H_ice` (m) and `mbal` (m/yr) in-place at every (xc[i], yc[j]) point.

`H_ice` and `mbal` must be 2D arrays of size `(length(xc), length(yc))`.
`xc` and `yc` are cell-centre coordinates in **metres**.

`time` is the elapsed time **after** the analytical reference
(`t = 0` corresponds to the initial dome at the start of the
simulation, **not** an absolute Halfar `t₀`). The internal `t₀`
shift is added inside the formula via the (R0, H0) dome scale.

`lambda = 0` reproduces BUELER-B (pure Halfar decay, no mass balance).
`lambda > 0` reproduces BUELER-C (with analytical `mbal = λ H / t`).

Port of `yelmo/tests/ice_benchmarks.f90:167 bueler_test_BC`.
"""
function bueler_test_BC!(H_ice::AbstractMatrix{<:Real},
    mbal::AbstractMatrix{<:Real},
    xc::AbstractVector{<:Real},
    yc::AbstractVector{<:Real},
    time::Real;
    R0::Real=750.0, H0::Real=3600.0,
    lambda::Real=0.0, n::Real=3.0,
    A::Real=1e-16,
    rho_ice::Real=910.0, g::Real=9.81)
    Nx, Ny = length(xc), length(yc)
    size(H_ice) == (Nx, Ny) ||
        error("bueler_test_BC!: H_ice has shape $(size(H_ice)); expected ($Nx, $Ny).")
    size(mbal) == (Nx, Ny) ||
        error("bueler_test_BC!: mbal has shape $(size(mbal)); expected ($Nx, $Ny).")

    # Convert R0 from km → m to match the Cartesian xc/yc grid.
    R0_m = Float64(R0) * 1e3

    # Halfar similarity exponents.
    α = (2.0 - (n + 1.0) * lambda) / (5.0 * n + 3.0)
    β = (1.0 + (2.0 * n + 1.0) * lambda) / (5.0 * n + 3.0)
    γ = bueler_gamma(A, n, rho_ice, g)
    t0 = (β / γ) * ((2.0 * n + 1.0) / (n + 1.0))^n *
         (R0_m^(n + 1) / H0^(2.0 * n + 1.0))
    t1 = Float64(time) + t0

    @inbounds for j in 1:Ny, i in 1:Nx
        r = sqrt(Float64(xc[i])^2 + Float64(yc[j])^2)
        # Halfar profile, Bueler 2005 Eqs. 10–11. The `max(0, ⋯)`
        # clamps the flank to zero outside the moving margin.
        fac = max(0.0, 1.0 - ((t1 / t0)^(-β) * r / R0_m)^((n + 1.0) / n))
        H_ice[i, j] = H0 * (t1 / t0)^(-α) * fac^(n / (2.0 * n + 1.0))
        mbal[i, j] = (lambda / t1) * H_ice[i, j]
    end

    return H_ice, mbal
end

# ----------------------------------------------------------------------
# HalfarDomeBenchmark — concrete AbstractBenchmark implementation.
# ----------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Halfar dome benchmark (Bueler et al. 2005).

# Fields:
 - `variant` : Benchmark variant.
    - :B = pure Halfar decay (no mass balance, `λ = 0`). The default. or :C)
    - `:C` — Halfar + analytical mass balance `mbal = λ H / t`. Requires `lambda` to be passed explicitly.
 - `xc`, `yc` : cell-centre x/y coordinates (m)
 - `R0_km` : Halfar dome radius (km). Default 750 km.
 - `H0` : Halfar dome height (m). Default 3600 m.
 - `lambda` : mass-balance scale (dimensionless). Default 0.0 for variant :B; must be > 0 for variant :C.
 - `n` : Glen exponent. Default 3.
 - `A` : Glen flow factor (Pa⁻³ yr⁻¹). Default 1e-16 Pa⁻³ yr⁻¹.
 - `rho_ice` : ice density (kg/m³). Default 910 kg/m³.
 - `g` : gravitational acceleration (m/s²). Default 9.81 m/s².
"""
struct HalfarDomeBenchmark{T<:AbstractFloat} <: AbstractBenchmark
    variant::Symbol
    xc::Vector{T}
    yc::Vector{T}
    R0_km::T
    H0::T
    lambda::T
    n::T
    A::T
    rho_ice::T
    g::T
end

function HalfarDomeBenchmark(variant::Symbol;
    dx_km::Real,
    R0_km::Real=750.0,
    H0::Real=3600.0,
    lambda=nothing,
    n::Real=3.0,
    A::Real=1e-16,
    rho_ice::Real=910.0,
    g::Real=9.81)
    variant in (:B, :C) || error(
        "HalfarDomeBenchmark: variant must be :B or :C, got $(variant).")

    if variant === :B
        lambda === nothing || lambda == 0.0 ||
            error("HalfarDomeBenchmark(:B): variant B has lambda=0 by definition; do not pass lambda kwarg.")
        lam = 0.0
    else  # :C
        lambda === nothing &&
            error("HalfarDomeBenchmark(:C): variant C requires `lambda` keyword (mass-balance scale).")
        lambda > 0.0 ||
            error("HalfarDomeBenchmark(:C): lambda must be > 0; got $(lambda).")
        lam = Float64(lambda)
    end

    ratio = 2 * R0_km / dx_km
    isinteger(ratio) ||
        error("HalfarDomeBenchmark: 2*R0_km/dx_km must be integer (got $ratio).")
    Nx = Int(ratio) + 1
    half_m = R0_km * 1e3
    xc = collect(range(-half_m, half_m; length=Nx))
    yc = copy(xc)

    return HalfarDomeBenchmark(variant, xc, yc, Float64(R0_km), Float64(H0),
        lam, Float64(n), Float64(A),
        Float64(rho_ice), Float64(g))
end

"""
$(TYPEDSIGNATURES)

Analytical Halfar state at time `t`. Returns a NamedTuple with:

  - `xc`, `yc`     : grid axes in metres (echoes `b.xc` / `b.yc`).
  - `H_ice`        : 2D analytical ice thickness from `bueler_test_BC!`.
  - `z_bed`        : 2D bedrock elevation (flat at zero).
  - `smb_ref`      : 2D analytical surface mass balance (zero for
                     `variant = :B`, `λ H / t` for `:C`).

The keys other than `xc` / `yc` map onto Yelmo schema variables and
are routed into the appropriate component group by the generic
`YelmoModel(::AbstractBenchmark, t)` constructor provided by the
`YelmoBenchmarks` package extension.
"""
function state(b::HalfarDomeBenchmark, t::Real)
    Nx = length(b.xc)
    Ny = length(b.yc)
    H = zeros(Nx, Ny)
    smb = zeros(Nx, Ny)
    bueler_test_BC!(H, smb, b.xc, b.yc, Float64(t);
        R0=b.R0_km,
        H0=b.H0,
        lambda=b.lambda,
        n=b.n,
        A=b.A,
        rho_ice=b.rho_ice,
        g=b.g)
    z_bed = zeros(Nx, Ny)
    return (xc=b.xc, yc=b.yc,
        H_ice=H, z_bed=z_bed, smb_ref=smb)
end

"""
$(TYPEDSIGNATURES)

Serialize the analytical Halfar state at each `t` in `times` to a
NetCDF restart at `path`. Single-time only for now — multi-time
fixtures (a `time` dimension with multiple snapshots) are deferred to
a future milestone.

Returns a 1-element `Vector{String}` containing `path`. The restart
schema matches the previously-committed `bueler_b_smoke__t1000.nc`
(uniform 11-point ice / 5-point rock sigma axes).
"""
function write_fixture!(b::HalfarDomeBenchmark, path::AbstractString;
    times::AbstractVector{<:Real}=[0.0])
    length(times) == 1 ||
        error("write_fixture!(HalfarDomeBenchmark, …): multi-time fixtures " *
              "deferred to a future milestone (got $(length(times)) times).")
    t = Float64(first(times))
    s = state(b, t)

    vars = (
        ("H_ice", s.H_ice, "m", "Ice thickness (analytical Halfar)"),
        ("smb_ref", s.smb_ref, "m/yr", "Surface mass balance (analytical)"),
        ("z_bed", s.z_bed, "m", "Bedrock elevation (flat)"),
    )
    attrs = (
        "benchmark" => "BUELER-$(string(b.variant))",
        "solution_type" => "analytical-halfar",
        "time_yr" => t,
        "R0_km" => b.R0_km,
        "H0_m" => b.H0,
        "lambda" => b.lambda,
        "n_glen" => b.n,
        "A_Pa-3yr-1" => b.A,
    )
    return _write_restart!(path, b.xc, b.yc,
        _DEFAULT_ZETA_AC, _DEFAULT_ZETA_ROCK_AC,
        vars, attrs)
end

# ----------------------------------------------------------------------
# Private helpers for the closed-form Halfar depth-averaged velocity.
#
# Derivation: for the BUELER-B (lambda = 0) Halfar dome (Bueler et al.
# 2005, J. Glaciol. 51 (173), Eqs. 9-11; original similarity solution
# Halfar 1981, JGR 86 (C11), pp. 11065-11072):
#
#   alpha = (2 - (n+1) lambda) / (5 n + 3),
#   beta  = (1 + (2 n + 1) lambda) / (5 n + 3),
#   gamma = 2 A (rho_i g)^n / (n + 2),
#   t0    = (beta / gamma) * ((2 n + 1)/(n + 1))^n * R0^(n+1) / H0^(2 n + 1).
#
# Following the Yelmo Fortran convention in `bueler_test_BC` (and
# matching Bueler 2005 Eq. 11 with the time variable shifted to
# `t -> t + t0`), the Halfar profile is parameterised by
#
#   u(r, t) = (t1/t0)^(-beta) * r / R0,    with t1 = t + t0,
#
# and for u < 1
#
#   H(r, t) = H0 * (t1/t0)^(-alpha) * [1 - u^((n+1)/n)]^(n/(2 n + 1)),
#
# giving (by direct r-differentiation, q = n/(2n+1), p = (n+1)/n,
# u_r = (t1/t0)^(-beta)/R0)
#
#   dH/dr(r, t) = -H0 * (t1/t0)^(-alpha) * (t1/t0)^(-beta) / R0 *
#                  (n+1)/(2 n + 1) * u^(1/n) *
#                  [1 - u^((n+1)/n)]^(n/(2 n + 1) - 1).
#
# At u >= 1 the dome margin has been crossed: H = dH/dr = 0.
#
# The depth-averaged SIA velocity for a flat-bed Halfar dome is
#
#   bar(u_r) = -gamma * H^(n+1) * |dH/dr|^(n-1) * dH/dr   (radial),
#
# converted to Cartesian as bar(u_x) = (x/r) * bar(u_r), bar(u_y) =
# (y/r) * bar(u_r). At r = 0 the symmetry forces bar(u) = 0.
# ----------------------------------------------------------------------

@inline function _halfar_exponents(n::Real, lambda::Real)
    α = (2.0 - (Float64(n) + 1.0) * Float64(lambda)) / (5.0 * Float64(n) + 3.0)
    β = (1.0 + (2.0 * Float64(n) + 1.0) * Float64(lambda)) / (5.0 * Float64(n) + 3.0)
    return α, β
end

@inline function _halfar_t0(b::HalfarDomeBenchmark)
    R0_m = b.R0_km * 1e3
    α, β = _halfar_exponents(b.n, b.lambda)
    γ = bueler_gamma(b.A, b.n, b.rho_ice, b.g)
    t0 = (β / γ) * ((2.0 * b.n + 1.0) / (b.n + 1.0))^b.n *
         (R0_m^(b.n + 1) / b.H0^(2.0 * b.n + 1.0))
    return α, β, γ, t0, R0_m
end

# Closed-form (H, dH/dr) at radius `r` (m) and time `t` (yr) for the
# Halfar dome described by `b`. Returns `(0.0, 0.0)` outside the moving
# margin (u >= 1).
function _halfar_HR_dHdr(b::HalfarDomeBenchmark, r::Real, t::Real)
    α, β, γ, t0, R0_m = _halfar_t0(b)
    t1 = Float64(t) + t0
    # Yelmo Fortran convention (matches `bueler_test_BC!` above):
    # `u = (t1/t0)^(-β) * r / R0`. As t increases, (t1/t0)^(-β)
    # decreases (β > 0), so the margin radius r_margin = R0 (t1/t0)^β
    # increases — the Halfar dome spreads with time.
    u = (t1 / t0)^(-β) * Float64(r) / R0_m
    if u >= 1.0
        return (0.0, 0.0)
    end
    H = b.H0 * (t1 / t0)^(-α) * (1.0 - u^((b.n + 1.0) / b.n))^(b.n / (2.0 * b.n + 1.0))
    prefac = -b.H0 * (t1 / t0)^(-α) * (t1 / t0)^(-β) / R0_m *
             (b.n + 1.0) / (2.0 * b.n + 1.0)
    dHdr = prefac * u^(1.0 / b.n) *
           (1.0 - u^((b.n + 1.0) / b.n))^(b.n / (2.0 * b.n + 1.0) - 1.0)
    return (H, dHdr)
end

# Convenience: closed-form dH/dr only (used by the unit test that
# cross-checks against numerical centred-differences of H).
_halfar_dHdr_closed(b::HalfarDomeBenchmark, r::Real, t::Real) =
    _halfar_HR_dHdr(b, r, t)[2]

"""
$(TYPEDSIGNATURES)

Closed-form depth-averaged Halfar velocity at time `t` (years) for the
Halfar dome described by `b`. Returns face-staggered 2D arrays:

  - `ux_bar` of shape `(Nx + 1, Ny)` — values at x-faces, layout
    matching `interior(y.dyn.ux_bar)[:, :, 1]` for an `XFaceField` on a
    `(Nx, Ny)`-Center grid.
  - `uy_bar` of shape `(Nx, Ny + 1)` — values at y-faces, matching
    `interior(y.dyn.uy_bar)[:, :, 1]` for a `YFaceField`.

The first/last x-face row (`ux_bar[1, :]`, `ux_bar[Nx+1, :]`) and the
first/last y-face row (`uy_bar[:, 1]`, `uy_bar[:, Ny+1]`) are left at
zero — these face cells lie outside the dome margin at all benchmark
resolutions and are not part of any margin-masked error metric.

Derivation: see the comment block above this method. References:

  - Halfar, P. 1981. On the dynamics of the ice sheets. JGR 86 (C11),
    pp. 11065-11072.
  - Bueler, E. et al. 2005. Exact solutions and verification of
    numerical models for isothermal ice sheets. J. Glaciol. 51 (173),
    pp. 291-306. Eqs. 9-11 and Section 4.
"""
function analytical_velocity(b::HalfarDomeBenchmark, t::Real)
    Nx, Ny = length(b.xc), length(b.yc)
    ux_bar = zeros(Nx + 1, Ny)
    uy_bar = zeros(Nx, Ny + 1)

    α, β, γ, t0, R0_m = _halfar_t0(b)

    # ux_bar at face position (xc[i] + xc[i+1])/2, yc[j].
    for j in 1:Ny, i in 1:(Nx-1)
        x_f = 0.5 * (b.xc[i] + b.xc[i+1])
        y_f = b.yc[j]
        r = sqrt(x_f^2 + y_f^2)
        H, dHdr = _halfar_HR_dHdr(b, r, t)
        if r < 1e-9 || H == 0.0
            ux_bar[i+1, j] = 0.0
        else
            grad_mag = abs(dHdr)
            ux_bar[i+1, j] = -γ * H^(b.n + 1.0) * grad_mag^(b.n - 1.0) *
                             (x_f / r) * dHdr
        end
    end

    # uy_bar at face position xc[i], (yc[j] + yc[j+1])/2.
    for j in 1:(Ny-1), i in 1:Nx
        x_f = b.xc[i]
        y_f = 0.5 * (b.yc[j] + b.yc[j+1])
        r = sqrt(x_f^2 + y_f^2)
        H, dHdr = _halfar_HR_dHdr(b, r, t)
        if r < 1e-9 || H == 0.0
            uy_bar[i, j+1] = 0.0
        else
            grad_mag = abs(dHdr)
            uy_bar[i, j+1] = -γ * H^(b.n + 1.0) * grad_mag^(b.n - 1.0) *
                             (y_f / r) * dHdr
        end
    end

    return (ux_bar, uy_bar)
end
