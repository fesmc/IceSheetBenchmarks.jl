# Regenerate the pre-rendered PNGs that the documentation pages embed.
#
# Run manually whenever a benchmark spec changes:
#
#   julia --project=docs docs/make_figs.jl
#
# Figures land in docs/src/assets/figs/ and are checked into the repo.

using IceSheetBenchmarks
using CairoMakie
using Printf

const FIG_DIR = joinpath(@__DIR__, "src", "assets", "figs")
mkpath(FIG_DIR)

CairoMakie.activate!(type = "png")

savefig(name, fig) = save(joinpath(FIG_DIR, name), fig; px_per_unit = 2)

# ----------------------------------------------------------------------
# Shared helpers
# ----------------------------------------------------------------------

"Map a 2D field over a Cartesian (km) grid as a heatmap with a colorbar."
function heatmap_panel!(ax, xc_km, yc_km, F; colormap = :viridis,
                        colorrange = nothing, hide_y = false)
    hm = if colorrange === nothing
        heatmap!(ax, xc_km, yc_km, F; colormap = colormap)
    else
        heatmap!(ax, xc_km, yc_km, F; colormap = colormap, colorrange = colorrange)
    end
    ax.xlabel = "x (km)"
    ax.ylabel = hide_y ? "" : "y (km)"
    ax.aspect = DataAspect()
    hide_y && hideydecorations!(ax; grid = false)
    return hm
end

# ----------------------------------------------------------------------
# Bueler — Halfar dome
# ----------------------------------------------------------------------

function fig_bueler_B()
    b = HalfarDomeBenchmark(:B; dx_km = 25.0)
    s0   = state(b, 0.0)
    sLat = state(b, 5_000.0)

    fig = Figure(size = (1100, 460))
    ax1 = Axis(fig[1, 1]; title = "BUELER-B  H_ice  (t = 0 yr)")
    ax2 = Axis(fig[1, 2]; title = "BUELER-B  H_ice  (t = 5 000 yr)")

    crange = (0.0, b.H0)
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3
    heatmap_panel!(ax1, xc_km, yc_km, s0.H_ice;
                   colormap = :ice, colorrange = crange)
    hm = heatmap_panel!(ax2, xc_km, yc_km, sLat.H_ice;
                        colormap = :ice, colorrange = crange, hide_y = true)
    Colorbar(fig[1, 3], hm; label = "H_ice (m)")
    savefig("bueler_B.png", fig)
end

function fig_bueler_C()
    b = HalfarDomeBenchmark(:C; dx_km = 25.0, lambda = 5.0)
    t = 1_000.0
    s = state(b, t)

    fig = Figure(size = (1100, 460))
    ax1 = Axis(fig[1, 1]; title = @sprintf("BUELER-C  H_ice  (t = %.0f yr)", t))
    ax2 = Axis(fig[1, 2]; title = @sprintf("BUELER-C  smb_ref  (t = %.0f yr)", t))
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3

    hm1 = heatmap_panel!(ax1, xc_km, yc_km, s.H_ice; colormap = :ice)
    Colorbar(fig[1, 1, Right()], hm1; label = "H_ice (m)")
    hm2 = heatmap_panel!(ax2, xc_km, yc_km, s.smb_ref;
                        colormap = :balance, hide_y = true)
    Colorbar(fig[1, 2, Right()], hm2; label = "smb_ref (m/yr)")
    savefig("bueler_C.png", fig)
end

# ----------------------------------------------------------------------
# ISMIP-HOM C
# ----------------------------------------------------------------------

function fig_hom_c()
    b = HOMCBenchmark(:C; L_km = 80.0, dx_km = 1.0)
    s = state(b, 0.0)

    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3
    beta  = [IceSheetBenchmarks._hom_c_beta(b, b.xc[i], b.yc[j])
             for i in eachindex(b.xc), j in eachindex(b.yc)]

    fig = Figure(size = (1100, 460))
    ax1 = Axis(fig[1, 1]; title = "ISMIP-HOM C  z_bed  (sloping)")
    ax2 = Axis(fig[1, 2]; title = "ISMIP-HOM C  β (basal friction)")
    hm1 = heatmap_panel!(ax1, xc_km, yc_km, s.z_bed; colormap = :terrain)
    Colorbar(fig[1, 1, Right()], hm1; label = "z_bed (m)")
    hm2 = heatmap_panel!(ax2, xc_km, yc_km, beta;
                        colormap = :viridis, hide_y = true)
    Colorbar(fig[1, 2, Right()], hm2; label = "β (Pa yr m⁻¹)")
    savefig("hom_c.png", fig)
end

# ----------------------------------------------------------------------
# Trough F17
# ----------------------------------------------------------------------

function fig_trough()
    b = TroughBenchmark(:F17; dx_km = 4.0)
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3
    z_bed = [IceSheetBenchmarks._trough_f17_zbed(xc_km[i], yc_km[j],
                                                 b.fc_km, b.dc_m, b.wc_km)
             for i in eachindex(xc_km), j in eachindex(yc_km)]

    fig = Figure(size = (1000, 380))
    ax = Axis(fig[1, 1]; title = "Trough F17  z_bed",
              xlabel = "x (km)", ylabel = "y (km)")
    ax.aspect = DataAspect()
    hm = heatmap!(ax, xc_km, yc_km, z_bed; colormap = :balance)
    Colorbar(fig[1, 2], hm; label = "z_bed (m)")
    savefig("trough.png", fig)
end

# ----------------------------------------------------------------------
# EISMINT-1 moving margin
# ----------------------------------------------------------------------

function fig_eismint_moving()
    b = EISMINT1MovingBenchmark()
    s = state(b, 0.0)
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3

    fig = Figure(size = (640, 540))
    ax = Axis(fig[1, 1]; title = "EISMINT-1 moving margin  smb_ref")
    ax.aspect = DataAspect()
    hm = heatmap!(ax, xc_km, yc_km, s.smb_ref; colormap = :balance,
                  colorrange = (-b.smb_max, b.smb_max))
    ax.xlabel = "x (km)"; ax.ylabel = "y (km)"
    Colorbar(fig[1, 2], hm; label = "smb_ref (m/yr)")
    savefig("eismint_moving.png", fig)
end

# ----------------------------------------------------------------------
# MISMIP3D
# ----------------------------------------------------------------------

function fig_mismip3d()
    b = MISMIP3DBenchmark(:Stnd; dx_km = 16.0)
    s = state(b, 0.0)
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3

    fig = Figure(size = (1100, 360))
    ax1 = Axis(fig[1, 1]; title = "MISMIP3D Stnd  z_bed  (reverse slope)")
    ax2 = Axis(fig[1, 2]; title = "MISMIP3D Stnd  H_ice  (initial slab)")
    ax1.xlabel = "x (km)"; ax1.ylabel = "y (km)"; ax1.aspect = DataAspect()
    ax2.xlabel = "x (km)"; ax2.ylabel = "y (km)"; ax2.aspect = DataAspect()

    hm1 = heatmap!(ax1, xc_km, yc_km, s.z_bed; colormap = :terrain)
    Colorbar(fig[1, 1, Right()], hm1; label = "z_bed (m)")
    hm2 = heatmap!(ax2, xc_km, yc_km, s.H_ice; colormap = :ice)
    Colorbar(fig[1, 2, Right()], hm2; label = "H_ice (m)")
    savefig("mismip3d.png", fig)
end

# ----------------------------------------------------------------------
# CalvingMIP — circular bowl (exp1/exp2) and Thule bowl (exp3/exp4)
# ----------------------------------------------------------------------

function fig_calvingmip_circular()
    b = CalvingMIPBenchmark(:exp1; dx_km = 25.0)
    s = state(b, 0.0)
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3

    fig = Figure(size = (720, 620))
    ax = Axis(fig[1, 1]; title = "CalvingMIP (exp1/exp2)  z_bed  — circular bowl")
    ax.xlabel = "x (km)"; ax.ylabel = "y (km)"; ax.aspect = DataAspect()
    hm = heatmap!(ax, xc_km, yc_km, s.z_bed; colormap = :balance)
    Colorbar(fig[1, 2], hm; label = "z_bed (m)")
    savefig("calvingmip_circular.png", fig)
end

function fig_calvingmip_thule()
    b = CalvingMIPBenchmark(:exp3; dx_km = 25.0)
    s = state(b, 0.0)
    xc_km = b.xc ./ 1e3
    yc_km = b.yc ./ 1e3

    fig = Figure(size = (720, 620))
    ax = Axis(fig[1, 1]; title = "CalvingMIP (exp3/exp4)  z_bed  — Thule bowl")
    ax.xlabel = "x (km)"; ax.ylabel = "y (km)"; ax.aspect = DataAspect()
    hm = heatmap!(ax, xc_km, yc_km, s.z_bed; colormap = :balance)
    Colorbar(fig[1, 2], hm; label = "z_bed (m)")
    savefig("calvingmip_thule.png", fig)
end

# ----------------------------------------------------------------------
# Drive
# ----------------------------------------------------------------------

function main()
    @info "Rendering figures into $FIG_DIR"
    fig_bueler_B()
    fig_bueler_C()
    fig_hom_c()
    fig_trough()
    fig_eismint_moving()
    fig_mismip3d()
    fig_calvingmip_circular()
    fig_calvingmip_thule()
    @info "Done."
end

main()
