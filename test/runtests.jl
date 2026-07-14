using IceSheetBenchmarks
using NCDatasets
using Test

const ISB = IceSheetBenchmarks

# Cell area (m²) for a uniform axis pair, used in volume checks.
cell_area(b) = (b.xc[2] - b.xc[1]) * (b.yc[2] - b.yc[1])

@testset "IceSheetBenchmarks" begin

    # ------------------------------------------------------------------
    @testset "Bueler / Halfar" begin
        @testset "constructor validation" begin
            @test_throws ErrorException HalfarDomeBenchmark(:X; dx_km = 25.0)
            @test_throws ErrorException HalfarDomeBenchmark(:C; dx_km = 25.0)            # needs lambda
            @test_throws ErrorException HalfarDomeBenchmark(:B; dx_km = 25.0, lambda = 5.0)  # B rejects lambda
            @test_throws ErrorException HalfarDomeBenchmark(:B; dx_km = 7.0)             # non-integer grid
        end

        b = HalfarDomeBenchmark(:B; dx_km = 25.0)
        Nx, Ny = length(b.xc), length(b.yc)
        @test Nx == Ny == Int(2 * b.R0_km / 25.0) + 1

        s0 = state(b, 0.0)
        @test size(s0.H_ice) == (Nx, Ny)
        @test all(s0.H_ice .>= 0.0)
        @test all(s0.z_bed .== 0.0)
        @test all(s0.smb_ref .== 0.0)                 # variant :B has no smb
        @test s0.H_ice[(Nx+1)÷2, (Ny+1)÷2] ≈ maximum(s0.H_ice)  # summit at centre

        @testset "volume conservation (Halfar B is source-free)" begin
            V0 = sum(state(b, 0.0).H_ice) * cell_area(b)
            V1 = sum(state(b, 2_000.0).H_ice) * cell_area(b)
            @test isapprox(V0, V1; rtol = 0.02)
        end

        @testset "analytical velocity" begin
            ux, uy = analytical_velocity(b, 100.0)
            @test size(ux) == (Nx + 1, Ny)
            @test size(uy) == (Nx, Ny + 1)
            @test all(isfinite, ux) && all(isfinite, uy)
            # Radial symmetry: outward flow flips sign across the summit.
            ic = (Nx + 1) ÷ 2
            @test ux[ic + 2, ic] ≈ -ux[ic - 1, ic] atol = 1e-6
        end

        @testset "closed-form dH/dr vs finite differences" begin
            t = 500.0
            h = 1.0e2
            for r in (1.0e5, 3.0e5, 5.0e5)
                Hp = ISB._halfar_HR_dHdr(b, r + h, t)[1]
                Hm = ISB._halfar_HR_dHdr(b, r - h, t)[1]
                fd = (Hp - Hm) / (2h)
                closed = ISB._halfar_dHdr_closed(b, r, t)
                @test isapprox(fd, closed; rtol = 1e-3, atol = 1e-6)
            end
        end

        @testset "variant :C has analytical smb" begin
            bc = HalfarDomeBenchmark(:C; dx_km = 25.0, lambda = 5.0)
            sc = state(bc, 1_000.0)
            @test any(sc.smb_ref .> 0.0)
            @test all(isfinite, sc.smb_ref)
        end
    end

    # ------------------------------------------------------------------
    @testset "ISMIP-HOM C" begin
        @test_throws ErrorException HOMCBenchmark(:A)
        @test_throws ErrorException HOMCBenchmark(:C; L_km = 80.0, dx_km = 7.0)  # non-integer

        b = HOMCBenchmark(:C; L_km = 80.0, dx_km = 2.0)
        Nx, Ny = length(b.xc), length(b.yc)
        @test Nx == Ny == 40
        s = state(b, 0.0)
        @test all(s.H_ice .== b.H)                    # uniform slab
        @test all(s.z_sl .< 0.0)                      # forced grounded
        # bed decreases linearly in x.
        @test s.z_bed[1, 1] > s.z_bed[end, 1]
        # β oscillates around β0 with amplitude β_amp.
        βmid = ISB._hom_c_beta(b, b.L_km * 1e3 / 4, b.L_km * 1e3 / 4)
        @test βmid > b.beta0
    end

    # ------------------------------------------------------------------
    @testset "Trough F17" begin
        @test_throws ErrorException TroughBenchmark(:X)
        b = TroughBenchmark(:F17; dx_km = 4.0)
        s = state(b, 0.0)
        @test size(s.z_bed) == (length(b.xc), length(b.yc))
        @test all(s.H_ice .== 0.0)
        @test all(s.z_bed .>= -720.0)                 # clamped to deep floor
        @test all(s.T_srf .≈ b.Tsrf_const + 273.15)   # °C → K
        @test_throws ErrorException state(b, 100.0)    # only t = 0
    end

    # ------------------------------------------------------------------
    @testset "EISMINT-1 moving margin" begin
        b = EISMINT1MovingBenchmark()
        s = state(b, 0.0)
        @test all(s.H_ice .== 0.0)
        @test all(s.z_bed .== 0.0)
        @test maximum(s.smb_ref) ≈ b.smb_max          # cap saturates at summit
        @test any(s.smb_ref .< 0.0)                    # negative beyond R_el
        @test_throws ErrorException state(b, 10.0)
        @test eismint_moving_smb(b, b.x_summit, b.y_summit) ≈ b.smb_max
    end

    # ------------------------------------------------------------------
    @testset "MISMIP3D" begin
        @test_throws ErrorException MISMIP3DBenchmark(:Foo)
        b = MISMIP3DBenchmark(:Stnd; dx_km = 16.0)
        Nx, Ny = length(b.xc), length(b.yc)
        @test (Nx, Ny) == (51, 7)                      # odd Ny → centreline row
        s = state(b, 0.0)
        @test s.z_bed[1, 1] > s.z_bed[end, 1]          # reverse slope
        @test all(s.mask_ice[Nx, :] .== 0.0)           # eastern calving column
        @test all(x -> x == 0.0 || x == b.H0, s.H_ice) # slab or empty
        @test_throws ErrorException state(b, 1.0)
    end

    # ------------------------------------------------------------------
    @testset "CalvingMIP" begin
        @test_throws ErrorException CalvingMIPBenchmark(:exp9)

        @test calvmip_bed_circular(0.0, 0.0) ≈ 900.0   # bowl centre = Bc
        @test calvmip_bed_circular(800e3, 0.0) ≈ -2000.0 atol = 1e-6  # rim = Bl

        b = CalvingMIPBenchmark(:exp1; dx_km = 25.0)
        @test length(b.xc) == length(b.yc) == 64
        @test b.domain == :circular
        @test CalvingMIPBenchmark(:exp3).domain == :thule
        s = state(b, 0.0)
        @test all(s.H_ice .== 0.0)
        @test all(s.lsf .== 1.0)
        @test_throws ErrorException state(b, 5.0)

        @testset "exp1 pins the front at r_lim" begin
            Nx = length(b.xc); Ny = length(b.yc)
            cr_x = zeros(Nx + 1, Ny); cr_y = zeros(Nx, Ny + 1)
            u = fill(100.0, Nx + 1, Ny); v = fill(100.0, Nx, Ny + 1)
            calvmip_exp1!(cr_x, cr_y, u, v, nothing, nothing, nothing, 0.0;
                          xc = b.xc, yc = b.yc, r_lim = 750e3)
            # Deep-interior x-faces (both neighbours well inside) are zeroed.
            ic = Nx ÷ 2
            @test cr_x[ic, ic] == 0.0
        end

        @testset "exp2 oscillation is bounded" begin
            Nx = length(b.xc); Ny = length(b.yc)
            cr_x = zeros(Nx + 1, Ny); cr_y = zeros(Nx, Ny + 1)
            u = fill(50.0, Nx + 1, Ny); v = fill(50.0, Nx, Ny + 1)
            calvmip_exp2!(cr_x, cr_y, u, v, nothing, nothing, nothing, 250.0;
                          xc = b.xc, yc = b.yc)
            @test all(isfinite, cr_x) && all(isfinite, cr_y)
        end
    end

    # ------------------------------------------------------------------
    @testset "InitMIP" begin
        mktempdir() do dir
            regions = joinpath(dir, "REGIONS.nc")
            NCDataset(regions, "c") do ds
                defDim(ds, "xc", 4); defDim(ds, "yc", 3)
                xv = defVar(ds, "xc", Float64, ("xc",)); xv[:] = [0.0, 16.0, 32.0, 48.0]
                xv.attrib["units"] = "km"
                yv = defVar(ds, "yc", Float64, ("yc",)); yv[:] = [0.0, 16.0, 32.0]
                yv.attrib["units"] = "km"
            end
            b = InitMIPBenchmark(regions)
            @test b.xc[end] ≈ 48_000.0                 # km → m conversion
            @test state(b, 0.0) == (xc = b.xc, yc = b.yc)
            @test_throws ErrorException write_fixture!(b, joinpath(dir, "x.nc"))
        end
    end

    # ------------------------------------------------------------------
    @testset "analytical_velocity fallback" begin
        @test_throws ErrorException analytical_velocity(EISMINT1MovingBenchmark(), 0.0)
    end

    # ------------------------------------------------------------------
    @testset "parametric element type {T}" begin
        # Default construction is Float64.
        @test HalfarDomeBenchmark(:B; dx_km = 50.0) isa HalfarDomeBenchmark{Float64}
        @test EISMINT1MovingBenchmark() isa EISMINT1MovingBenchmark{Float64}

        # All-Float32 inputs give a Float32 benchmark with Float32 axes.
        b32 = HalfarDomeBenchmark(:B; dx_km = 50.0f0, R0_km = 750.0f0, H0 = 3600.0f0,
                              n = 3.0f0, A = 1.0f-16, rho_ice = 910.0f0, g = 9.81f0)
        @test b32 isa HalfarDomeBenchmark{Float32}
        @test eltype(b32.xc) == Float32

        e32 = EISMINT1MovingBenchmark(; dx_km = 50.0f0, L_km = 1500.0f0,
                                      R_el_km = 450.0f0, smb_max = 0.5f0,
                                      smb_grad = 0.01f0, T_srf_const = 270.0f0,
                                      Q_geo_const = 42.0f0, n_glen = 3.0f0,
                                      A_glen = 1.0f-16)
        @test e32 isa EISMINT1MovingBenchmark{Float32}
        @test eltype(e32.xc) == Float32 && typeof(e32.x_summit) == Float32

        # Mixed / integer inputs promote to a common float type (no error).
        @test MISMIP3DBenchmark(:Stnd; dx_km = 16) isa MISMIP3DBenchmark{Float64}
        @test CalvingMIPBenchmark(:exp1; dx_km = 25) isa CalvingMIPBenchmark{Float64}
    end

    # ------------------------------------------------------------------
    @testset "fixture round-trip" begin
        cases = (
            HalfarDomeBenchmark(:B; dx_km = 50.0),
            HOMCBenchmark(:C; L_km = 80.0, dx_km = 4.0),
            TroughBenchmark(:F17; dx_km = 8.0),
            EISMINT1MovingBenchmark(),
            MISMIP3DBenchmark(:Stnd; dx_km = 16.0),
            CalvingMIPBenchmark(:exp1; dx_km = 50.0),
        )
        mktempdir() do dir
            for (i, b) in enumerate(cases)
                path = joinpath(dir, "fixture_$i.nc")
                out = write_fixture!(b, path; times = [0.0])
                @test out == [path]
                @test isfile(path)
                s = state(b, 0.0)
                NCDataset(path) do ds
                    @test ds.dim["xc"] == length(b.xc)
                    @test ds.dim["yc"] == length(b.yc)
                    @test ds["xc"][:] ≈ b.xc ./ 1e3     # stored in km
                    @test ds["H_ice"][:, :] ≈ s.H_ice
                    @test haskey(ds.attrib, "benchmark")
                end
            end
        end
    end

end
