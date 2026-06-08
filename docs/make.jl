using Documenter
using IceSheetBenchmarks

DocMeta.setdocmeta!(IceSheetBenchmarks, :DocTestSetup,
                    :(using IceSheetBenchmarks); recursive=true)

makedocs(;
    modules  = [IceSheetBenchmarks],
    authors  = "Alexander Robinson <alexander.robinson@awi.de> and contributors",
    sitename = "IceSheetBenchmarks.jl",
    format   = Documenter.HTML(;
        canonical    = "https://fesmc.github.io/IceSheetBenchmarks.jl",
        edit_link    = "main",
        assets       = String[],
        prettyurls   = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home"        => "index.md",
        "Interface"   => "interface.md",
        "Benchmarks"  => [
            "Bueler (Halfar dome)"      => "benchmarks/bueler.md",
            "ISMIP-HOM C"               => "benchmarks/hom_c.md",
            "Trough (F17)"              => "benchmarks/trough.md",
            "EISMINT-1 moving margin"   => "benchmarks/eismint_moving.md",
            "MISMIP3D"                  => "benchmarks/mismip3d.md",
            "CalvingMIP"                => "benchmarks/calvingmip.md",
            "InitMIP"                   => "benchmarks/initmip.md",
        ],
        "API reference" => "api.md",
    ],
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(;
    repo        = "github.com/fesmc/IceSheetBenchmarks.jl",
    devbranch   = "main",
    push_preview = true,
)
