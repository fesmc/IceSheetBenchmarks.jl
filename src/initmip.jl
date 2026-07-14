# InitMIPBenchmark — real-data ice-sheet benchmark (Greenland or Antarctica).
#
# This is a thin grid-spec struct: it reads xc/yc from a REGIONS.nc
# file and exposes them so the YelmoBenchmarks extension can build the
# Oceananigans grid. All physical fields (topography, climate, GHF) are
# loaded separately in run.jl via init_topo_load! / data_load!.

Base.@kwdef struct InitMIPBenchmark{T<:AbstractFloat} <: AbstractBenchmark
    xc::Vector{T}   # cell-centre x [m]
    yc::Vector{T}   # cell-centre y [m]
end

"""
$(TYPEDSIGNATURES)

Read `xc`/`yc` coordinates from a REGIONS NetCDF file and construct the
benchmark. Coordinates are converted to metres if the file stores them
in kilometres. Works for any real-data domain (e.g. GRL-16KM,
ANT-32KM).
"""
function InitMIPBenchmark(regions_nc::AbstractString)
    NCDataset(regions_nc) do ds
        xc = Vector{Float64}(ds["xc"][:])
        yc = Vector{Float64}(ds["yc"][:])
        xu = lowercase(strip(get(ds["xc"].attrib, "units", "")))
        yu = lowercase(strip(get(ds["yc"].attrib, "units", "")))
        (xu == "km" || xu == "kilometers") && (xc .*= 1e3)
        (yu == "km" || yu == "kilometers") && (yc .*= 1e3)
        return InitMIPBenchmark(xc, yc)
    end
end

"""
$(TYPEDSIGNATURES)

Grid-only state `(xc, yc)`. InitMIP uses real-data forcing (topography,
climate, geothermal flux) that the host loads via its own data I/O, so
no physical fields are produced here at any `t`.
"""
function state(b::InitMIPBenchmark, ::Real)
    return (xc = b.xc, yc = b.yc)
end

"""
$(TYPEDSIGNATURES)

Not supported: InitMIP is a real-data benchmark with no analytical
state to serialise. Construct it from a `REGIONS.nc` file and let the
host load the physical fields directly. Throws an informative error.
"""
function write_fixture!(::InitMIPBenchmark, path::AbstractString;
                        times::AbstractVector{<:Real} = [0.0])
    error("write_fixture!(InitMIPBenchmark, …): InitMIP has no analytical " *
          "fixture — it is a real-data benchmark. Load the physical fields " *
          "from the source dataset via the host's data I/O instead.")
end
