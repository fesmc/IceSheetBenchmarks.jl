# InitMIP

[`InitMIPBenchmark`](@ref) is a thin spec wrapper around the InitMIP
intercomparison (Goelzer et al. 2018): it reads the grid axes from a
host-provided `REGIONS.nc` NetCDF file (typically the standard
Greenland or Antarctic grids at 16 km or 32 km) and leaves all
physical fields — topography, climate, geothermal flux — to be loaded
by the host's own data I/O.

Because InitMIP uses real-data forcings rather than idealised
analytical fields, there is no spec figure here. The benchmark exists
to give InitMIP runs a uniform `AbstractBenchmark` handle alongside
the idealised benchmarks.

## Constructor

```julia
InitMIPBenchmark(regions_nc::AbstractString)
```

`regions_nc` is a path to a NetCDF file containing 1-D coordinate
variables `xc` and `yc` (in metres). The constructor reads only those
axes; nothing else.

## State

```julia
b = InitMIPBenchmark("/path/to/REGIONS.nc")
s = state(b, 0.0)   # (xc = b.xc, yc = b.yc)
```

All other fields are the host's responsibility.

## Helpers

```@docs
InitMIPBenchmark
```

## References

- Goelzer, H. et al. (2018). *Design and results of the ice sheet model
  initialisation experiments initMIP-Greenland: an ISMIP6
  intercomparison.* The Cryosphere, 12, 1433–1460.
