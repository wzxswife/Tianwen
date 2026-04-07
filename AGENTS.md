# Tianwen Repository Agent Guidelines

Julia codebase for space physics data analysis — MAVEN spacecraft and Tianwen-1 mission magnetic field & plasma data.

## Project Structure

```
Tianwen/
├── src/
│   ├── scripts/              # Reusable modules
│   │   ├── MAVEN_load.jl     # MAVEN data loading (CDF, STS, Fortran binary, JLD2)
│   │   ├── MAVEN_plot.jl     # MAVEN plotting (VDF, PAD, heatmap, orbit)
│   │   ├── MAVEN_SWIA.jl     # SWIA instrument data
│   │   ├── MAVEN_SWEA.jl     # SWEA instrument data
│   │   ├── MAVEN_STATIC.jl   # STATIC instrument data
│   │   ├── TW_load.jl        # Tianwen-1 data loading
│   │   ├── TW_plot.jl        # Tianwen-1 plotting (bowshock, wavelet, orbit)
│   │   ├── MVA_plot.jl       # Minimum Variance Analysis
│   │   ├── wave_caculate.jl  # Wave analysis (PSD, MVA, dominant frequency)
│   │   └── TCWavelet.jl      # Custom wavelet transform (Torrence & Compo)
│   ├── DemoMAVEN.jl          # MAVEN entry point
│   ├── DemoTW.jl             # Tianwen-1 entry point
│   └── *.jl / *.ipynb        # Standalone analysis scripts
├── Project/
│   ├── IonBeam/              # Ion beam VDF analysis
│   └── TianwenData/          # Tianwen-specific processing (bowshock, wave stats, polar)
├── Data/                     # Raw data (gitignored)
├── Results/                  # Output plots/CSVs (gitignored)
├── Doc/                      # Reference docs, PDFs, skeleton tables
└── workflows/                # (gitignored)
```

## Running Scripts

```bash
# No Project.toml exists — scripts use global Julia environment
julia src/DemoMAVEN.jl
julia src/DemoTW.jl

# Run individual script
julia src/scripts/MVA_plot.jl
```

## Critical Gotchas

### No Project.toml
This repo has **no** `Project.toml`. All packages must be installed in the global environment. Do NOT add `--project=.` or try to instantiate.

### Hardcoded Path in TW_load.jl
`TW_load.jl` line 5: `root_path = "E:/Tianwen-1/"` — hardcoded Windows path. If running on a different machine, this must be updated.

### Data Directory Case Mismatch
`.gitignore` lists `results` (lowercase) but actual directory is `Results/`. Same for `Data/`. Data is expected at:
- MAVEN: `Data/MAVEN/` (STS files)
- Tianwen: `Data/32Hz/` (DAT files)

### MAVEN_load.jl Quirk
Lines 8-18 have `using` statements **outside** the `module MAVEN_load` block. The module then re-declares its own `using` inside. This is intentional — the outer imports are for top-level test usage.

## Data Patterns

### Data Container
All data loading returns `Dict{Symbol, Any}()` with these standard keys:
```julia
data = Dict{Symbol,Any}(
    :epoch => times,           # DateTime vector
    :B => B,                   # N×3 magnetic field matrix (nT)
    :position => position,     # N×3 position matrix (km, MSO coords)
    :B_total => B_total,       # N-element total field magnitude
    :data_load_flag => true,   # false on load failure
)
```

### Position Normalization
Position is loaded in km, normalized to Mars radii **after** loading:
```julia
data[:position] ./= 3390.0   # or 3393.5 depending on context
```
Mars radius constants used across codebase: `3390.0`, `3393.5`, `3389.5` — pick the one matching the module you're in.

### Coordinate System
All positions in **MSO** (Mars Solar Orbital): X sunward, Y opposite orbital motion, Z completes right-hand system.

### Time Range Slicing
Standard pattern to extract a time window:
```julia
time_range = shock_time - Dates.Minute(1) .+ Dates.Minute(1) .* range(0, 3)
mag_data = find_avail_data(data, time_range, [:epoch, :B])
```
`find_avail_data` (in `TW_load.jl`) also filters NaN rows from position and B.

### File Formats
| Format | Used For |
|--------|----------|
| `.sts` | MAVEN mag L2 (fixed-width text, skip header to "END" line) |
| `.dat` | Tianwen MOMAG (space-delimited, skip 19 header lines) |
| Fortran binary | MAVEN mag L3, STATIC d1/c6 data |
| `.cdf` | MAVEN SWEA, SWIA, KP data |
| `.jld2` | Pre-processed NGIMS, KP L3 data |

## Plotting (CairoMakie)

- Headless backend — no display needed, use `save("output.png", fig)`
- Always activate: `CairoMakie.activate!()`
- Use LaTeXStrings for math labels: `L"B_{\mathrm{min}}"`
- Use `limits!()` and `hidespines!()` for clean axes
- Named parameters: `lines!(ax, x, y; color=:blue, linewidth=2)`

## Dependencies

Core packages (must be in global environment):
- `CairoMakie`, `GeometryBasics`, `LaTeXStrings` — plotting
- `TimesDates`, `Dates`, `DataFrames` — time/data handling
- `DSP`, `Wavelets`, `FFTW` — signal processing
- `LinearAlgebra`, `Statistics` — math
- `CommonDataFormat` — CDF file reading
- `FortranFiles` — Fortran binary reading
- `JLD2` — Julia data serialization
- `Rotations` — quaternion/rotation handling
- `Optim`, `SpecialFunctions` — wavelet significance testing
- `DataInterpolations`, `DelaunayTriangulation`, `ProgressMeter`, `ColorTypes` — MAVEN plotting

## Code Conventions

- Chinese comments for docstrings (matches existing codebase)
- `snake_case` for functions, `PascalCase` for modules/types
- Error handling: return `Dict(:data_load_flag => false)` rather than throwing
- Preallocate arrays with `Matrix{Float32}(undef, n, m)` + `@inbounds` loops
- Broadcasting for element-wise ops: `sqrt.(x)`

## Notes

- No test suite exists — validate with demo scripts
- Data files are gitignored; download from NASA/USTC portals (see README.md links)
- Scripts use relative paths: `joinpath(@__DIR__, "..", "Data", ...)`
