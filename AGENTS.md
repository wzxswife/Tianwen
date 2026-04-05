# Tianwen Repository Agent Guidelines

This repository contains Julia code for space physics data analysis, focusing on MAVEN spacecraft and Tianwen-1 mission magnetic field and plasma data.

## Project Structure

```
Tianwen/
├── src/                      # Main entry points and demo scripts
│   ├── scripts/              # Reusable modules (load, plot, calculate)
│   │   ├── MAVEN_load.jl     # MAVEN data loading module
│   │   ├── TW_load.jl        # Tianwen data loading module
│   │   ├── MAVEN_plot.jl     # MAVEN plotting functions
│   │   ├── TW_plot.jl        # Tianwen plotting functions
│   │   ├── MVA_plot.jl       # Minimum Variance Analysis plotting
│   │   ├── wave_caculate.jl  # Wave analysis calculations
│   │   └── TCWavelet.jl       # Wavelet transform
│   ├── DemoMAVEN.jl           # MAVEN demo script
│   └── DemoTW.jl              # Tianwen demo script
├── Project/                   # Sub-projects
│   ├── IonBeam/              # Ion beam analysis
│   └── TianwenData/          # Tianwen data processing
├── data/                      # Data files (gitignored)
└── Results/                   # Output results
```

## Running Scripts

```bash
# Run a demo script directly
julia src/DemoMAVEN.jl

# Run with specific data path
julia --project=. src/scripts/MVA_plot.jl
```

## Code Style Guidelines

### Module Organization
- Use `module` to define namespaces (e.g., `MAVEN_load`, `TW_plot`)
- Use `export` to define public API
- Use `include()` for loading dependent modules
- Use `.ModuleName` syntax when accessing submodules

```julia
module MAVEN_load
using Dates, DataFrames
export load_mag_l2, load_cdf

include("helper_functions.jl")
# ... module content ...
end
```

### Data Structures
- Use `Dict{Symbol, Any}()` for data containers (consistent with CDF/spacescraft data formats)
- Access fields with Symbol keys: `data[:epoch]`, `data[:B]`
- Include `:data_load_flag` to indicate successful loading

```julia
data = Dict{Symbol,Any}(
    :epoch => times,
    :B => B,
    :position => position,
    :data_load_flag => true
)
```

### Function Naming
- Use snake_case for function names: `load_mag_l2`, `plot_wavelet`
- Use PascalCase for types and modules
- Be descriptive but concise

### Type Annotations
- Annotate function parameters with types when beneficial:
```julia
function load_mag_l2(file::String)
function plot_MVA(ax1, ax2, ax3, mag_wave::Matrix{Float64}; smooth_window::Int=10)
```
- Use parametric types when appropriate: `Vector{T}` where T is a type parameter

### Imports
- Standard library imports first, then external packages
- Import only what's needed to avoid namespace pollution

```julia
using Dates
using Statistics
using LinearAlgebra
using CairoMakie, GeometryBasics
```

### Constants
- Define physical constants at module level in UPPER_SNAKE_CASE:
```julia
const Rm = 3390.0  # Mars radius in km
const Me = 9.109e-31  # Electron mass in kg
```

### Plotting (CairoMakie)
- Activate CairoMakie before plotting: `CairoMakie.activate!()`
- Use named parameters for clarity:
```julia
lines!(ax, x, y; color=:blue, linewidth=2, label="data")
scatter!(ax, x, y; markersize=10, color=:red)
```
- Use LaTeXStrings for mathematical labels: `xlabel = L"B_{min}"`
- Use `limits!()` to set axis limits, `hidespines!()` for clean axes

### Error Handling
- Return `Dict(:data_load_flag => false)` on load failure rather than throwing
- Use try-catch for file operations:
```julia
try
    cdf_data = CDFDataset(file)
catch e
    println("Error loading: ", file)
    return Dict(:data_load_flag => false)
end
```

### Documentation
- Use Chinese comments for function descriptions (matches existing codebase)
- Document input parameters and return values
```julia
"""
    load_mag_l2(file::String)
加载MAVEN磁力计L2数据
输入: file - 文件路径
返回: Dict包含:epoch, B, position, data_load_flag
"""
```

### Code Formatting
- 4 spaces for indentation (Julia standard)
- No trailing whitespace
- Maximum line length: 120 characters
- Use `begin`/`end` blocks for multi-line expressions when needed
- Use broadcasting (`.`) for element-wise operations: `sqrt.(x)`

## Common Patterns

### Time Range Selection
```julia
time_range = DateTime(date, Time(7, 30, 00)) .+ Dates.Second(10) .* range(0, 4)
```

### Finding Available Data
```julia
mag_data = find_avail_data(data, time_range, [:epoch, :B])
```

### Preallocating Arrays
```julia
B = Matrix{Float32}(undef, nums, 3)
@inbounds for (i, line) in enumerate(lines)
    # process line
end
```

### Unit Conversions
- Position typically normalized to Mars radii (Rm = 3390 km)
- Magnetic field in nT
- Energy in eV
- Use defined constants: `const Rm = 3390.0`

## Dependencies

Key packages used in this codebase:
- `CairoMakie` - Publication-quality plotting
- `GeometryBasics` - Geometric data structures
- `LaTeXStrings` - LaTeX formatting in plots
- `Dates` - Date/time handling
- `DataFrames` - Tabular data
- `DSP`, `Wavelets` - Signal processing
- `LinearAlgebra`, `Statistics` - Math operations
- `CommonDataFormat` - CDF file handling
- `FortranFiles` - Reading FORTRAN binary data

## Notes

- No formal test suite exists; validate changes with demo scripts
- Data files are gitignored; download from NASA/USTC data portals
- Scripts use relative paths with `joinpath(@__DIR__, "..", "data", ...)`
