"""
    PhysicsCoordTransforms

A zero-dependency module providing coordinate system transformations
following the ISO 80000-2 (physics) convention.

## Conventions

### 2D Polar
- `(r, ϕ)` where `r ≥ 0`, `ϕ ∈ (-π, π]` is the angle from the +x axis.

### 3D Spherical
- `(r, θ, ϕ)` where:
  - `r ≥ 0` is the radial distance,
  - `θ ∈ [0, π]` is the polar angle (colatitude) from the +z axis,
  - `ϕ ∈ (-π, π]` is the azimuthal angle from the +x axis toward +y.

### 3D Cylindrical
- `(ρ, ϕ, z)` where:
  - `ρ ≥ 0` is the radial distance in the xy-plane,
  - `ϕ ∈ (-π, π]` is the azimuthal angle from the +x axis,
  - `z` is the height along the +z axis.
"""
module PhysicsCoordTransforms

export
    # 2D
    polar_from_cartesian,
    cartesian_from_polar,
    # 3D
    spherical_from_cartesian,
    cartesian_from_spherical,
    spherical_from_cylindrical,
    cylindrical_from_spherical,
    cartesian_from_cylindrical,
    cylindrical_from_cartesian

# ============================================================================
# Helper: safe type promotion and conversion
# ============================================================================
@inline _clamp_cos(x::T) where {T<:AbstractFloat} = clamp(x, -one(T), one(T))

# ============================================================================
# 2D Transformations
# ============================================================================

"""
    polar_from_cartesian()

Singleton transformation: Cartesian `(x, y)` → Polar `(r, ϕ)`.
Returns a named tuple `(r=, ϕ=)`.
"""
struct polar_from_cartesian end

function (::polar_from_cartesian)(x::Real, y::Real)
    xₜ, yₜ = promote(x, y)
    r = hypot(xₜ, yₜ)
    ϕ = atan(yₜ, xₜ)
    return (r=r, ϕ=ϕ)
end

"""
    cartesian_from_polar()

Singleton transformation: Polar `(r, ϕ)` → Cartesian `(x, y)`.
Returns a tuple `(x, y)`.
"""
struct cartesian_from_polar end

function (::cartesian_from_polar)(r::Real, ϕ::Real)
    rₜ, ϕₜ = promote(r, ϕ)
    c, s = sincos(ϕₜ)
    return (rₜ * c, rₜ * s)
end

# ============================================================================
# 3D: Spherical ↔ Cartesian
# ============================================================================

"""
    spherical_from_cartesian()

Singleton transformation: Cartesian `(x, y, z)` → Spherical `(r, θ, ϕ)`.
Physics convention: θ = polar angle from +z, ϕ = azimuth from +x.
Returns a named tuple `(r=, θ=, ϕ=)`.
"""
struct spherical_from_cartesian end

function (::spherical_from_cartesian)(x::Real, y::Real, z::Real)
    xₜ, yₜ, zₜ = promote(x, y, z)
    r = sqrt(xₜ^2 + yₜ^2 + zₜ^2)
    if iszero(r)
        nan = convert(typeof(r), NaN)
        return (r=nan, θ=nan, ϕ=nan)
    end
    θ = acos(_clamp_cos(zₜ / r))
    ϕ = atan(yₜ, xₜ)
    return (r=r, θ=θ, ϕ=ϕ)
end

"""
    cartesian_from_spherical()

Singleton transformation: Spherical `(r, θ, ϕ)` → Cartesian `(x, y, z)`.
Physics convention. Returns a tuple `(x, y, z)`.
"""
struct cartesian_from_spherical end

function (::cartesian_from_spherical)(r::Real, θ::Real, ϕ::Real)
    rₜ, θₜ, ϕₜ = promote(r, θ, ϕ)
    sinθ, cosθ = sincos(θₜ)
    sinϕ, cosϕ = sincos(ϕₜ)
    return (rₜ * sinθ * cosϕ,
            rₜ * sinθ * sinϕ,
            rₜ * cosθ)
end

# ============================================================================
# 3D: Cylindrical ↔ Cartesian
# ============================================================================

"""
    cylindrical_from_cartesian()

Singleton transformation: Cartesian `(x, y, z)` → Cylindrical `(ρ, ϕ, z)`.
Returns a named tuple `(ρ=, ϕ=, z=)`.
"""
struct cylindrical_from_cartesian end

function (::cylindrical_from_cartesian)(x::Real, y::Real, z::Real)
    xₜ, yₜ, zₜ = promote(x, y, z)
    ρ = hypot(xₜ, yₜ)
    ϕ = atan(yₜ, xₜ)
    return (ρ=ρ, ϕ=ϕ, z=zₜ)
end

"""
    cartesian_from_cylindrical()

Singleton transformation: Cylindrical `(ρ, ϕ, z)` → Cartesian `(x, y, z)`.
Returns a tuple `(x, y, z)`.
"""
struct cartesian_from_cylindrical end

function (::cartesian_from_cylindrical)(ρ::Real, ϕ::Real, z::Real)
    ρₜ, ϕₜ, zₜ = promote(ρ, ϕ, z)
    sinϕ, cosϕ = sincos(ϕₜ)
    return (ρₜ * cosϕ, ρₜ * sinϕ, zₜ)
end

# ============================================================================
# 3D: Spherical ↔ Cylindrical
# ============================================================================

"""
    spherical_from_cylindrical()

Singleton transformation: Cylindrical `(ρ, ϕ, z)` → Spherical `(r, θ, ϕ)`.
Note: the azimuthal angle ϕ is shared between both systems.
Returns a named tuple `(r=, θ=, ϕ=)`.
"""
struct spherical_from_cylindrical end

function (::spherical_from_cylindrical)(ρ::Real, ϕ_cyl::Real, z::Real)
    ρₜ, ϕₜ, zₜ = promote(ρ, ϕ_cyl, z)
    r = hypot(ρₜ, zₜ)
    if iszero(r)
        nan = convert(typeof(r), NaN)
        return (r=nan, θ=nan, ϕ=nan)
    end
    θ = acos(_clamp_cos(zₜ / r))
    return (r=r, θ=θ, ϕ=ϕₜ)
end

"""
    cylindrical_from_spherical()

Singleton transformation: Spherical `(r, θ, ϕ)` → Cylindrical `(ρ, ϕ, z)`.
Note: the azimuthal angle ϕ is shared between both systems.
Returns a named tuple `(ρ=, ϕ=, z=)`.
"""
struct cylindrical_from_spherical end

function (::cylindrical_from_spherical)(r::Real, θ::Real, ϕ_sph::Real)
    rₜ, θₜ, ϕₜ = promote(r, θ, ϕ_sph)
    sinθ, cosθ = sincos(θₜ)
    return (ρ=rₜ * sinθ, ϕ=ϕₜ, z=rₜ * cosθ)
end

end # module PhysicsCoordTransforms