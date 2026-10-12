#===============================================================================================
Spherical spline interpolation based on:

Blick, C. & Freeden, W. (2011). Spherical Spline Application to Radio Occultation Data. Journal of Geodetic Science, 1(4), pp. 379-395. Retrieved 20 Nov. 2018, from doi:10.2478/v10156-011-0017-8.

Blick, C. (2012). Spherical Spline Approximation and its Application to Radio Occultation Data for Climate Monitoring. 2nd International Radio Occultation Workshop 2012.

Usage

Given direction vectors stored as a 3×N matrix (each column is a unit vector) or an N×3 matrix,
and a vector of real values associated with those directions, return an interpolating spline
solution (which passes through all the values in the given directions).

S = spherical_spline(directions, values)
y = S(new_direction)

All direction vectors must be normalized to unit vectors.
===============================================================================================#
module SphericalSplinesOptim
export spherical_spline

using LinearAlgebra
using Statistics
using StaticArrays

"""
    SphericalSpline

A callable type representing a spherical spline interpolation solution.
"""
struct SphericalSpline{T<:AbstractFloat}
    c::T
    a::Vector{T}
    directions::Vector{SVector{3,T}}  # internal storage: always Vector of SVectors
end

"""
    (S::SphericalSpline)(direction::AbstractVector)

Evaluate the spline at a given 3D direction vector.
"""
function (S::SphericalSpline)(direction::AbstractVector)
    @boundscheck length(direction) == 3 || throw(DimensionMismatch("Direction must be a 3D vector"))
    result = S.c
    @inbounds for k in eachindex(S.a)
        result -= S.a[k] * Green(direction, S.directions[k])
    end
    return result
end

"""
    spherical_spline(directions, values)

Construct a spherical spline interpolant. `directions` can be either:
- a 3×N matrix (each column is a direction vector)
- an N×3 matrix (each row is a direction vector)

`values` must be a vector of length N.
"""
function spherical_spline(directions::AbstractMatrix{T}, values::AbstractVector{T}) where {T<:AbstractFloat}
    # Normalize to internal format: Vector{SVector{3,T}}
    sz = size(directions)
    if sz[1] == 3 && sz[2] != 3
        # Assume 3×N format
        dirs = [SVector{3,T}(directions[:, j]) for j in 1:sz[2]]
    elseif sz[2] == 3
        # Assume N×3 format
        dirs = [SVector{3,T}(directions[i, :]) for i in 1:sz[1]]
    else
        throw(ArgumentError("Input directions must be either 3×N or N×3, got $(size(directions))"))
    end

    N = length(dirs)
    length(values) == N || throw(DimensionMismatch("Number of values ($(length(values))) does not match number of directions ($N)"))

    if N ≤ 3
        return SphericalSpline{T}(mean(values), T[], dirs)
    end

    # Build Gram matrix (symmetric, diagonal = 1/(4π))
    G = zeros(T, N, N)
    inv4pi = one(T) / (T(4) * T(π))
    @inbounds for n in 1:N
        an = dirs[n]
        G[n, n] = inv4pi
        for m in (n+1):N
            g = Green(an, dirs[m])
            G[n, m] = g
            G[m, n] = g
        end
    end

    # Solve system: G * a = c*ones(N) - values
    F = cholesky!(Symmetric(G), RowMaximum(); check = true)

    # Compute u = G^{-1} * values, v = G^{-1} * ones(N)
    u = F \ values
    v = F \ ones(T, N)

    c = sum(u) / sum(v)
    a = c .* v .- u

    return SphericalSpline{T}(c, a, dirs)
end

"""
    Green(a, b)

Green's function for the iterated Beltrami operator on the sphere.
Both arguments are 3D direction vectors.
"""
@inline function Green(a::AbstractVector{T}, b::AbstractVector{T}) where {T<:AbstractFloat}
    d = dot(a, b)
    inv4pi = one(T) / (T(4) * T(π))
    pi24 = T(π) / T(24)
    pi2_12 = T(π)^2 / T(12)
    l2 = log(T(2))

    if d > one(T) - eps(T)
        return inv4pi
    elseif d < -one(T) + eps(T)
        return inv4pi - pi24
    end

    if abs(d) < eps(T)
        return (one(T) - pi2_12 - T(0.5) * l2^2) * inv4pi
    end

    t1 = log1p(-d) * (log1p(d) - l2)
    t2 = Li2((one(T) - d) / T(2))
    t3 = l2^2 - l2 * log1p(d)

    return (one(T) - t1 - t2 - t3) * inv4pi
end

"""
    Li2(x)

Compute the dilogarithm (Spence's function) to machine precision.
Uses the series expansion with reflection formula for rapid convergence.
"""
function Li2(x::T; tol = eps(T)) where {T<:AbstractFloat}
    if x > T(0.5)
        return T(π)^2 / T(6) - log(x) * log1p(-x) - Li2(one(T) - x)
    end

    s = zero(T)
    term = x
    k = 1
    while true
        s += term
        ratio = T(k) / T(k + 1)
        term *= x * ratio * ratio
        k += 1
        if abs(term) <= tol * abs(s)
            break
        end
    end
    return s
end

end # module