"""
    CubicSplineCache{T<:AbstractFloat}

Precomputed LU factorization for natural cubic spline interpolation on a fixed
time grid.  Stores the interval widths, their reciprocals, and the Thomas
factors for the tridiagonal second-derivative system.

After constructing the cache once from the knot vector `t`, repeated
interpolation of new data vectors `dₖ` requires only an O(N) forward/backward
sweep to compute the second derivatives `d̈ₖ`, plus O(1) per evaluation point.

# Fields
- `h`:   interval widths `h[j] = t[j+1]−t[j]`, length N−1
- `h⁻¹`: reciprocals `1/h[j]`, length N−1
- `l`:   Thomas sub-diagonal factors `lᵢ = h[i+1]/uᵢ` for i=1..N−3
- `u⁻¹`: reciprocals of the modified diagonal `1/uᵢ` for i=1..N−2
"""
struct CubicSplineCache{T<:Real}
    h::Vector{T}
    h⁻¹::Vector{T}
    l::Vector{T}
    u⁻¹::Vector{T}
end

"""
    CubicSplineCache(t)

Precompute the Thomas LU factorization of the natural-cubic-spline tridiagonal
system for knot vector `t`.  `t` must be strictly increasing with at least 4
elements.

# Natural spline system

For N knots, the N−2 interior second derivatives `d̈[2..N-1]` satisfy the
symmetric tridiagonal system

    h[j-1]*d̈[j-1] + 2(h[j-1]+h[j])*d̈[j] + h[j]*d̈[j+1] = rₖ[j]

with `d̈[1]=d̈[N]=0` (natural boundary conditions) and

    rₖ[j] = 6*( (d[j+1]−d[j])/h[j] − (d[j]−d[j-1])/h[j-1] ).

Denoting the system `Ad̈ₖ = rₖ`, the Thomas LU factorization gives

    A = LU,  L unit lower bidiagonal with factors lᵢ,
              U upper bidiagonal with diagonal uᵢ and super-diagonal sᵢ = h[i+1].

# Thomas forward sweep

With the cache built, the forward sweep for a new data vector d is

    z[2] = r[2] * u⁻¹[1]
    z[j] = ( r[j] − h[j-1] * z[j-1] ) * u⁻¹[j-1]   for j = 3..N-1

and the backward sweep recovers d̈ from z:

    d̈[N-1] = z[N-1]
    d̈[j]   = z[j] − l[j-1] * d̈[j+1]                for j = N-2:-1:2

# Evaluation

On interval [t[j], t[j+1]] at offset τ = t_query − t[j]:

    c = (d̈[j+1]−d̈[j]) * h⁻¹[j] / 6
    b = d̈[j] / 2
    a = h⁻¹[j]*(d[j+1]−d[j]) − h[j]/6*(2d̈[j]+d̈[j+1])
    S = d[j] + τ*(a + τ*(b + τ*c))
"""
function CubicSplineCache(t::AbstractVector{T}) where {T<:Real}
    N = length(t)
    N ≥ 4 || throw(ArgumentError("CubicSplineCache requires at least 4 knots, got $N"))
    h = Vector{T}(undef, N - 1)
    for j ∈ 1:(N - 1)
        h[j] = t[j + 1] - t[j]
        h[j] > 0 || throw(
            ArgumentError(
                "Knot vector must be strictly increasing: t[$j]=$(t[j]) ≥ t[$(j+1)]=$(t[j+1])",
            ),
        )
    end
    h⁻¹ = inv.(h)
    u⁻¹ = Vector{T}(undef, N - 2)
    l = Vector{T}(undef, N - 3)
    u⁻¹[1] = inv(2 * (h[1] + h[2]))
    for k ∈ 2:(N - 2)
        l[k - 1] = h[k] * u⁻¹[k - 1]
        u⁻¹[k] = inv(2 * (h[k] + h[k + 1]) - l[k - 1] * h[k])
    end
    return CubicSplineCache{T}(h, h⁻¹, l, u⁻¹)
end

"""
    spline_eval(dⱼ, dⱼ₊₁, d̈ⱼ, d̈ⱼ₊₁, hⱼ, h⁻¹ⱼ, τ)

Evaluate the natural cubic spline on interval j at offset τ = t_query − t[j],
given the knot values `dⱼ`, `dⱼ₊₁`, the second derivatives `d̈ⱼ`, `d̈ⱼ₊₁`,
the interval width `hⱼ`, and its reciprocal `h⁻¹ⱼ`.  Horner form.
"""
@inline function spline_eval(dⱼ, dⱼ₊₁, d̈ⱼ, d̈ⱼ₊₁, hⱼ, h⁻¹ⱼ, τ)
    c = (d̈ⱼ₊₁ - d̈ⱼ) * h⁻¹ⱼ / 6
    b = d̈ⱼ / 2
    a = h⁻¹ⱼ * (dⱼ₊₁ - dⱼ) - hⱼ / 6 * (2 * d̈ⱼ + d̈ⱼ₊₁)
    return dⱼ + τ * (a + τ * (b + τ * c))
end
