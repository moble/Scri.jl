# Quadrature on the sphere, for objectives that should be true L² norms.
#
# TODO: Gauss–Legendre quadrature belongs in SphericalFunctions; move it there and use it
# from here once it exists.  Be careful about the numerics, though, because it (and how it
# interacts with sphere_quadrature) can go bad near the poles or the equator depending on
# the method.  Maybe just look into FastGaussQuadrature.jl (though it's just for Float64,
# IIRC) or similar. (2026-10-07)

@doc raw"""
    gauss_legendre(n, T=Float64)

Return the nodes `x` (ascending) and weights `w` of the `n`-point Gauss–Legendre rule on
``[-1, 1]``, which integrates every polynomial of degree at most ``2n-1`` exactly:
``∫_{-1}^{1} p(x)\, dx = ∑_i w_i\, p(x_i)``.

The nodes are found by Newton's method on the Legendre polynomial ``P_n``, evaluated by its
three-term recurrence, starting from Tricomi's approximation ``\cos(π(4i-1)/(4n+2))``.  This
works in any floating-point type `T`, such as `Double64` or `BigFloat`, where an eigenvalue
method would need a generic eigensolver.  The weights are ``2/((1-x_i^2) P_n'(x_i)^2)``.
"""
function gauss_legendre(n::Int, ::Type{T}=Float64) where {T<:AbstractFloat}
    if n < 1
        throw(ArgumentError("Gauss–Legendre needs at least one node; got n=$n"))
    end
    x = Vector{T}(undef, n)
    w = Vector{T}(undef, n)
    for i ∈ 1:((n + 1) ÷ 2)
        z = cos(T(π) * (4i - 1) / (4n + 2))
        # Newton's method converges quadratically from this start; one more step after the
        # correction falls below roundoff makes the last digits as good as they can be.
        for _ ∈ 1:100
            p, dp = legendre_and_derivative(n, z)
            Δ = p / dp
            z -= Δ
            abs(Δ) ≤ 4eps(T) && break
        end
        p, dp = legendre_and_derivative(n, z)
        z -= p / dp
        _, dp = legendre_and_derivative(n, z)
        wᵢ = 2 / ((1 - z^2) * dp^2)
        x[i], x[n + 1 - i] = -z, z
        w[i], w[n + 1 - i] = wᵢ, wᵢ
    end
    return x, w
end

"""
    legendre_and_derivative(n, x)

Return ``P_n(x)`` and ``P_n'(x)``, from the recurrence ``k P_k = (2k-1) x P_{k-1} - (k-1)
P_{k-2}`` and the identity ``(x^2-1) P_n' = n (x P_n - P_{n-1})``, which holds away from
``x = ±1``, where Gauss–Legendre nodes never lie.
"""
function legendre_and_derivative(n::Int, x::T) where {T}
    p₀, p₁ = one(T), x
    if n == 0
        return p₀, zero(T)
    end
    # n == 1 intentionally drops through to the final statement
    for k ∈ 2:n
        p₀, p₁ = p₁, ((2k - 1) * x * p₁ - (k - 1) * p₀) / k
    end
    return p₁, n * (x * p₁ - p₀) / (x^2 - 1)
end

@doc raw"""
    sphere_quadrature(ℓₘₐₓ, T=Float64)

Return rotors `R` and weights `w` of a quadrature on the sphere that integrates exactly the
product of any two fields of the same spin weight that are band-limited to ``ℓ ≤ ℓₘₐₓ``:

```math
∫ f\, \bar{g}\, dΩ = ∑_i w_i\, f(R_i)\, \bar{g}(R_i) = ∑_{ℓ,m} f_{ℓ,m}\, \bar{g}_{ℓ,m}.
```

The grid is the product of ``ℓₘₐₓ+1`` Gauss–Legendre nodes in ``\cos θ``, from
[`gauss_legendre`](@ref), and ``2ℓₘₐₓ+1`` equally spaced values of ``φ``, about twice as
many points as there are modes.  That suffices because such a product is a polynomial of
degree at most ``2ℓₘₐₓ`` in ``\cos θ`` times ``e^{ikφ}`` with ``|k| ≤ 2ℓₘₐₓ``.  The rotors
are `from_spherical_coordinates(θ, φ)`; no node lies at a pole.
"""
function sphere_quadrature(ℓₘₐₓ::Integer, ::Type{T}=Float64) where {T<:AbstractFloat}
    if ℓₘₐₓ < 0
        throw(ArgumentError("ℓₘₐₓ must be nonnegative; got ℓₘₐₓ=$ℓₘₐₓ"))
    end
    x, wₓ = gauss_legendre(ℓₘₐₓ + 1, T)
    Nˣ = length(x)
    Nᵠ = 2ℓₘₐₓ + 1
    R = Vector{Rotor{T}}(undef, Nˣ * Nᵠ)
    w = Vector{T}(undef, Nˣ * Nᵠ)
    let π = T(π)
        for (j, (xⱼ, wⱼ)) ∈ enumerate(zip(x, wₓ))
            θⱼ, w̃ⱼ = acos(xⱼ), wⱼ * 2π / Nᵠ    # once per ring
            for k ∈ 0:(Nᵠ - 1)
                i = Nᵠ * (j - 1) + k + 1
                R[i] = from_spherical_coordinates(θⱼ, 2π * k / Nᵠ)
                w[i] = w̃ⱼ
            end
        end
    end
    return R, w
end
