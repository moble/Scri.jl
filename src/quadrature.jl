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
    one_minus_node(n, z)

Return ``1 - z`` for a node ``z ≥ 0`` of the `n`-point Gauss–Legendre rule, to full relative
accuracy.  A node near ``z = 1`` is known only to an absolute accuracy of about ``ε``, so
the difference ``1 - z`` there has a relative error of about ``ε/(1 - z)``, which is
``ε/θ^2`` in terms of the colatitude ``θ`` of the node.  As [Bogaert_2014](@citet) points
out, the nodes cluster quadratically near ``z = ±1``, so they are better computed in ``θ``,
or in ``1 - z``, than in ``z`` itself.  For ``z ≥ 1/2``, ``v = 1 - z`` is therefore found
afresh by Newton's method in ``v``, starting from ``1 - z``.  The Legendre
polynomial is evaluated there by Reinsch's form of the recurrence, which works with the
differences ``d_k = P_k - P_{k-1}``,

```math
d_k = \frac{(k-1)\, d_{k-1} - (2k-1)\, v\, P_{k-1}}{k},
\qquad
P_k = P_{k-1} + d_k,
```

starting from ``P_0 = 1`` and ``d_1 = -v``, so that nothing is lost to cancellation near
``z = 1``.  The derivative follows from the identity of
[`legendre_and_derivative`](@ref Scri.legendre_and_derivative), as ``dP_n/dv = n (d_n - v
P_n)/(v (2 - v))``.  For ``z < 1/2``, ``1 - z`` is already accurate and is returned as it is.
"""
function one_minus_node(n::Int, z::T) where {T}
    z < 1 / 2 && return 1 - z
    v = 1 - z  # exact, by Sterbenz's lemma, but no more accurate than z
    # Newton's method converges quadratically from an error of about ε, so two steps
    # suffice; the third makes the last digits as good as they can be.
    for _ ∈ 1:3
        d = -v              # d₁ = P₁ - P₀
        p = one(T) + d      # P₁
        for k ∈ 2:n
            d = ((k - 1) * d - (2k - 1) * v * p) / k
            p += d
        end
        v -= p * v * (2 - v) / (n * (d - v * p))
    end
    return v
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
are those of `from_spherical_coordinates(θ, φ)`; no node lies at a pole.

The rotors are formed from ``\sin(θ/2)`` and ``\cos(θ/2)``, which are ``\sqrt{v/2}`` and
``\sqrt{1 - v/2}`` for a node with ``\cos θ = 1 - v``, rather than from ``θ = \arccos x``.
Near a pole, ``\arccos x`` has a relative error of about ``ε/θ^2``, and the small component
of a rotor there would inherit it: about ``26ε`` at ``ℓₘₐₓ = 32`` and ``24000ε`` at ``ℓₘₐₓ =
1000``.  With ``v`` from [`one_minus_node`](@ref Scri.one_minus_node), and the roles of the
two reversed at the southern nodes, the relative error is about ``ε`` at every node.  The
same holds in azimuth, where ``\sin(φ/2)`` and ``\cos(φ/2)`` are formed with `sinpi` from
exact ratios of integers no larger than ``1/2``, because either of them formed from a
rounded ``φ`` would lose digits where it is small.
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
            # Once per ring: the sine and cosine of half the colatitude, from v = 1 - |xⱼ|,
            # which is 1 - cos θ for a northern node and 1 + cos θ for a southern one.
            v = one_minus_node(Nˣ, abs(xⱼ))
            a, b = sqrt(v / 2), sqrt(1 - v / 2)
            sθ, cθ = xⱼ ≥ 0 ? (a, b) : (b, a)
            w̃ⱼ = wⱼ * 2π / Nᵠ
            for k ∈ 0:(Nᵠ - 1)
                i = Nᵠ * (j - 1) + k + 1
                # sin(φ/2) and cos(φ/2) at φ = 2πk/Nᵠ, each as sinpi of an exact ratio of
                # integers no larger than 1/2, so that neither loses digits where it is small.
                sφ = sinpi(T(min(k, Nᵠ - k)) / Nᵠ)
                cφ = sinpi(T(Nᵠ - 2k) / 2Nᵠ)
                # The components of `from_spherical_coordinates(θ, φ)`
                R[i] = Quaternionic.rotor(cθ * cφ, -sθ * sφ, sθ * cφ, cθ * sφ)
                w[i] = w̃ⱼ
            end
        end
    end
    return R, w
end
