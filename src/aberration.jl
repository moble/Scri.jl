@doc raw"""
    aberration(RRₚᵢ, v⃗; emitted=true)

Transform from the rotor `RRₚᵢ` in the boosted frame to the corresponding rotor in the rest
frame, given the boost velocity vector `v⃗`.  This is derived in Appendix C of
[Boyle_2015](@cite).  It is also described in Penrose-Rindler Vol. 1, around Eq. (1.3.5).

A rotor on the sphere encodes not only a point (direction) but also an orientation of the
tangent frame at that point — the full bundle on which spin-weighted functions are defined.
This function implements the complete transformation of that bundle under a boost: both the
direction aberration and the rotation of the tangent frame (the spin-weight phase factor).

The input `RRₚᵢ = R * Rₚᵢ` is the product of the overall frame rotation `R` and the pixel
rotor `Rₚᵢ`.  The rotation is applied **before** the boost correction, matching the standard
BMS transformation order (translate → rotate → boost).  The boost correction `B′` is then
computed from the rotated direction `n̂′ = RRₚᵢ(𝐤)`.

All rotors — input and output — are expressed in the same standard spatial coordinate basis;
the "boosted frame" label refers to which physical direction the input rotor points to, not
to any change of coordinates.

## The `emitted` keyword

The `emitted` keyword selects which null infinity the field lives on:

- `emitted = true` (default): **future null infinity** ℐ⁺ — outgoing radiation emitted by
  the source.  A wave propagating at angle `Θ` in the rest frame appears at the *larger*
  angle `Θ̑ > Θ` in the boosted frame.
- `emitted = false`: **past null infinity** ℐ⁻ — incoming radiation received by the
  observer.  A wave arriving from direction `n̂` in the rest frame appears shifted *toward*
  the boost axis in the boosted frame, so `Θ̑ < Θ`.

Penrose-Rindler work on the past celestial sphere (`emitted = false`), which reverses the
sign compared to the `emitted = true` convention used throughout this package; hence the
sign difference noted at Eq. (1.3.5) there.

## Formula

The function returns `B′ * RRₚᵢ`, where `B′` is the aberration rotor.  `B′` is computed as
`exp(n̂′ × v⃗ ⋅ δ/2)` where `δ = (Θ̑−Θ)/|n̂′ × v⃗|` absorbs the normalization.  Here `Θ̑` is
the polar angle of `n̂′` relative to `v⃗` in the boosted frame and `Θ` is the corresponding
angle in the rest frame, related by

```math
\tan(Θ/2) = e^{-εφ} \tan(Θ̑/2),
\quad φ = \operatorname{atanh}(β),
\quad ε = \begin{cases}+1 & \text{emitted} \\ -1 & \text{received}\end{cases}
```

`Θ̑` is computed via `atan(|n̂′ × v⃗|, v⃗ ⋅ n̂′)` (a scale-invariant `atan2` that requires
no division by `|v⃗|`), keeping the computation well-conditioned.

## Numerical stability

When `β sinΘ̑ = |n̂′ × v⃗|` is small — either because `β → 0` or because `n̂′` is nearly
(anti-)parallel to `v⃗` — the exact formula is ill-conditioned.  In that regime we use
```math
\exp\!\left(\frac{n̂′ \times v⃗}{|n̂′ × v⃗|}\,\frac{Θ̑ - Θ}{2}\right)
= \cos\!\left(\frac{Θ̑ - Θ}{2}\right)
+ \frac{n̂′ \times v⃗}{|n̂′ × v⃗|}\,\sin\!\left(\frac{Θ̑ - Θ}{2}\right),
```
and expand the `cos` term and the ratio ``\sin((Θ̑-Θ)/2)\,/\,(β \sin Θ̑)`` as Taylor series
in `εβ` to fifth order.  Because the error in a degree-5 truncation is ``O((β \sin Θ̑)^6)``,
the Taylor branch is used when ``β \sin Θ̑ < \epsilon^{1/3}``; the resulting error is
``O(\epsilon^2)``, well below floating-point precision.

The implementation is compatible with arbitrary floating-point types, including dual numbers
for automatic differentiation.  Branch selection is based on `value(β sinΘ̑)` so that the
derivative of the chosen branch is always evaluated, with no discontinuity in derivatives at
the branch boundary.
"""
function aberration(RRₚᵢ, v⃗; emitted::Bool=true)
    # Equation numbers refer to Appendix C of Boyle (2015).
    β = absvec(v⃗)
    φ = atanh(β)
    ε = emitted ? 1 : -1  # +1 for ℐ⁺ (outgoing), -1 for ℐ⁻ (incoming)
    εβ = ε * β            # signed β; ε = -1 flips the rapidity φ → -φ
    n̂′ = RRₚᵢ(𝐤)  # direction in the boosted frame corresponding to this pixel

    # The quaternion product of two pure vectors p, q satisfies pq = -(p⋅q) + p×q,
    # giving the dot product in the scalar part and the cross product in the vector part.
    n̂′v⃗ = n̂′ * v⃗
    n̂′xv⃗ = QuatVec(n̂′v⃗)  # extracts the vector (cross-product) part; zeros the scalar part

    μ = absvec(n̂′xv⃗)        # |n̂′ × v⃗| = β sinΘ̑; vanishes when n̂′ ∥ ±v⃗
    # atan(β sinΘ̑, β cosΘ̑) = Θ̑; the β factors cancel so this is independent of |v⃗|
    Θ̑ = atan(μ, -n̂′v⃗.w)    # Eq. (C6): n̂′v⃗.w = -(n̂′⋅v⃗) = -β cosΘ̑
    sinΘ̑, cosΘ̑ = sincos(Θ̑)
    T = typeof(sinΘ̑)  # float type, which may be a Dual for AD
    ϵ = value(eps(T))

    # Use the Taylor expansion when β sinΘ̑ is too small for the exact formula.
    # The 5th-order expansion has error O((β sinΘ̑)^6), so the threshold ∛ϵ gives O(ϵ²).
    B′ = if value(β * sinΘ̑) < ∛ϵ
        # Taylor expansion of cos((Θ̑−Θ)/2) about β=0, evaluated at εβ.
        # Flipping ε negates φ (i.e., β → -β), so even powers are unchanged and odd
        # powers pick up a sign; using εβ in the Horner form handles this automatically.
        cosΔΘ╱2 = (
            1 +
            εβ * (
                0 +
                εβ * (
                    -sinΘ̑^2 / 8 +
                    εβ * (
                        sinΘ̑^2 * cosΘ̑ / 8 +
                        εβ * (
                            5 * (3sinΘ̑^2 - 4) * sinΘ̑^2 / 128 +
                            εβ * (10 - 7sinΘ̑^2) * sinΘ̑^2 * cosΘ̑ / 64
                        )
                    )
                )
            )
        )

        # Taylor expansion of sin((Θ̑−Θ)/2) / (β sinΘ̑) about β=0.
        # For ε = -1 the ratio is negative (leading term -1/2); factoring out ε and using
        # εβ inside restores a uniform +1/2 leading term in the series body.
        # Multiplying by n̂′xv⃗ (which carries the explicit factor of β sinΘ̑) gives the
        # full vector term sin((Θ̑−Θ)/2) * (n̂′×v⃗)/|n̂′×v⃗| without any division by μ.
        sinΔΘ╱2╱βsinΘ̑ =
            ε * (
                1 +
                εβ * (
                    -cosΘ̑ / 2 +
                    εβ * (
                        (1 - 3sinΘ̑^2 / 4) / 2 +
                        εβ * (
                            (-5cosΘ̑^2 - 1) * cosΘ̑ / 16 +
                            εβ * (
                                (35sinΘ̑^4 / 16 - 19sinΘ̑^2 / 4 + 3) / 8 +
                                εβ * (-61cosΘ̑^4 + 34cosΘ̑^2 + 27) * cosΘ̑ / 768
                            )
                        )
                    )
                )
            ) / 2

        # exp[n̂′×v⃗/|n̂′×v⃗| * (Θ̑−Θ)/2] to fifth order in |n̂′xv⃗|
        Q = cosΔΘ╱2 + n̂′xv⃗ * sinΔΘ╱2╱βsinΘ̑
        Rotor{basetype(Q)}(Q)
    else
        Θ = 2atan(exp(-ε * φ) * tan(Θ̑ / 2))  # Eq. (C7): rest-frame polar angle; ε flips sign
        exp((n̂′xv⃗/μ) * (Θ̑-Θ)/2)              # Eq. (C8): sign of (Θ̑-Θ) carries ε
    end

    return B′ * RRₚᵢ
end
