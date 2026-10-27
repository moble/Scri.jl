"""
    aberration(R′ₚᵢ, Λ)

Transform from the pixel rotor `R′ₚᵢ` in the boosted frame to the corresponding rotor in the
rest frame, given the precomputed Lorentz rotor `Λ = Boost(ℐ*v⃗) * R`, which already folds in
the boost velocity `v⃗`, the overall frame rotation `R`, and the choice of null infinity
through `ℐ = ±1` — see [Future and past null infinity](@ref scri_pm_conventions).

The result is the ``K`` factor of the Iwasawa ``KAN`` decomposition of the full Lorentz
transformation ``Λ`` `* R′ₚᵢ`: it points along the same null direction for the rest observer
as that transformation produces for the boosted observer, and carries the tangent-frame
orientation (spin phase) along with it — exactly the information a spin-weighted field
needs.  The ``A`` factor is accounted for by the conformal factor ``κ``, and the ``N``
factor by component mixing, so only ``K`` is computed here, which is cheaper than the full
[`Quaternionic.KAN`](@extref Quaternionic :jl:function:`Quaternionic.KAN`): projecting with
the idempotent ``u₊ = (1+𝐭𝐳)/2`` gives ``ℂℜ(Λu₊) = ½e^{φₐ/2}𝐑_K``, which normalizes to
``𝐑_K``.  See the "Reversing the order" section of [The Lorentz Group](@ref) page for why
``K`` is the right object, and [the `Quaternionic` documentation](@extref Quaternionic
:std:label:`iwasawa-kan`) for the extraction.

The sign `ℐ` enters only through the boost velocity baked into `Λ` by the caller.  The
``u₊`` extraction is tied to the null vector ``ℓ ∝ (1, +k̂)`` where the pixel rotor takes
``𝐳 → k̂``, which matches the ℐ⁺ labeling directly.  At ℐ⁻ the labeling is antipodal — pixel
``k̂`` labels the null vector ``(1, −k̂)`` — and conjugating by spatial parity turns
``B(v⃗)`` acting on ``(1, −k̂)`` into ``B(−v⃗)`` acting on the label ``(1, +k̂)``, while
leaving the parity-invariant rotation factor ``K`` unchanged.  Thus both cases reduce to
building `Λ` with `Boost(ℐ*v⃗)` and using the same extraction here.

This computation is purely algebraic — no angles, branches, or transcendental functions —
and globally nonsingular, since ``ℂℜ(Λu₊)`` has norm ``½e^{φₐ/2} > 0`` everywhere.  An
independent derivation via explicit angles on the sphere, used as the test oracle, is
described on the [Aberration of Gravitational Waves](@ref "Aberration of Gravitational
Waves") documentation page.
"""
function aberration(R′ₚᵢ::Rotor, Λ)
    ΛR′ₚᵢ2u₊ = Λ * R′ₚᵢ * (1 + im * 𝐤)  # Ignore the 1/2 in u₊ since it's normalized away
    return Rotor(ℂreal(ΛR′ₚᵢ2u₊))
end

@testitem "aberration: identity at β=0" tags = [:unit, :fast] setup = [AberrationSetup] begin
    using Scri: aberration
    import Quaternionic: Rotor, QuatVec, Boost, components
    using .AberrationSetup: FloatTypes

    # With no boost, Λ = Boost(ℐ*0) = 1, so Λ*R′ₚᵢ = R′ₚᵢ is purely real, ℂℜ{Λ R′ₚᵢ u₊} =
    # R′ₚᵢ/2, and the K factor is R′ₚᵢ itself, for either sign of ℐ.
    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        let π = T(π)
            v⃗_zero = QuatVec(zero(T), 0, 0)
            Λ = Boost(ℐ * v⃗_zero)
            for R′ₚᵢ ∈ [
                Rotor(one(T), 0, 0, 0),
                Rotor(cos(π/4), sin(π/4), 0, 0),
                Rotor(cos(π/4), 0, sin(π/4), 0),
                Rotor(cos(π/4), 0, 0, sin(π/4)),
                Rotor(cos(π/4), sin(π/4), 0, 0) * Rotor(cos(π/3), 0, sin(π/3), 0),
            ]
                @test components(aberration(R′ₚᵢ, Λ)) ≈ components(R′ₚᵢ) atol = 4eps(T)
            end
        end
    end
end

@testitem "aberration: pole invariance — no tangent rotation along boost axis" tags = [
    :unit, :fast
] setup = [AberrationSetup] begin
    using Scri: aberration
    import Quaternionic: Rotor, QuatVec, Boost, components
    using .AberrationSetup: FloatTypes

    # North pole: identity rotor sends 𝐤 → 𝐤 = v̂ for a z-axis boost.
    # South pole: rotation by π about x sends 𝐤 → −𝐤 (anti-aligned).
    # In both cases the boost is along the pixel's null ray, so it lies entirely in the
    # AN factor of the KAN decomposition, and the K factor is the input rotor unchanged.
    # This holds for both ℐ = ±1 since flipping v⃗ preserves (anti-)alignment.  With no
    # frame rotation, Λ = Boost(ℐ*v⃗).
    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        ϵ = eps(T)
        R_north = Rotor(one(T), 0, 0, 0)
        R_south = Rotor(0, one(T), 0, 0)
        for β ∈ T.([0.1, 0.5, 0.9])
            v⃗ = QuatVec(0, 0, β)
            Λ = Boost(ℐ * v⃗)
            @test components(aberration(R_north, Λ)) ≈ components(R_north) atol=4ϵ
            @test components(aberration(R_south, Λ)) ≈ components(R_south) atol=4ϵ
        end
    end
end
