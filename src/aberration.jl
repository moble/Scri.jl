"""
    aberration(RRₚᵢ, v⃗, ℐ=+1)

Transform from the rotor `RRₚᵢ` in the boosted frame to the corresponding rotor in the rest
frame, given the boost velocity vector `v⃗`, with `ℐ = ±1` selecting future or past null
infinity — see [Future and past null infinity](@ref scri_pm_conventions).

The result is the ``K`` factor of the Iwasawa ``KAN`` decomposition of the full Lorentz
transformation ``Λ =`` `Boost(ℐ*v⃗) * RRₚᵢ`: it points along the same null direction for the
rest observer as ``Λ`` produces for the boosted observer, and carries the tangent-frame
orientation (spin phase) along with it — exactly the information a spin-weighted field
needs.  The ``A`` factor is accounted for by the conformal factor ``κ``, and the ``N``
factor by component mixing, so only ``K`` is computed here, which is cheaper than the full
[`Quaternionic.KAN`](@extref Quaternionic :jl:function:`Quaternionic.KAN`): projecting with
the idempotent ``u₊ = (1+𝐭𝐳)/2`` gives ``ℂℜ(Λu₊) = ½e^{φₐ/2}𝐑_K``, which normalizes to
``𝐑_K``.  See the "Reversing the order" section of [The Lorentz Group](@ref) page for why
``K`` is the right object, and [the `Quaternionic` documentation](@extref Quaternionic
:std:label:`iwasawa-kan`) for the extraction.

The sign `ℐ` enters only through the boost velocity.  The ``u₊`` extraction is tied to the
null vector ``ℓ ∝ (1, +n̂)`` after the pixel rotor carries ``𝐳 → n̂``, which matches the ℐ⁺
labeling directly.  At ℐ⁻ the labeling is antipodal — pixel ``n̂`` labels the null vector
``(1, −n̂)`` — and conjugating by spatial parity turns ``B(v⃗)`` acting on ``(1, −n̂)`` into
``B(−v⃗)`` acting on the label ``(1, +n̂)``, while leaving the parity-invariant rotation
factor ``K`` unchanged.  Thus both cases reduce to `Boost(ℐ*v⃗)` with the same extraction.

This computation is purely algebraic — no angles, branches, or transcendental functions
beyond those in `Boost` — and globally nonsingular, since ``ℂℜ(Λu₊)`` has norm
``½e^{φₐ/2} > 0`` everywhere.  An independent derivation via explicit angles on the
sphere, used as the test oracle, is described in the "Aberration of Gravitational Waves"
page of the documentation.
"""
function aberration(RRₚᵢ::Rotor{T}, v⃗, ℐ::Int=+1) where {T}
    u₊ = (1 + im * 𝐤) / T(2)
    B = Boost(ℐ * v⃗)
    Λ = B * RRₚᵢ
    λ = Λ * u₊
    ρ = ℂreal(λ)
    return Rotor(ρ)
end

@testitem "aberration: identity at β=0" tags = [:unit, :fast] setup = [AberrationSetup] begin
    using Scri: aberration
    import Quaternionic: Rotor, QuatVec, components
    using .AberrationSetup: FloatTypes

    # With no boost, Λ = RRₚᵢ is purely real, so ℂℜ{Λu₊} = RRₚᵢ/2 and the K factor is
    # RRₚᵢ itself, for either sign of ℐ.
    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        let π = T(π)
            v⃗_zero = QuatVec(zero(T), 0, 0)
            for R ∈ [
                Rotor(one(T), 0, 0, 0),
                Rotor(cos(π/4), sin(π/4), 0, 0),
                Rotor(cos(π/4), 0, sin(π/4), 0),
                Rotor(cos(π/4), 0, 0, sin(π/4)),
                Rotor(cos(π/4), sin(π/4), 0, 0) * Rotor(cos(π/3), 0, sin(π/3), 0),
            ]
                @test components(aberration(R, v⃗_zero, ℐ)) ≈ components(R) atol = 4eps(T)
            end
        end
    end
end

@testitem "aberration: pole invariance — no tangent rotation along boost axis" tags = [
    :unit, :fast
] setup = [AberrationSetup] begin
    using Scri: aberration
    import Quaternionic: Rotor, QuatVec, components
    using .AberrationSetup: FloatTypes

    # North pole: identity rotor sends 𝐤 → 𝐤 = v̂ for a z-axis boost.
    # South pole: rotation by π about x sends 𝐤 → −𝐤 (anti-aligned).
    # In both cases the boost is along the pixel's null ray, so it lies entirely in the
    # AN factor of the KAN decomposition, and the K factor is the input rotor unchanged.
    # This holds for both ℐ = ±1 since flipping v⃗ preserves (anti-)alignment.
    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        ϵ = eps(T)
        R_north = Rotor(one(T), 0, 0, 0)
        R_south = Rotor(0, one(T), 0, 0)
        for β ∈ T.([0.1, 0.5, 0.9])
            v⃗ = QuatVec(0, 0, β)
            @test components(aberration(R_north, v⃗, ℐ)) ≈ components(R_north) atol=4ϵ
            @test components(aberration(R_south, v⃗, ℐ)) ≈ components(R_south) atol=4ϵ
        end
    end
end
