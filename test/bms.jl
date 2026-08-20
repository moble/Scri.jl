# Metamorphic tests for the BMS group (src/bms.jl).
#
# These suites test the composition law *without assuming its derivation*: every expected
# value is produced by an independent oracle — the raw 4-vector action from Quaternionic,
# Minkowski geometry of translations, or Wigner-D matrices — never by the composition law
# itself.  The shared helpers (and the same philosophy) are found in the `BMSTestSetup`
# module defined alongside the implementation in src/bms.jl.

@testitem "BMS: action-consistency of composition" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: random_bms, random_direction, act_ref
    using Quaternionic: components

    # THE metamorphic test: acting with g₂∘g₁ must equal acting with g₁ and then g₂,
    # where the action is computed from raw geometry only (`act_ref`); `compose` is the
    # only Scri group code involved.  This checks every sign, factor, and evaluation
    # point in the composition law at once.
    rng = Random.Xoshiro(1414)

    # Rotation-only Lorentz parts: `compose` is exact, so tolerances are tight.
    for _ ∈ 1:4
        g₁ = random_bms(rng, Float64; ℓₘₐₓ=3, βmax=0)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=0)
        h = g₂ * g₁
        for c_α ∈ (+1, -1), _ ∈ 1:4
            t = 2randn(rng)
            k̂ = random_direction(rng, Float64)
            t₁, k̂₁ = act_ref(g₁.Λ, g₁.α, t, k̂; c_α)
            t₂, k̂₂ = act_ref(g₂.Λ, g₂.α, t₁, k̂₁; c_α)
            tₕ, k̂ₕ = act_ref(h.Λ, h.α, t, k̂; c_α)
            @test abs(tₕ - t₂) < 1e-10
            @test maximum(abs, components(k̂ₕ - k̂₂)) < 1e-12
        end
    end

    # With boosts the composed supertranslation is band-limited only approximately;
    # mild boosts and a generous ℓₘₐₓ keep the truncation tail far below tolerance.
    for _ ∈ 1:4
        g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
        h = Scri.compose(g₂, g₁; ℓₘₐₓ=12)
        for c_α ∈ (+1, -1), _ ∈ 1:4
            t = 2randn(rng)
            k̂ = random_direction(rng, Float64)
            t₁, k̂₁ = act_ref(g₁.Λ, g₁.α, t, k̂; c_α)
            t₂, k̂₂ = act_ref(g₂.Λ, g₂.α, t₁, k̂₁; c_α)
            tₕ, k̂ₕ = act_ref(h.Λ, h.α, t, k̂; c_α)
            @test abs(tₕ - t₂) < 1e-7
            @test maximum(abs, components(k̂ₕ - k̂₂)) < 1e-12
        end
    end
end

@testitem "BMS: action-consistency of composition on ℐ⁻" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: random_bms, random_direction, act_ref
    using Quaternionic: components

    # The metamorphic composition test, repeated on ℐ⁻ (ℐ = -1): acting with g₂∘g₁ must
    # equal acting with g₁ then g₂, where the reference action is computed from raw geometry
    # with the antipodal section 𝐤 = (1, -k̂).  The elements include the ℐ⁻ representation
    # as their `ℐ` type parameter (attached to the same raw modes the oracle reads), and
    # plain `compose` must thread the ℐ⁻ conformal factor and direction map through
    # correctly, independently of the derivation.
    rng = Random.Xoshiro(3141)
    at_ℐ⁻(g) = Scri.BMS(g.Λ, g.α; ℐ=-1)

    # Rotation-only Lorentz parts: `compose` is exact, so tolerances are tight.
    for _ ∈ 1:4
        g₁ = random_bms(rng, Float64; ℓₘₐₓ=3, βmax=0)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=0)
        h = Scri.compose(at_ℐ⁻(g₂), at_ℐ⁻(g₁))
        @test Scri.ℐ(h) == -1
        for c_α ∈ (+1, -1), _ ∈ 1:4
            t = 2randn(rng)
            k̂ = random_direction(rng, Float64)
            t₁, k̂₁ = act_ref(g₁.Λ, g₁.α, t, k̂; c_α, ℐ=-1)
            t₂, k̂₂ = act_ref(g₂.Λ, g₂.α, t₁, k̂₁; c_α, ℐ=-1)
            tₕ, k̂ₕ = act_ref(h.Λ, h.α, t, k̂; c_α, ℐ=-1)
            @test abs(tₕ - t₂) < 1e-10
            @test maximum(abs, components(k̂ₕ - k̂₂)) < 1e-12
        end
    end

    # With boosts the composed supertranslation is band-limited only approximately; mild
    # boosts and a generous ℓₘₐₓ keep the truncation tail far below tolerance.
    for _ ∈ 1:4
        g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
        h = Scri.compose(at_ℐ⁻(g₂), at_ℐ⁻(g₁); ℓₘₐₓ=12)
        for c_α ∈ (+1, -1), _ ∈ 1:4
            t = 2randn(rng)
            k̂ = random_direction(rng, Float64)
            t₁, k̂₁ = act_ref(g₁.Λ, g₁.α, t, k̂; c_α, ℐ=-1)
            t₂, k̂₂ = act_ref(g₂.Λ, g₂.α, t₁, k̂₁; c_α, ℐ=-1)
            tₕ, k̂ₕ = act_ref(h.Λ, h.α, t, k̂; c_α, ℐ=-1)
            @test abs(tₕ - t₂) < 1e-7
            @test maximum(abs, components(k̂ₕ - k̂₂)) < 1e-12
        end
    end
end

@testitem "BMS: Poincaré closure — conjugated translations" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_lorentz, tol
    using Quaternionic: QuatVec, absvec, components

    # Conjugating the spacetime translation by (δt, δx⃗) with a Lorentz transformation
    # must give the translation by the Lorentz-transformed 4-vector — computed here with
    # Quaternionic's 4-vector action, never with the composition law.  The composed
    # integrand is *exactly* band-limited to ℓ ≤ 1 (the 1/κ factor cancels the Möbius
    # denominator), so this is sharp: any error in signs, factors, or evaluation points
    # of `compose` shows up far above the eps-scale tolerance.
    rng = Random.Xoshiro(1515)
    for T ∈ FloatTypes
        for _ ∈ 1:3
            δt = randn(rng, T)
            δx⃗ = QuatVec(randn(rng, T), randn(rng, T), randn(rng, T))
            T_d = BMS{T}(; time_translation=δt, space_translation=δx⃗)
            Λ = random_lorentz(rng, T; βmax=1//2)
            gΛ = BMS(Λ, zeros(Complex{T}, 1))

            c = gΛ * T_d * inv(gΛ)
            d′ = Λ(T[δt, vec(δx⃗)...])
            expected = BMS{T}(; time_translation=d′[1], space_translation=d′[2:4])
            scale = max(one(T), abs(δt) + absvec(δx⃗))
            @test c ≈ expected atol = tol(T, 3) * scale

            # Even with ℓ-headroom, no ℓ ≥ 2 modes appear.
            c₄ = Scri.compose(Scri.compose(gΛ, T_d), inv(gΛ); ℓₘₐₓ=4, ℓʷ=9)
            @test maximum(abs, c₄.α[5:end]) < tol(T, 9) * scale
        end
    end
end

@testitem "BMS: Doppler factors from boost conjugation" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: α_value, tol
    using Quaternionic: QuatVec

    # Special case of Poincaré closure with analytic pole values: conjugating a pure time
    # translation by a z-boost yields a translation whose values at the poles are the
    # Doppler factors e^{±η} δt.
    for η ∈ (0.3, 0.6, 1.0)
        β = tanh(η)
        δt = 1.75
        g = BMS{Float64}(; boost_velocity=QuatVec(0.0, 0.0, β))
        T_δt = BMS{Float64}(; time_translation=δt)
        c = g * T_δt * inv(g)

        # Expected: translation by d′ = Λ(δt, 0⃗), with Λ = lorentz(g) acting as a
        # 4-vector map; for this passive convention d′ = γδt(1, -v⃗).
        γ = 1 / √(1 - β^2)
        expected = BMS{Float64}(;
            time_translation=γ * δt, space_translation=QuatVec(0.0, 0.0, -γ * δt * β)
        )
        @test c ≈ expected atol = tol(Float64, 3) * δt * γ * 2

        # Pole values are the Doppler factors.
        @test α_value(c.α, QuatVec(0.0, 0.0, +1.0)) ≈ exp(+η) * δt rtol = 1e-12
        @test α_value(c.α, QuatVec(0.0, 0.0, -1.0)) ≈ exp(-η) * δt rtol = 1e-12

        # The Lorentz part cancels exactly up to rotor-product roundoff.
        @test Scri.is_identity_rotor(c.Λ) ||
            maximum(abs, Quaternionic.components(c.Λ - one(c.Λ))) < 1e-14
    end
end

@testitem "BMS: κ-direction worked example from the docs" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: random_direction, α_value, tol
    using Quaternionic: QuatVec, absvec
    using LinearAlgebra: dot

    # The "Which way around is 𝐾?" note in docs/src/80-details/30-bms_group.md works
    # through the conjugation Λ⁻¹ ∘ δt ∘ Λ for a pure boost with velocity v⃗, claiming
    # the result is the pure supertranslation δt/κ = γ δt (1 - v⃗⋅k̂) — equivalently the
    # spacetime translation by the 4-vector Λ⁻¹δ = γ δt (1; v⃗).  Pin all three faces of
    # that worked example, with κ taken from the package's own `conformal_factor`.
    rng = Random.Xoshiro(2424)
    for _ ∈ 1:5
        v⃗ = (0.7rand(rng)) * random_direction(rng, Float64)
        β = absvec(v⃗)
        γ = 1 / √(1 - β^2)
        δt = randn(rng)
        g = BMS{Float64}(; boost_velocity=v⃗)  # the note's Λ: κ = 1/(γ(1 - v⃗⋅k̂))
        T_δt = BMS{Float64}(; time_translation=δt)
        scale = max(1.0, γ * abs(δt))

        c = inv(g) * T_δt * g  # the note's Λ⁻¹ ∘ δt ∘ Λ

        # (1) 4-vector face: the result is the translation by Λ⁻¹δ, computed with
        # Quaternionic's 4-vector action, which should equal γ δt (1; v⃗).
        d′ = inv(Scri.lorentz(g))([δt, 0.0, 0.0, 0.0])
        @test d′ ≈ γ * δt .* [1.0, vec(v⃗)...] atol = 1e-13 * scale
        expected = BMS{Float64}(; time_translation=d′[1], space_translation=d′[2:4])
        @test c ≈ expected atol = tol(Float64, 3) * scale

        # (2) Pointwise faces: the supertranslation is δt/κ — with κ the package's own
        # conformal factor of the note's Λ — and equals γ δt (1 - v⃗⋅k̂).
        for _ ∈ 1:5
            k̂ = random_direction(rng, Float64)
            αc = α_value(c.α, k̂)
            @test αc ≈ δt / Scri.conformal_factor(g, k̂) atol = 1e-12 * scale
            @test αc ≈ γ * δt * (1 - dot(vec(v⃗), vec(k̂))) atol = 1e-12 * scale
        end
    end
end

@testitem "BMS: rotation conjugation matches Wigner-D" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup:
        random_rotation,
        random_direction,
        random_α,
        rotate_modes,
        rotor_pointing,
        α_value,
        α_eval,
        tol
    using Quaternionic: Lorentz, components

    # Conjugating a supertranslation by a pure rotation rotates the function on the
    # sphere: (Q,0)∘(1,α)∘(Q,0)⁻¹ = (1, α∘Q⁻¹).  The rotated modes are computed with
    # SphericalFunctions' Wigner-D matrices — an oracle completely independent of the
    # composition machinery.  Everything is exactly band-limited here, so tolerances are
    # eps-scale.
    rng = Random.Xoshiro(1717)
    ℓ = 3
    for _ ∈ 1:5
        Q = random_rotation(rng, Float64)
        α = random_α(rng, Float64, ℓ)
        gQ = BMS(Lorentz(Q), zeros(ComplexF64, 1))
        s = BMS{Float64}(; supertranslation=α)

        c = gQ * s * inv(gQ)
        # `compose`'s default output resolution is one ℓ higher than its inputs' (to make
        # room for the ℓ≤1 factor 1/κ₁), so pad the oracle's modes with zeros to match.
        α_D = rotate_modes(α, Q, ℓ)
        α_D = [α_D; zeros(ComplexF64, length(c.α) - length(α_D))]
        @test maximum(abs, c.α - α_D) < tol(Float64, 2ℓ + 1)
        @test Scri.is_identity_rotor(c.Λ) ||
            maximum(abs, components(c.Λ - one(c.Λ))) < 1e-14

        # Pointwise pin of the `rotate_modes` convention itself: f′(k̂) = f(Q⁻¹ k̂).
        for _ ∈ 1:3
            k̂ = random_direction(rng, Float64)
            @test α_value(α_D, k̂) ≈ α_value(α, conj(Q)(k̂)) atol = 1e-11
        end

        # Rotations cannot mix ℓ: the per-ℓ power spectrum is invariant, and the ℓ=0
        # mode is untouched.
        for ℓ′ ∈ 0:ℓ
            block = (ℓ′ ^ 2 + 1):((ℓ′ + 1) ^ 2)
            @test sum(abs2, c.α[block]) ≈ sum(abs2, α[block]) rtol = 1e-12
        end
        @test c.α[1] ≈ α[1] atol = 1e-13
    end
end

@testitem "BMS: group axioms" tags = [:validation, :fast] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: random_bms, tol

    rng = Random.Xoshiro(1818)

    # Rotation-only elements: everything is exactly band-limited, so the axioms hold to
    # roundoff at the default resolutions.
    for _ ∈ 1:4
        g = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=0)
        ϵ = tol(Float64, 5)
        @test g * inv(g) ≈ one(BMS{Float64}) atol = ϵ
        @test inv(g) * g ≈ one(BMS{Float64}) atol = ϵ
        @test inv(inv(g)) ≈ g atol = ϵ
        g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=0)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=1, βmax=0)
        @test (g * g₁) * g₂ ≈ g * (g₁ * g₂) atol = ϵ
    end

    # Boosted elements: `inv` and `compose` truncate unbounded bandwidth, so the axioms hold
    # only as the resolution grows; with mild boosts and ℓ-headroom the truncation tail sits
    # far below these tolerances.
    for _ ∈ 1:3
        g = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=3//20)
        g⁻¹ = inv(g; ℓₘₐₓ=12)
        @test Scri.compose(g, g⁻¹; ℓₘₐₓ=12) ≈ one(BMS{Float64}) atol = 1e-7
        @test Scri.compose(g⁻¹, g; ℓₘₐₓ=12) ≈ one(BMS{Float64}) atol = 1e-7

        g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=3//20)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=3//20)
        lhs = Scri.compose(Scri.compose(g, g₁; ℓₘₐₓ=8), g₂; ℓₘₐₓ=8)
        rhs = Scri.compose(g, Scri.compose(g₁, g₂; ℓₘₐₓ=8); ℓₘₐₓ=8)
        @test lhs ≈ rhs atol = 1e-6
    end
end

@testitem "BMS: normal subgroup and φ homomorphism" tags = [:validation, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: random_lorentz, random_α, random_direction, ray_map, α_value
    using Quaternionic: components

    # 𝒮 is a normal subgroup: conjugating a pure supertranslation by a Lorentz
    # transformation gives a pure supertranslation, namely φ(Λ)(α) = α∘Λ⁻¹ / κ_{Λ⁻¹} —
    # which we evaluate pointwise from raw geometry, not from the composition law.
    rng = Random.Xoshiro(1919)
    for _ ∈ 1:4
        Λ = random_lorentz(rng, Float64; βmax=3//10)
        α = random_α(rng, Float64, 2)
        gΛ = BMS(Λ, zeros(ComplexF64, 1))
        s = BMS{Float64}(; supertranslation=α)
        c = Scri.compose(Scri.compose(gΛ, s), inv(gΛ); ℓₘₐₓ=10)

        @test Scri.is_identity_rotor(c.Λ) ||
            maximum(abs, components(c.Λ - one(c.Λ))) < 1e-14
        for _ ∈ 1:5
            k̂ = random_direction(rng, Float64)
            κ⁻, k̂⁻ = ray_map(inv(Λ), k̂)
            @test α_value(c.α, k̂) ≈ α_value(α, k̂⁻) / κ⁻ atol = 3e-6
        end
    end

    # φ is a homomorphism: conjugating by g₂g₁ equals conjugating by g₁ then by g₂.
    conjugate(gΛ, s; ℓₘₐₓ) = Scri.compose(Scri.compose(gΛ, s), inv(gΛ; ℓₘₐₓ); ℓₘₐₓ)
    for _ ∈ 1:3
        Λ₁ = random_lorentz(rng, Float64; βmax=1//5)
        Λ₂ = random_lorentz(rng, Float64; βmax=1//5)
        α = random_α(rng, Float64, 2)
        s = BMS{Float64}(; supertranslation=α)
        g₁ = BMS(Λ₁, zeros(ComplexF64, 1))
        g₂ = BMS(Λ₂, zeros(ComplexF64, 1))
        g₂₁ = BMS(Λ₂ * Λ₁, zeros(ComplexF64, 1))
        c_nested = conjugate(g₂, conjugate(g₁, s; ℓₘₐₓ=10); ℓₘₐₓ=10)
        c_direct = conjugate(g₂₁, s; ℓₘₐₓ=10)
        @test c_nested ≈ c_direct atol = 1e-5
    end
end

@testitem "BMS: ℒ is not normal — conjugation by supertranslations" tags = [
    :validation, :fast
] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: random_lorentz, random_α, random_direction, ray_map, α_value, tol
    using Quaternionic: QuatVec, Boost, Lorentz
    using LinearAlgebra: norm

    # The complement of the normal-subgroup statement: conjugating a Lorentz
    # transformation by a supertranslation, (1,α)(Λ,0)(1,−α) = (Λ, α∘Λ/κ_Λ − α), leaves
    # the Lorentz part untouched but generically acquires a supertranslation part — so ℒ
    # is *not* a normal subgroup of BMS.
    rng = Random.Xoshiro(2323)
    conjugate(s, gΛ; ℓₘₐₓ) = Scri.compose(Scri.compose(s, gΛ; ℓₘₐₓ), inv(s); ℓₘₐₓ)
    for _ ∈ 1:4
        Λ = random_lorentz(rng, Float64; βmax=3//10)
        α = random_α(rng, Float64, 2)
        s = BMS{Float64}(; supertranslation=α)
        gΛ = BMS(Λ, zeros(ComplexF64, 1))
        c = conjugate(s, gΛ; ℓₘₐₓ=10)

        # Conjugation by 𝒮 cannot change the image in the quotient BMS/𝒮 ≅ ℒ: the
        # Lorentz part is preserved *exactly* (the ±1 factors short-circuit).
        @test c.Λ == Λ

        # The supertranslation part matches α∘Λ/κ_Λ − α, evaluated from raw geometry.
        for _ ∈ 1:5
            k̂ = random_direction(rng, Float64)
            κ, k̂′ = ray_map(Λ, k̂)
            @test α_value(c.α, k̂) ≈ α_value(α, k̂′) / κ - α_value(α, k̂) atol = 3e-6
        end
    end

    # Sharp positive case: constant α = δt conjugating a pure z-boost.  Then
    # α∘Λ/κ_Λ − α = δt(1/κ_Λ − 1) = δt(γ − 1) + δt γ β k̂ᶻ — an exact ℓ ≤ 1 closed form,
    # bounded well away from zero, so the conjugate is definitely not in ℒ.
    η = 0.6
    β = tanh(η)
    γ = 1 / √(1 - β^2)
    δt = 1.25
    s = BMS{Float64}(; time_translation=δt)
    gB = BMS(Boost(η, QuatVec(0.0, 0.0, 1.0)), zeros(ComplexF64, 1))
    c = conjugate(s, gB; ℓₘₐₓ=4)
    @test c.Λ == gB.Λ
    expected = Scri.translation_modes(
        Float64, δt * (γ - 1), QuatVec(0.0, 0.0, -δt * γ * β), 4
    )
    @test maximum(abs, c.α - expected) < tol(Float64, 9) * δt * γ
    @test norm(c.α) > δt * (γ - 1)                  # ⇒ c ∉ ℒ
    @test maximum(abs, c.α[5:end]) < tol(Float64, 9) * δt * γ  # exactly ℓ ≤ 1

    # Boundary case showing the claim is about *generic* elements: if α is invariant
    # under Λ (here: constant α, pure rotation, so κ ≡ 1 and α∘Λ = α), the conjugate
    # falls back into ℒ.
    R = BMSTestSetup.random_rotation(rng, Float64)
    gR = BMS(Lorentz(R), zeros(ComplexF64, 1))
    c = conjugate(BMS{Float64}(; time_translation=δt), gR; ℓₘₐₓ=4)
    @test c.Λ == gR.Λ
    @test maximum(abs, c.α) < tol(Float64, 9) * δt
end

@testitem "BMS: convergence with the working bandwidth ℓʷ" tags = [:validation, :slow] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: random_bms
    using LinearAlgebra: norm

    rng = Random.Xoshiro(2020)
    g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=3//10)
    g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=3//10)

    # Boosted composition converges as the working bandwidth grows.
    ref = Scri.compose(g₂, g₁; ℓₘₐₓ=6, ℓʷ=40).α
    errs = [norm(Scri.compose(g₂, g₁; ℓₘₐₓ=6, ℓʷ).α - ref) for ℓʷ ∈ (7, 13, 25)]
    @test errs[1] > 1e-10            # the test is nontrivial: ℓʷ=7 really does alias
    @test errs[2] < errs[1]
    @test errs[3] < errs[2]
    @test errs[3] < errs[1] / 100

    # Rotation-only composition is exact already at the minimal working bandwidth.
    h₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=0)
    h₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=0)
    @test norm(
        Scri.compose(h₂, h₁; ℓₘₐₓ=2, ℓʷ=2).α - Scri.compose(h₂, h₁; ℓₘₐₓ=2, ℓʷ=20).α
    ) < 1e-12
end

@testitem "BMS: multi-precision sweep" tags = [:validation, :slow] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup:
        FloatTypes,
        tol,
        random_rotation,
        random_α,
        random_direction,
        ray_map,
        α_value,
        rotate_modes
    using Quaternionic: QuatVec, Rotor, Lorentz
    using DoubleFloats: Double64

    # Rotation conjugation is exactly band-limited, so its error is pure roundoff and
    # must scale with eps(T).  This catches any hidden Float64 literal in the pipeline
    # (grids, ₛ𝐘, rotors, …).  NOTE: the Wigner-D reference is computed in BigFloat and
    # cast, because `D_matrices` itself loses half its digits at Double64 (≈ 3e-17, i.e.
    # ~√eps(Double64), as of SphericalFunctions 2.2.9) — whereas golden_ratio grids, ₛ𝐘,
    # and the lu solve are all clean at eps(Double64).
    err = Dict{DataType,Float64}()
    rng₀ = Random.Xoshiro(2121)
    Q64 = random_rotation(rng₀, Float64)
    α64 = random_α(rng₀, Float64, 2)
    α_D_ref = rotate_modes(Complex{BigFloat}.(α64), Rotor{BigFloat}(Q64), 2)
    for T ∈ FloatTypes
        Q = Rotor{T}(Q64)
        α = Complex{T}.(α64)
        gQ = Scri.BMS(Lorentz(Q), zeros(Complex{T}, 1))
        s = Scri.BMS{T}(; supertranslation=α)
        c = gQ * s * inv(gQ)
        # Pad the oracle to `compose`'s (one-ℓ-higher) default output resolution.
        α_D = [α_D_ref; zeros(Complex{BigFloat}, length(c.α) - length(α_D_ref))]
        err[T] = Float64(maximum(abs, ComplexF64.(c.α .- α_D)))
        @test err[T] < tol(T, 5)
    end
    @test err[Double64] < err[Float64] / 1e8
    @test err[BigFloat] < err[Float64] / 1e8

    # Boosted composition vs the pointwise oracle: the error is dominated by the
    # (precision-independent) ℓ-truncation, so it must merely stay below the same small
    # bound for every precision.
    for T ∈ (Float64, Double64, BigFloat)
        rng = Random.Xoshiro(2222)  # same seeds for every T
        v⃗64 = 0.2rand(rng) * random_direction(rng, Float64)
        α₁64 = random_α(rng, Float64, 2)
        α₂64 = random_α(rng, Float64, 2)
        g₁ = Scri.BMS{T}(;
            boost_velocity=QuatVec{T}(v⃗64), supertranslation=Complex{T}.(α₁64)
        )
        g₂ = Scri.BMS{T}(; frame_rotation=one(Rotor{T}), supertranslation=Complex{T}.(α₂64))
        h = Scri.compose(g₂, g₁; ℓₘₐₓ=8, ℓʷ=9)
        for _ ∈ 1:5
            k̂ = random_direction(rng, T)
            κ₁, k̂′ = ray_map(g₁.Λ, k̂)
            expected = α_value(g₁.α, k̂) + α_value(g₂.α, k̂′) / κ₁
            @test abs(α_value(h.α, k̂) - expected) < 1e-5
        end
    end
end
