@testmodule AberrationSetup begin
    using DoubleFloats: Double64
    const FloatTypes = (Float32, Float64, Double64, BigFloat)
end

@testmodule AberrationOracle begin
    # The superseded implementation of `aberration`, from Appendix C of Boyle (2015),
    # retained verbatim as an independent oracle for the KAN-based implementation in
    # src/aberration.jl.  It works with explicit angles on the sphere — see the "Aberration
    # of Gravitational Waves" page of the documentation for its derivation — and needs a
    # Taylor branch where the angle formulas are ill-conditioned.  That completely different
    # structure is what makes it a truly independent check.
    #
    # The `emitted` keyword corresponds to the primary implementation's `ℐ` argument:
    # `emitted = true` ⟺ `ℐ = +1` (ℐ⁺), `emitted = false` ⟺ `ℐ = -1` (ℐ⁻).
    using Quaternionic: Rotor, QuatVec, 𝐤, absvec, basetype, value

    function aberration(RRₚᵢ, v⃗; emitted::Bool=true)
        # Equation numbers refer to Appendix C of Boyle (2015).
        β = absvec(v⃗)
        φ = atanh(β)
        ε = emitted ? 1 : -1  # +1 for ℐ⁺ (outgoing), -1 for ℐ⁻ (incoming)
        εβ = ε * β            # signed β; ε = -1 flips the rapidity φ → -φ
        k̂′ = RRₚᵢ(𝐤)  # direction in the boosted frame corresponding to this pixel

        # The quaternion product of two pure vectors p, q satisfies pq = -(p⋅q) + p×q,
        # giving the dot product in the scalar part and the cross product in the vector part.
        k̂′v⃗ = k̂′ * v⃗
        k̂′xv⃗ = QuatVec(k̂′v⃗)  # extracts the vector (cross-product) part; zeros the scalar part

        μ = absvec(k̂′xv⃗)        # |k̂′ × v⃗| = β sinΘ̑; vanishes when k̂′ ∥ ±v⃗
        # atan(β sinΘ̑, β cosΘ̑) = Θ̑; the β factors cancel so this is independent of |v⃗|
        Θ̑ = atan(μ, -k̂′v⃗.w)    # Eq. (C6): k̂′v⃗.w = -(k̂′⋅v⃗) = -β cosΘ̑
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

            # Taylor expansion of sin((Θ̑−Θ)/2) / (β sinΘ̑) about β=0.  For ε = -1 the ratio
            # is negative (leading term -1/2); factoring out ε and using εβ inside restores
            # a uniform +1/2 leading term in the series body.  Multiplying by k̂′xv⃗ (which
            # includes the explicit factor of β sinΘ̑) gives the full vector term
            # sin((Θ̑−Θ)/2) * (k̂′×v⃗)/|k̂′×v⃗| without any division by μ.
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

            # exp[k̂′×v⃗/|k̂′×v⃗| * (Θ̑−Θ)/2] to fifth order in |k̂′xv⃗|
            Q = cosΔΘ╱2 + k̂′xv⃗ * sinΔΘ╱2╱βsinΘ̑
            Rotor{basetype(Q)}(Q)
        else
            Θ = 2atan(exp(-ε * φ) * tan(Θ̑ / 2))  # Eq. (C7): rest-frame polar angle; ε flips sign
            exp((k̂′xv⃗/μ) * (Θ̑-Θ)/2)              # Eq. (C8): sign of (Θ̑-Θ) implicitly includes ε
        end

        return B′ * RRₚᵢ
    end
end

@testitem "aberration: KAN K factor matches Boyle (2015) oracle" tags = [:unit, :validation] setup = [
    AberrationSetup, AberrationOracle
] begin
    import Quaternionic: Rotor, QuatVec, Boost, components, absvec
    using .AberrationSetup: FloatTypes
    using Random: Xoshiro

    # The primary implementation extracts the K factor of the Iwasawa KAN decomposition;
    # the oracle works with explicit angles on the sphere.  The two derivations are
    # entirely independent, so agreement pins down both the direction map and — because
    # we compare full rotors, not just directions — the tangent-frame (spin-phase) part
    # of the transformation, for both signs of ℐ.  Both implementations are continuous
    # in β and agree exactly at β = 0, so there is no double-cover sign ambiguity: any
    # sign flip would be a real discrepancy, not test noise.
    rng = Xoshiro(42)
    randrotor(T) = Rotor(T.(randn(rng, 4))...)
    function randdirection(T)
        v = QuatVec(T.(randn(rng, 3))...)
        return v / absvec(v)
    end

    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        ϵ = eps(T)
        emitted = (ℐ == 1)
        rotors = [randrotor(T) for _ ∈ 1:5]
        for β ∈ T.([1e-8, 1e-3, 0.1, 0.5, 0.9, 0.99])
            v⃗ = β * randdirection(T)
            for R ∈ rotors
                @test components(Scri.aberration(R, Boost(ℐ * v⃗))) ≈
                    components(AberrationOracle.aberration(R, v⃗; emitted)) atol = 50eps(T)
            end
        end
        # Near-pole configurations (k̂′ nearly parallel to ±v⃗).  Here the oracle takes its
        # Taylor branch (β·sinΘ̑ < ∛ϵ), which is a series in β alone — accurate for small β,
        # but wrong when β is large and only sinΘ̑ is small (e.g., at δ = 1e-7, β = 0.99, ℐ⁺,
        # the rotor x-component should be e⁻ᵠ·δ/2 ≈ 3.544e-9; the KAN K factor gets this
        # right, while the Taylor branch returns 9.25e-10).  So evaluate the oracle at
        # BigFloat, where β·sinΘ̑ ≫ ∛eps(BigFloat) keeps it in its well-conditioned exact
        # branch; the T-valued inputs embed exactly.
        for δ ∈ T.([1e-7, 1e-3]), β ∈ T.([1e-3, 0.5, 0.99]), θ₀ ∈ (δ, T(π) - δ)
            R = Rotor(cos(θ₀ / 2), sin(θ₀ / 2), 0, 0)  # pixel at angle θ₀ from ẑ
            v⃗ = QuatVec(zero(T), 0, β)
            oracle = AberrationOracle.aberration(
                Rotor{BigFloat}(R), QuatVec{BigFloat}(v⃗); emitted
            )
            @test components(Scri.aberration(R, Boost(ℐ * v⃗))) ≈ T.(components(oracle)) atol =
                50ϵ
        end
    end
end

@testitem "aberration: is the K factor of a genuine KAN factorization" tags = [
    :unit, :validation
] setup = [AberrationSetup] begin
    import Random
    import Quaternionic:
        Rotor, Quaternion, QuatVec, Boost, KAN, components, normalize, randn
    using .AberrationSetup: FloatTypes

    # The docstring's central claim is that `aberration` returns the ``K`` factor of the
    # Iwasawa ``KAN`` decomposition of `Λ * R′ₚᵢ`.  The other tests here verify
    # *consequences* of that claim (the direction map, the spin phase against the Boyle
    # (2015) oracle, the group identities).  This one verifies the claim itself, and does
    # so without trusting `Quaternionic.KAN`'s internals — which compute the same
    # idempotent projection, so bare agreement with `KAN`'s first return value would be a
    # tautology rather than a test.
    #
    # Instead we take `A` and `N` from `KAN` merely as *candidate* factors and check the
    # defining properties directly:
    #
    #   1. `K * A * N` reconstructs `Λ * R′ₚᵢ`,
    #   2. `K` is a real rotor — a pure rotation, so `K ∈ Spin(3)`,
    #   3. `A = coshφₐ╱2 + sinhφₐ╱2 𝐭𝐳` — a boost along the preferred axis, so `A ∈ A`,
    #   4. `N` fixes the null direction `ℓ` — components `(1, ζ, imζ, 0)`, so `N ∈ N`.
    #
    # Iwasawa's theorem says such a factorization is unique, so 1–4 together pin `K` to
    # be *the* ``K`` factor.  This also fails loudly if `Quaternionic` ever changes its
    # preferred time axis or null direction, which the Scri-local oracle cannot detect.
    #
    # β is capped at 0.999: the reconstruction in (1) involves cancellation between
    # ultra-relativistic factors and is conditioned like γ, not like `aberration` itself
    # (whose behavior out to β = 1 - 64eps is covered by the unit-norm tests below).
    rng = Random.Xoshiro(4646)
    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        ϵ = eps(T)
        for β ∈ T.([1//10^6, 1//1000, 1//2, 9//10, 999//1000]),
            R′ₚᵢ ∈ (one(Rotor{T}), Rotor{T}(1, T(1)/10^7, 0, 0), randn(rng, Rotor{T}))

            R = randn(rng, Rotor{T})  # frame rotation, as folded in by the caller
            v⃗ = β * normalize(randn(rng, QuatVec{T}))
            Λ = Boost(ℐ * v⃗) * R
            K = Scri.aberration(R′ₚᵢ, Λ)
            _, A, N = KAN(Λ * R′ₚᵢ)

            # (1) K A N reconstructs the full transformation.  Compare as `Quaternion`s so
            # that no constructor silently renormalizes the product away.
            ΛR′ₚᵢ = Λ * R′ₚᵢ
            scale = maximum(abs, components(ΛR′ₚᵢ))
            @test maximum(abs, components(Quaternion(K * A * N) - Quaternion(ΛR′ₚᵢ))) <
                100ϵ * scale

            # (2) K is a rotation: real components, unit norm.
            @test all(isreal, components(K))
            @test abs(sum(abs2, components(K)) - 1) < 8ϵ

            # (3) A is a boost along 𝐭𝐳 = im𝐤: scalar part real, no 𝐢/𝐣, 𝐤 part imaginary.
            @test abs(imag(A[1])) < 8ϵ
            @test abs(A[2]) < 8ϵ
            @test abs(A[3]) < 8ϵ
            @test abs(real(A[4])) < 8ϵ

            # (4) N is a null rotation fixing ℓ: components (1, ζ, imζ, 0).
            @test abs(N[1] - 1) < 8ϵ
            @test abs(N[4]) < 8ϵ
            @test abs(N[3] - im * N[2]) < 8ϵ
        end
    end
end

@testitem "aberration: geometric sign — equatorial pixel: cosΘ = ℐβ" tags = [
    :unit, :validation, :fast
] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, Boost, 𝐤
    using .AberrationSetup: FloatTypes

    # R maps 𝐤 → 𝐢: rotation by π/2 about y.  This is an equatorial pixel (Θ̑ = π/2) for
    # a z-axis boost.  Setting cosΘ′ = 0 in the aberration formula gives cosΘ = ℐβ exactly:
    # ℐ = +1 (ℐ⁺): rest-frame direction is in the northern hemisphere (cosΘ = +β).
    # ℐ = -1 (ℐ⁻): rest-frame direction is in the southern hemisphere (cosΘ = -β).
    for T ∈ FloatTypes, ℐ ∈ (-1, +1)
        let π = T(π)
            R = Rotor(cos(π/4), 0, sin(π/4), 0)
            for β ∈ T.([0.1, 0.3, 0.5, 0.7, 0.9])
                v⃗ = QuatVec(0, 0, β)
                R_rest = Scri.aberration(R, Boost(ℐ * v⃗))
                k̂_rest = R_rest(𝐤)
                # For pure vectors p, q: (p*q).w = −(p·q), so k̂_rest · ẑ = −(k̂_rest * 𝐤).w.
                cos_Θ = -(k̂_rest * 𝐤).w
                @test cos_Θ ≈ ℐ * β atol = 4eps(T)
                @test ℐ * cos_Θ > 0  # ℐ=+1: northern; ℐ=-1: southern hemisphere
            end
        end
    end
end

@testitem "aberration: round-trip with inverse gives identity" tags = [:unit, :validation] setup = [
    AberrationSetup
] begin
    import Quaternionic: Rotor, QuatVec, Boost, 𝐤, absvec, components, ×̂
    using .AberrationSetup: FloatTypes

    # boosted_rotor(v⃗, R) is an alternative implementation of the ℐ = -1 direction:
    # it uses acos+exp(φ) explicitly, so it cross-validates the numerics of aberration
    # at a slightly looser tolerance.
    function boosted_rotor(v⃗, R)
        β = absvec(v⃗)
        β < 4eps(typeof(β)) && return R
        v̂ = v⃗ / β
        k̂ = R(𝐤)
        φ = atanh(β)
        cos_Θ = clamp(-(k̂ * v̂).w, -one(β), one(β))
        Θ = acos(cos_Θ)
        Θ′ = 2atan(exp(φ) * tan(Θ / 2))
        return exp(((Θ - Θ′) / 2) * (k̂×̂v̂)) * R
    end

    for T ∈ FloatTypes
        let π = T(π)
            rotors = [
                Rotor(one(T), 0, 0, 0),
                Rotor(cos(π/4), sin(π/4), 0, 0),
                Rotor(cos(π/4), 0, sin(π/4), 0),
                Rotor(cos(π/4), sin(π/4), 0, 0) * Rotor(cos(π/3), 0, sin(π/3), 0),
            ]
            for β ∈ T.([0.1, 0.5, 0.9]), R ∈ rotors
                v⃗ = QuatVec(T(0.6) * β, T(-0.8) * β, 0)
                # Cross-check against alternative implementation (looser tolerance due to acos).
                R_boost = boosted_rotor(v⃗, R)
                @test components(Scri.aberration(R_boost, Boost(v⃗))) ≈ components(R) atol =
                    100eps(T)
                # ℐ⁺ and ℐ⁻ are mutual inverses in both directions: the two differ only by
                # v⃗ → -v⃗, and K(B(-v⃗) K(B(v⃗) R)) = R by uniqueness of the Iwasawa
                # decomposition (B(-v⃗) cancels B(v⃗) and the leftover AN factor is absorbed).
                for ℐ ∈ (-1, +1)
                    R_first = Scri.aberration(R, Boost(ℐ * v⃗))
                    @test components(Scri.aberration(R_first, Boost(-ℐ * v⃗))) ≈
                        components(R) atol = 4eps(T)
                end
            end
        end
    end
end

@testitem "aberration: azimuthal symmetry — z-rotation commutes with z-boost" tags = [
    :unit, :validation, :fast
] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, Boost, components
    using .AberrationSetup: FloatTypes

    # For boost along z, any rotation Rz about z preserves v⃗, so
    # aberration(Rz * R, v⃗) = Rz * aberration(R, v⃗).  This follows from the equivariance
    # of the K factor under rotations about the boost axis.
    for T ∈ FloatTypes
        let π = T(π)
            R_base = Rotor(cos(π/4), 0, sin(π/4), 0)  # equatorial pixel
            for β ∈ T.([0.3, 0.7])
                v⃗ = QuatVec(0, 0, β)
                for α ∈ [π/6, π/3, 2π/3]
                    Rz = Rotor(cos(α/2), 0, 0, sin(α/2))  # rotation by α about z
                    lhs = Scri.aberration(Rz * R_base, Boost(v⃗))
                    rhs = Rz * Scri.aberration(R_base, Boost(v⃗))
                    @test components(lhs) ≈ components(rhs) atol = 4eps(T)
                end
            end
        end
    end
end

@testitem "aberration: small-β precision — Float64 matches BigFloat" tags = [
    :unit, :validation, :fast
] begin
    import Quaternionic: Rotor, QuatVec, Boost, components

    # The old implementation needed a 5th-order Taylor branch when β·sinΘ̑ was small; the
    # KAN K-factor extraction is purely algebraic and globally nonsingular, so no branch
    # is needed.  Verify by comparing Float64 against BigFloat evaluation in exactly the
    # regime the old branch existed for — tiny β, and pixels at or near the poles.  All
    # inputs are exactly representable at Float64, so the BigFloat evaluation is the same
    # mathematical function at higher precision.
    for ℐ ∈ (-1, +1),
        β64 ∈ [1e-18, sqrt(eps(Float64)), cbrt(eps(Float64)), 1e-3],
        θ64 ∈ [0.0, 1e-8, 1e-4, 0.5, π/2]  # pixel angle from the boost (z) axis

        θ = BigFloat(θ64)
        R_big = Rotor(cos(θ / 2), sin(θ / 2), 0, 0)
        v⃗_big = QuatVec(zero(BigFloat), 0, BigFloat(β64))
        result_big = Scri.aberration(R_big, Boost(ℐ * v⃗_big))
        result_f64 = Scri.aberration(
            Rotor{Float64}(R_big), Boost(ℐ * QuatVec{Float64}(v⃗_big))
        )
        @test components(Rotor{Float64}(result_big)) ≈ components(result_f64) atol = 4eps()
    end
end

@testitem "aberration: two-boost composition matches vR decomposition (Wigner rotation)" tags = [
    :unit, :validation
] setup = [AberrationSetup] begin
    import Quaternionic
    import Quaternionic: Boost, Rotor, QuatVec, components
    using .AberrationSetup: FloatTypes

    # Applying aberration twice with velocities v⃗₁ then v⃗₂ is equivalent to applying aberration once
    # with the effective boost and pre-rotating by the Wigner rotation:
    #   aberration(aberration(R, v⃗₁), v⃗₂) ≈ aberration(R_Wigner * R, v⃗_eff)
    # where (v⃗_eff, R_Wigner) = vR(L₂ * L₁).
    # vR is not exported from Quaternionic; access as Quaternionic.vR.
    for T ∈ FloatTypes
        let π = T(π)
            s = T(1) / sqrt(T(2))   # 1/√2, exact in T
            boost_pairs = [
                (T(0.5), T[1, 0, 0], T(0.4), T[0, 1, 0]),   # x then y
                (T(0.3), T[s, s, 0], T(0.6), T[0, s, s]),   # diagonals in xy and yz
                (T(0.8), T[0, 0, 1], T(0.5), T[s, 0, s]),   # z then xz-plane
            ]
            rotors = [
                Rotor(one(T), 0, 0, 0),
                Rotor(cos(π/4), sin(π/4), 0, 0),
                Rotor(cos(π/4), 0, sin(π/4), 0),
                Rotor(cos(π/4), sin(π/4), 0, 0) * Rotor(cos(π/3), 0, sin(π/3), 0),
            ]
            for (η₁, k̂₁, η₂, k̂₂) ∈ boost_pairs
                L₁ = Boost(η₁, k̂₁)
                L₂ = Boost(η₂, k̂₂)
                v⃗_eff, R_Wigner = Quaternionic.vR(L₂ * L₁)
                v⃗₁ = tanh(η₁) * QuatVec(k̂₁[1], k̂₁[2], k̂₁[3])
                v⃗₂ = tanh(η₂) * QuatVec(k̂₂[1], k̂₂[2], k̂₂[3])
                for R ∈ rotors
                    lhs = Scri.aberration(Scri.aberration(R, Boost(v⃗₁)), Boost(v⃗₂))
                    rhs = Scri.aberration(R_Wigner * R, Boost(v⃗_eff))
                    @test components(lhs) ≈ components(rhs) atol = 10eps(T)
                end
            end
        end
    end
end

@testitem "aberration: unit-norm preservation at extreme β" tags = [:unit, :fast] setup = [
    AberrationSetup
] begin
    using Scri: aberration
    import Random
    using Quaternionic: Rotor, QuatVec, Boost, components, normalize, randn
    using .AberrationSetup: FloatTypes

    # The K factor is a unit rotor by construction (ℂℜ(Λu₊) normalized), but the
    # normalization passes through floating-point arithmetic; check it survives extreme
    # boosts and arbitrary pixel directions, including near-pole pixels where the
    # numerically vulnerable small-β·sinΘ̑ structure appears.
    rng = Random.Xoshiro(4242)
    for T ∈ (Float32, Float64), ℐ ∈ (-1, +1)
        for β ∈ T.([1e-8, 1e-3, 0.5, 0.99, 1 - 64eps(T)])
            for R ∈ (
                randn(rng, Rotor{T}),
                one(Rotor{T}),                             # pixel at the pole
                Rotor{T}(1, 1e-7, 0, 0),                   # pixel barely off the pole
            )
                v⃗ = β * normalize(randn(rng, QuatVec{T}))
                K = aberration(R, Boost(ℐ * v⃗))
                @test abs(sum(abs2, components(K)) - 1) ≤ 8eps(T)
            end
        end
    end
end

@testitem "aberration: robustness as β → 1" tags = [:unit, :fast] setup = [AberrationSetup] begin
    using Scri: aberration
    using Quaternionic: Rotor, QuatVec, Boost, components
    using .AberrationSetup: FloatTypes

    # The fragile paths of angle-based implementations (atanh(β) → ∞; exp(-φ)tan(Θ̑/2) → 0
    # near the pole) must remain finite here.  Test extreme speeds at equatorial,
    # near-pole, and mid-latitude pixels for both ℐ signs.
    pixels = [
        Rotor{Float64}(1, 0, 0, 0),                  # pole
        Rotor{Float64}(1, 1e-10, 0, 0),              # essentially at the pole
        Rotor{Float64}(cos(π/8), sin(π/8), 0, 0),    # mid-latitude
        Rotor{Float64}(cos(π/4), 0, sin(π/4), 0),    # equator
    ]
    for β ∈ (0.9, 0.99, 0.999, 0.9999), ℐ ∈ (-1, +1)
        v⃗ = QuatVec(0.0, 0.0, β)
        for R ∈ pixels
            K = aberration(R, Boost(ℐ * v⃗))
            c = components(K)
            @test all(isfinite, c)
            @test abs(sum(abs2, c) - 1) ≤ 16eps()
        end
    end
end

@testitem "aberration: general spatial-rotation covariance" tags = [
    :unit, :validation, :fast
] setup = [AberrationSetup] begin
    using Scri: aberration
    import Random
    using Quaternionic: Rotor, QuatVec, Boost, components, normalize, randn
    using .AberrationSetup: FloatTypes

    # For ANY rotation Rf (not just about the boost axis),
    #     aberration(Rf·R, Rf(v⃗), ℐ) = Rf · aberration(R, v⃗, ℐ).
    # The azimuthal-symmetry test is the special case Rf = Rz; this stronger version
    # catches sign errors in the rotation axis that azimuthal symmetry misses.
    rng = Random.Xoshiro(4343)
    for T ∈ (Float32, Float64, BigFloat), ℐ ∈ (-1, +1)
        for _ ∈ 1:5
            R = randn(rng, Rotor{T})
            Rf = randn(rng, Rotor{T})
            β = T(9//10) * rand(rng, T)
            v⃗ = β * normalize(randn(rng, QuatVec{T}))
            lhs = aberration(Rf * R, Boost(ℐ * Rf(v⃗)))
            rhs = Rf * aberration(R, Boost(ℐ * v⃗))
            err = min(
                maximum(abs, components(lhs - rhs)), maximum(abs, components(lhs + rhs))
            )
            @test err < 60eps(T)
        end
    end
end

@testitem "aberration: general aberration-angle formula" tags = [:unit, :validation, :fast] setup = [
    AberrationSetup
] begin
    using Scri: aberration
    using Quaternionic: Rotor, QuatVec, Boost, components, 𝐤
    using LinearAlgebra: dot

    # The full formula cosΘ = (cosΘ̑ + ℐβ)/(1 + ℐβ cosΘ̑) on a grid of pixel angles Θ̑
    # (measured from the boost axis in the boosted frame) and several speeds — extending
    # the equatorial (Θ̑ = π/2) special case tested elsewhere.  Boost along z, pixels in
    # the xz-plane; the rest-frame direction is K(𝐤).
    for β ∈ (0.1, 0.5, 0.9), ℐ ∈ (-1, +1)
        v⃗ = QuatVec(0.0, 0.0, β)
        for Θ̑ ∈ (π/6, π/4, π/3, π/2, 2π/3, 3π/4)
            R = Rotor{Float64}(cos(Θ̑ / 2), 0, sin(Θ̑ / 2), 0)  # rotate 𝐤 by Θ̑ about y
            K = aberration(R, Boost(ℐ * v⃗))
            k̂_rest = K(𝐤)
            cosΘ = dot([k̂_rest.x, k̂_rest.y, k̂_rest.z], [0, 0, 1.0])
            expected = (cos(Θ̑) + ℐ * β) / (1 + ℐ * β * cos(Θ̑))
            @test abs(cosΘ - expected) < 1e-13
        end
    end
end

@testitem "aberration: collinear boost composition" tags = [:unit, :validation] setup = [
    AberrationSetup
] begin
    using Scri: aberration
    import Random
    using Quaternionic: Rotor, QuatVec, Boost, components, randn
    using .AberrationSetup: FloatTypes

    # Two successive parallel boosts equal one boost at the relativistically composed
    # speed: β = (β₁ + β₂)/(1 + β₁β₂).  (For non-collinear boosts a Wigner rotation
    # appears; that case is tested separately.)
    rng = Random.Xoshiro(4444)
    v̂ = QuatVec(0.0, 0.0, 1.0)
    for (β₁, β₂) ∈ ((0.1, 0.2), (0.5, 0.3), (0.9, 0.09), (0.99, -0.5)), ℐ ∈ (-1, +1)
        β = (β₁ + β₂) / (1 + β₁ * β₂)
        for _ ∈ 1:4
            R = randn(rng, Rotor{Float64})
            lhs = aberration(aberration(R, Boost(ℐ * β₁ * v̂)), Boost(ℐ * β₂ * v̂))
            rhs = aberration(R, Boost(ℐ * β * v̂))
            err = min(
                maximum(abs, components(lhs - rhs)), maximum(abs, components(lhs + rhs))
            )
            @test err < 1e-13
        end
    end
end

@testitem "aberration: ForwardDiff derivatives are finite and continuous" tags = [
    :unit, :fast
] setup = [AberrationSetup] begin
    using Scri: aberration
    import ForwardDiff
    import Random
    using Quaternionic: Rotor, QuatVec, Boost, components, randn

    # `aberration` is purely algebraic with no branches, so it should differentiate
    # cleanly — including at small β, where angle-based implementations switch to Taylor
    # branches.  Check every rotor component's derivative with respect to β against
    # central finite differences.
    rng = Random.Xoshiro(4545)
    R = randn(rng, Rotor{Float64})
    for ℐ ∈ (-1, +1), i ∈ 1:4
        f(β) = components(aberration(R, Boost(ℐ * QuatVec(zero(β), zero(β), β))))[i]
        for β₀ ∈ (0.5, 1e-5, 0.0)
            d = ForwardDiff.derivative(f, β₀)
            @test isfinite(d)
            h = 1e-6
            fd = (f(β₀ + h) - f(β₀ - h)) / 2h
            @test abs(d - fd) < 1e-8 * max(1, abs(d))
        end
    end
end
