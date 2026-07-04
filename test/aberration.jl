@testmodule AberrationSetup begin
    using DoubleFloats: Double64
    const FloatTypes = (Float32, Float64, Double64, BigFloat)
end

@testmodule AberrationOracle begin
    # The superseded implementation of `aberration`, from Appendix C of Boyle (2015),
    # retained verbatim as an independent oracle for the KAN-based implementation in
    # src/aberration.jl.  It works with explicit angles on the sphere — see the
    # "Aberration of Gravitational Waves" page of the documentation for its derivation —
    # and needs a Taylor branch where the angle formulas are ill-conditioned.  That
    # completely different structure is what makes it a genuinely independent check.
    #
    # The `emitted` keyword corresponds to the primary implementation's `εᴵ` argument:
    # `emitted = true` ⟺ `εᴵ = +1` (ℐ⁺), `emitted = false` ⟺ `εᴵ = -1` (ℐ⁻).
    using Quaternionic: Rotor, QuatVec, 𝐤, absvec, basetype, value

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
end

@testitem "aberration: KAN K factor matches Boyle (2015) oracle" tags = [:unit, :validation] setup = [
    AberrationSetup, AberrationOracle
] begin
    import Quaternionic: Rotor, QuatVec, components, absvec
    using .AberrationSetup: FloatTypes
    using Random: Xoshiro

    # The primary implementation extracts the K factor of the Iwasawa KAN decomposition;
    # the oracle works with explicit angles on the sphere.  The two derivations are
    # entirely independent, so agreement pins down both the direction map and — because
    # we compare full rotors, not just directions — the tangent-frame (spin-phase) part
    # of the transformation, for both signs of εᴵ.  Both implementations are continuous
    # in β and agree exactly at β = 0, so there is no double-cover sign ambiguity: any
    # sign flip would be a real discrepancy, not test noise.
    rng = Xoshiro(42)
    randrotor(T) = Rotor(T.(randn(rng, 4))...)
    function randdirection(T)
        v = QuatVec(T.(randn(rng, 3))...)
        return v / absvec(v)
    end

    for T ∈ FloatTypes, εᴵ ∈ (-1, +1)
        ϵ = eps(T)
        emitted = (εᴵ == 1)
        rotors = [randrotor(T) for _ ∈ 1:5]
        for β ∈ T.([1e-8, 1e-3, 0.1, 0.5, 0.9, 0.99])
            v⃗ = β * randdirection(T)
            for R ∈ rotors
                @test components(Scri.aberration(R, v⃗, εᴵ)) ≈
                    components(AberrationOracle.aberration(R, v⃗; emitted)) atol = 50eps(T)
            end
        end
        # Near-pole configurations (n̂′ nearly parallel to ±v⃗).  Here the oracle takes
        # its Taylor branch (β·sinΘ̑ < ∛ϵ), which is a series in β alone — accurate for
        # small β, but genuinely wrong when β is large and only sinΘ̑ is small (e.g. at
        # δ = 1e-7, β = 0.99, ℐ⁺, the rotor x-component should be e⁻ᵠ·δ/2 ≈ 3.544e-9;
        # the KAN K factor gets this right, while the Taylor branch returns 9.25e-10).
        # So evaluate the oracle at BigFloat, where β·sinΘ̑ ≫ ∛eps(BigFloat) keeps it in
        # its well-conditioned exact branch; the T-valued inputs embed exactly.
        for δ ∈ T.([1e-7, 1e-3]), β ∈ T.([1e-3, 0.5, 0.99]), θ₀ ∈ (δ, T(π) - δ)
            R = Rotor(cos(θ₀ / 2), sin(θ₀ / 2), 0, 0)  # pixel at angle θ₀ from ẑ
            v⃗ = QuatVec(zero(T), 0, β)
            oracle = AberrationOracle.aberration(
                Rotor{BigFloat}(R), QuatVec{BigFloat}(v⃗); emitted
            )
            @test components(Scri.aberration(R, v⃗, εᴵ)) ≈ T.(components(oracle)) atol = 50ϵ
        end
    end
end

@testitem "aberration: geometric sign — equatorial pixel: cosΘ = εᴵβ" tags = [
    :unit, :validation, :fast
] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, 𝐤
    using .AberrationSetup: FloatTypes

    # R maps 𝐤 → 𝐢: rotation by π/2 about y.  This is an equatorial pixel (Θ̑ = π/2) for
    # a z-axis boost.  Setting cosΘ′ = 0 in the aberration formula gives cosΘ = εᴵβ exactly:
    # εᴵ = +1 (ℐ⁺): rest-frame direction is in the northern hemisphere (cosΘ = +β).
    # εᴵ = -1 (ℐ⁻): rest-frame direction is in the southern hemisphere (cosΘ = -β).
    for T ∈ FloatTypes, εᴵ ∈ (-1, +1)
        let π = T(π)
            R = Rotor(cos(π/4), 0, sin(π/4), 0)
            for β ∈ T.([0.1, 0.3, 0.5, 0.7, 0.9])
                v⃗ = QuatVec(0, 0, β)
                R_rest = Scri.aberration(R, v⃗, εᴵ)
                n̂_rest = R_rest(𝐤)
                # For pure vectors p, q: (p*q).w = −(p·q), so n̂_rest · ẑ = −(n̂_rest * 𝐤).w.
                cos_Θ = -(n̂_rest * 𝐤).w
                @test cos_Θ ≈ εᴵ * β atol = 4eps(T)
                @test εᴵ * cos_Θ > 0  # εᴵ=+1: northern; εᴵ=-1: southern hemisphere
            end
        end
    end
end

@testitem "aberration: round-trip with inverse gives identity" tags = [:unit, :validation] setup = [
    AberrationSetup
] begin
    import Quaternionic: Rotor, QuatVec, 𝐤, absvec, components, ×̂
    using .AberrationSetup: FloatTypes

    # boosted_rotor(v⃗, R) is an alternative implementation of the εᴵ = -1 direction:
    # it uses acos+exp(φ) explicitly, so it cross-validates the numerics of aberration
    # at a slightly looser tolerance.
    function boosted_rotor(v⃗, R)
        β = absvec(v⃗)
        β < 4eps(typeof(β)) && return R
        v̂ = v⃗ / β
        n̂ = R(𝐤)
        φ = atanh(β)
        cos_Θ = clamp(-(n̂ * v̂).w, -one(β), one(β))
        Θ = acos(cos_Θ)
        Θ′ = 2atan(exp(φ) * tan(Θ / 2))
        return exp(((Θ - Θ′) / 2) * (n̂×̂v̂)) * R
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
                @test components(Scri.aberration(R_boost, v⃗)) ≈ components(R) atol =
                    100eps(T)
                # ℐ⁺ and ℐ⁻ are mutual inverses in both directions: the two differ only by
                # v⃗ → -v⃗, and K(B(-v⃗) K(B(v⃗) R)) = R by uniqueness of the Iwasawa
                # decomposition (B(-v⃗) cancels B(v⃗) and the leftover AN factor is absorbed).
                for εᴵ ∈ (-1, +1)
                    R_first = Scri.aberration(R, v⃗, εᴵ)
                    @test components(Scri.aberration(R_first, v⃗, -εᴵ)) ≈ components(R) atol =
                        4eps(T)
                end
            end
        end
    end
end

@testitem "aberration: azimuthal symmetry — z-rotation commutes with z-boost" tags = [
    :unit, :validation, :fast
] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, components
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
                    lhs = Scri.aberration(Rz * R_base, v⃗)
                    rhs = Rz * Scri.aberration(R_base, v⃗)
                    @test components(lhs) ≈ components(rhs) atol = 4eps(T)
                end
            end
        end
    end
end

@testitem "aberration: small-β precision — Float64 matches BigFloat" tags = [
    :unit, :validation, :fast
] begin
    import Quaternionic: Rotor, QuatVec, components

    # The old implementation needed a 5th-order Taylor branch when β·sinΘ̑ was small; the
    # KAN K-factor extraction is purely algebraic and globally nonsingular, so no branch
    # is needed.  Verify by comparing Float64 against BigFloat evaluation in exactly the
    # regime the old branch existed for — tiny β, and pixels at or near the poles.  All
    # inputs are exactly representable at Float64, so the BigFloat evaluation is the same
    # mathematical function at higher precision.
    for εᴵ ∈ (-1, +1),
        β64 ∈ [1e-18, sqrt(eps(Float64)), cbrt(eps(Float64)), 1e-3],
        θ64 ∈ [0.0, 1e-8, 1e-4, 0.5, π/2]  # pixel angle from the boost (z) axis

        θ = BigFloat(θ64)
        R_big = Rotor(cos(θ / 2), sin(θ / 2), 0, 0)
        v⃗_big = QuatVec(zero(BigFloat), 0, BigFloat(β64))
        result_big = Scri.aberration(R_big, v⃗_big, εᴵ)
        result_f64 = Scri.aberration(Rotor{Float64}(R_big), QuatVec{Float64}(v⃗_big), εᴵ)
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
            for (η₁, n̂₁, η₂, n̂₂) ∈ boost_pairs
                L₁ = Boost(η₁, n̂₁);
                L₂ = Boost(η₂, n̂₂)
                v⃗_eff, R_Wigner = Quaternionic.vR(L₂ * L₁)
                v⃗₁ = tanh(η₁) * QuatVec(n̂₁[1], n̂₁[2], n̂₁[3])
                v⃗₂ = tanh(η₂) * QuatVec(n̂₂[1], n̂₂[2], n̂₂[3])
                for R ∈ rotors
                    lhs = Scri.aberration(Scri.aberration(R, v⃗₁), v⃗₂)
                    rhs = Scri.aberration(R_Wigner * R, v⃗_eff)
                    @test components(lhs) ≈ components(rhs) atol = 10eps(T)
                end
            end
        end
    end
end
