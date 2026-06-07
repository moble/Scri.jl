@testmodule AberrationSetup begin
    using DoubleFloats: Double64
    const FloatTypes = (Float32, Float64, Double64, BigFloat)
end

@testitem "aberration: identity at β=0" tags = [:unit, :fast] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, components
    using .AberrationSetup: FloatTypes

    for T ∈ FloatTypes, ε ∈ (-1, +1)
        let π = T(π), emitted = (ε == 1)
            v⃗_zero = QuatVec(zero(T), 0, 0)
            for R ∈ [
                Rotor(one(T), 0, 0, 0),
                Rotor(cos(π/4), sin(π/4), 0, 0),
                Rotor(cos(π/4), 0, sin(π/4), 0),
                Rotor(cos(π/4), 0, 0, sin(π/4)),
                Rotor(cos(π/4), sin(π/4), 0, 0) * Rotor(cos(π/3), 0, sin(π/3), 0),
            ]
                @test components(Scri.aberration(R, v⃗_zero; emitted)) ≈ components(R) atol =
                    4eps(T)
            end
        end
    end
end

@testitem "aberration: pole invariance — no tangent rotation along boost axis" tags = [
    :unit, :fast
] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, components
    using .AberrationSetup: FloatTypes

    # North pole: identity rotor sends 𝐤 → 𝐤 = v̂ for a z-axis boost.
    # South pole: rotation by π about x sends 𝐤 → −𝐤 (anti-aligned).
    # In both cases n̂′ × v⃗ = 0, so B′ = 1 and aberration returns the input unchanged.
    # This holds for both ε = ±1 since the cross product is independent of the sign convention.
    for T ∈ FloatTypes, ε ∈ (-1, +1)
        let emitted = (ε == 1)
            R_north = Rotor(one(T), 0, 0, 0)
            R_south = Rotor(0, one(T), 0, 0)
            for β ∈ T.([0.1, 0.5, 0.9])
                v⃗ = QuatVec(0, 0, β)
                @test components(Scri.aberration(R_north, v⃗; emitted)) ≈
                    components(R_north) atol = 4eps(T)
                @test components(Scri.aberration(R_south, v⃗; emitted)) ≈
                    components(R_south) atol = 4eps(T)
            end
        end
    end
end

@testitem "aberration: geometric sign — equatorial pixel: cosΘ = εβ" tags = [
    :unit, :validation, :fast
] setup = [AberrationSetup] begin
    import Quaternionic: Rotor, QuatVec, 𝐤
    using .AberrationSetup: FloatTypes

    # R maps 𝐤 → 𝐢: rotation by π/2 about y.  This is an equatorial pixel (Θ̑ = π/2) for
    # a z-axis boost.  Setting cosΘ′ = 0 in the aberration formula gives cosΘ = εβ exactly:
    # ε = +1 (emitted/ℐ⁺): rest-frame direction is in the northern hemisphere (cosΘ = +β).
    # ε = -1 (received/ℐ⁻): rest-frame direction is in the southern hemisphere (cosΘ = -β).
    for T ∈ FloatTypes, ε ∈ (-1, +1)
        let π = T(π), emitted = (ε == 1)
            R = Rotor(cos(π/4), 0, sin(π/4), 0)
            for β ∈ T.([0.1, 0.3, 0.5, 0.7, 0.9])
                v⃗ = QuatVec(0, 0, β)
                R_rest = Scri.aberration(R, v⃗; emitted)
                n̂_rest = R_rest(𝐤)
                # For pure vectors p, q: (p*q).w = −(p·q), so n̂_rest · ẑ = −(n̂_rest * 𝐤).w.
                cos_Θ = -(n̂_rest * 𝐤).w
                @test cos_Θ ≈ ε * β atol = 4eps(T)
                @test ε * cos_Θ > 0  # ε=+1: northern; ε=-1: southern hemisphere
            end
        end
    end
end

@testitem "aberration: round-trip with inverse gives identity" tags = [:unit, :validation] setup = [
    AberrationSetup
] begin
    import Quaternionic: Rotor, QuatVec, 𝐤, absvec, components, ×̂
    using .AberrationSetup: FloatTypes

    # boosted_rotor(v⃗, R) is an alternative implementation of the emitted=false direction:
    # it uses acos+exp(φ) rather than the Horner-form Taylor branch, so it cross-validates
    # the numerics of aberration at a slightly looser tolerance.
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
                # emitted and received are mutual inverses in both directions.
                for ε ∈ (-1, +1)
                    emitted = (ε == 1)
                    R_first = Scri.aberration(R, v⃗; emitted)
                    @test components(Scri.aberration(R_first, v⃗; emitted=(!emitted))) ≈
                        components(R) atol = 4eps(T)
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
    # aberration(Rz * R, v⃗) = Rz * aberration(R, v⃗).  This follows from the equivariance of B′
    # under rotations about the boost axis.
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

@testitem "aberration:\n    Taylor branch matches exact formula near Float64/BigFloat threshold" tags = [
    :unit, :validation, :fast
] begin
    import Quaternionic: Rotor, QuatVec, components

    # For Float64, the Taylor branch is taken when β·sinΘ̑ < ∛eps(Float64) ≈ 6.06e-6.
    # BigFloat's threshold is ∛eps(BigFloat) ≈ 4e-26, so the same value of β·sinΘ̑ ≈ 5e-6
    # sits firmly in BigFloat's exact-formula branch while forcing Float64 into the Taylor
    # branch.  The two must agree to Float64 precision.
    #
    # Use an equatorial pixel (Θ̑ = π/2, sinΘ̑ = 1) so β·sinΘ̑ = β = target, with no extra
    # geometry to think about.  The rotor for a 90° rotation about y maps 𝐤 → 𝐢.

    cbrt_ε₆₄ = cbrt(eps(Float64))  # ≈ 6.06e-6

    for frac ∈ [0.1, 0.25, 0.5, 0.75, 0.95]
        let π = BigFloat(π)
            target = BigFloat(frac) * cbrt_ε₆₄
            R_big = Rotor(cos(π/4), 0, sin(π/4), 0)
            v⃗_big = QuatVec(0, 0, target)

            # BigFloat: β·sinΘ̑ ≫ ∛eps(BigFloat) → exact branch, high precision
            result_big = Scri.aberration(R_big, v⃗_big)

            # Float64: β·sinΘ̑ < ∛eps(Float64) → Taylor branch
            result_f64 = Scri.aberration(Rotor{Float64}(R_big), QuatVec{Float64}(v⃗_big))

            # Round the BigFloat result down to Float64 and compare component-wise.
            @test components(Rotor{Float64}(result_big)) ≈
                components(Rotor{Float64}(result_f64)) atol = 4eps(Float64)
        end
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
