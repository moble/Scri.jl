# Tests for impose_reality (src/utilities.jl)

@testmodule RealitySetup begin
    using DoubleFloats: Double64
    const FloatTypes = (Float32, Float64, Double64, BigFloat)
end

@testitem "impose_reality: known analytic values (ℓ=0…3)" tags = [:unit, :fast] begin
    # N=16 modes ordered by increasing ℓ, then m ∈ {−ℓ,…,+ℓ}.
    # Index formula: ℓ²+ℓ+m+1.  Input: αᵢₙ[k] = (2k−1) + 2k·i, εᵅ=1.
    #
    # output[i₊] = (αᵢₙ[i₊] + (−1)ᵐ conj(αᵢₙ[i₋])) / 2,  output[i₋] = (−1)ᵐ conj(output[i₊])
    #
    # ℓ=1 m=1 i=4 j=2:  ((7+8i) − (3−4i))/2 = 2+6i,   −2+6i
    # ℓ=2 m=1 i=8 j=6:  ((15+16i) − (11−12i))/2 = 2+14i,  −2+14i
    # ℓ=2 m=2 i=9 j=5:  ((17+18i) + (9−10i))/2 = 13+4i,   13−4i
    # ℓ=3 m=1 i=14 j=12: ((27+28i) − (23−24i))/2 = 2+26i,  −2+26i
    # ℓ=3 m=2 i=15 j=11: ((29+30i) + (21−22i))/2 = 25+4i,  25−4i
    # ℓ=3 m=3 i=16 j=10: ((31+32i) − (19−20i))/2 = 6+26i,  −6+26i
    αᵢₙ = [ComplexF64(2k - 1, 2k) for k ∈ 1:16]
    α = Scri.impose_reality(αᵢₙ, 3, 1)

    @test α[1]  ==  1 + 0im   # ℓ=0 m=0

    @test α[3]  ==  5 + 0im   # ℓ=1 m= 0
    @test α[4]  ==  2 + 6im   # ℓ=1 m=+1
    @test α[2]  == -2 + 6im   # ℓ=1 m=−1

    @test α[7]  == 13 + 0im   # ℓ=2 m= 0
    @test α[8]  ==  2 + 14im  # ℓ=2 m=+1
    @test α[6]  == -2 + 14im  # ℓ=2 m=−1
    @test α[9]  == 13 + 4im   # ℓ=2 m=+2
    @test α[5]  == 13 - 4im   # ℓ=2 m=−2

    @test α[13] == 25 + 0im   # ℓ=3 m= 0
    @test α[14] ==  2 + 26im  # ℓ=3 m=+1
    @test α[12] == -2 + 26im  # ℓ=3 m=−1
    @test α[15] == 25 + 4im   # ℓ=3 m=+2
    @test α[11] == 25 - 4im   # ℓ=3 m=−2
    @test α[16] ==  6 + 26im  # ℓ=3 m=+3
    @test α[10] == -6 + 26im  # ℓ=3 m=−3
end

@testitem "impose_reality: output satisfies reality condition" tags = [:unit, :fast] setup = [
    RealitySetup,
] begin
    import Random
    using .RealitySetup: FloatTypes

    rng = Random.Xoshiro(42)
    for L ∈ [1, 2, 3, 5, 8]
        α_big = randn(rng, BigFloat, L^2) .+ im .* randn(rng, BigFloat, L^2)
        for T ∈ FloatTypes
            α = Scri.impose_reality(Complex{T}.(α_big), L - 1, 1)
            for ℓ ∈ 0:(L - 1)
                @test imag(α[ℓ^2 + ℓ + 1]) == 0  # m=0 modes are real
                for m ∈ 1:ℓ
                    i = ℓ^2 + ℓ + m + 1
                    j = ℓ^2 + ℓ - m + 1
                    # α[j] is set exactly to (-1)^m*conj(α[i]) in the last line of the loop
                    @test α[j] == (-1)^m * conj(α[i])
                end
            end
        end
    end
end

@testitem "impose_reality: idempotent" tags = [:unit, :fast] setup = [RealitySetup] begin
    import Random
    using .RealitySetup: FloatTypes

    rng = Random.Xoshiro(99)
    for L ∈ [1, 2, 5, 8]
        α_big = randn(rng, BigFloat, L^2) .+ im .* randn(rng, BigFloat, L^2)
        for T ∈ FloatTypes
            αᵢₙ = Complex{T}.(α_big)
            α_once = Scri.impose_reality(αᵢₙ, L - 1, 1)
            α_twice = Scri.impose_reality(α_once, L - 1, 1)
            # The first call already satisfies the condition,
            # so the second call is an identity: arithmetic is exact.
            @test α_twice == α_once
        end
    end
end

@testitem "impose_reality: valid input is a fixed point" tags = [:unit, :fast] setup = [
    RealitySetup,
] begin
    import Random
    using .RealitySetup: FloatTypes

    rng = Random.Xoshiro(17)
    for L ∈ [1, 2, 4, 6]
        # Build a BigFloat array satisfying the condition exactly, then convert.
        # Negation and conj are exact in any float type, so the condition is
        # preserved through Complex{T}.(α_big) for every T.
        α_big = zeros(Complex{BigFloat}, L^2)
        for ℓ ∈ 0:(L - 1)
            α_big[ℓ^2 + ℓ + 1] = randn(rng, BigFloat)  # m=0: real
            for m ∈ 1:ℓ
                i = ℓ^2 + ℓ + m + 1
                j = ℓ^2 + ℓ - m + 1
                α_big[i] = randn(rng, BigFloat) + im * randn(rng, BigFloat)
                α_big[j] = (-1)^m * conj(α_big[i])
            end
        end
        for T ∈ FloatTypes
            α = Complex{T}.(α_big)
            @test Scri.impose_reality(α, L - 1, 1) == α
        end
    end
end

@testitem "impose_reality: output size, padding, and εᵅ scaling" tags = [:unit, :fast] setup = [
    RealitySetup,
] begin
    import Random
    using .RealitySetup: FloatTypes

    rng = Random.Xoshiro(63)
    for L ∈ [2, 4, 6]
        α_big = randn(rng, BigFloat, L^2) .+ im .* randn(rng, BigFloat, L^2)
        for T ∈ FloatTypes
            αᵢₙ = Complex{T}.(α_big)
            # Output length is (ℓₘₐₓ+1)² regardless of input length.
            ℓₘₐₓ = L + 2
            α_padded = Scri.impose_reality(αᵢₙ, ℓₘₐₓ, 1)
            @test length(α_padded) == (ℓₘₐₓ + 1)^2
            # Modes beyond the input ℓₘₐₓ are zero-padded.
            @test all(iszero, α_padded[(L^2 + 1):end])
            # εᵅ scales the output linearly (verified to be exact for real εᵅ).
            c = T(3) / T(2)
            @test Scri.impose_reality(αᵢₙ, L - 1, c) == c .* Scri.impose_reality(αᵢₙ, L - 1, 1)
        end
    end
end

@testitem "impose_reality: rejects non-square-length input" tags = [:unit, :fast] begin
    for N ∈ [2, 3, 5, 6, 7, 10]
        @test_throws AssertionError Scri.impose_reality(zeros(ComplexF64, N), 10, 1)
    end
end

# ── compute_t′ ────────────────────────────────────────────────────────────────

@testitem "compute_t′: β=0, only supertranslation shifts valid range" tags = [:unit, :fast] begin
    import Quaternionic: Rotor, QuatVec
    # With β=0 the rotors are irrelevant; γ=1 and k⁻¹=1 for every pixel.
    Rₚ = [Rotor(1.0, 0.0, 0.0, 0.0), Rotor(1.0, 0.0, 0.0, 0.0), Rotor(1.0, 0.0, 0.0, 0.0)]
    αₚ = [1.0, -0.5, 2.0]   # min = -0.5, max = 2.0
    v⃗  = QuatVec(0.0, 0.0, 0.0)
    t  = collect(range(-10.0, 10.0; length = 101))
    # εᵅ = -1:
    #   t′ₘᵢₙ = max_p(tₘᵢₙ - αₚ[p]) = tₘᵢₙ - min(αₚ) = -10 - (-0.5) = -9.5
    #   t′ₘₐₓ = min_p(tₘₐₓ - αₚ[p]) = tₘₐₓ - max(αₚ) = 10 - 2 = 8
    (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
    @test length(t′) == length(t)
    @test t′[begin] ≈ -9.5
    @test t′[end]   ≈  8.0
end

@testitem "compute_t′: α=0, boost along z contracts range by Doppler factor" tags = [:unit, :fast] begin
    import Quaternionic: Rotor, QuatVec
    β  = 0.5
    γ  = 1 / √(1 - β^2)
    # Pixel at +z (identity rotor): k⁻¹ = γ(1-β)  — blue-shifted, weak constraint
    # Pixel at -z (Rotor(0,1,0,0)): k⁻¹ = γ(1+β)  — red-shifted, binding constraint
    Rₚ = [Rotor(1.0, 0.0, 0.0, 0.0), Rotor(0.0, 1.0, 0.0, 0.0)]
    αₚ = [0.0, 0.0]
    v⃗  = QuatVec(0.0, 0.0, β)
    t  = collect(range(-10.0, 10.0; length = 101))
    (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
    # For tₘᵢₙ<0<tₘₐₓ both bounds are set by the anti-boost (-z) pixel:
    @test length(t′) == length(t)
    @test t′[begin] ≈ -10.0 / (γ * (1 + β))
    @test t′[end]   ≈  10.0 / (γ * (1 + β))
end

@testitem "compute_t′: inlined vdotn formula matches R(𝐤) direct application" tags = [
    :unit, :validation, :fast,
] begin
    import Quaternionic: RotorF64, QuatVecF64, 𝐤
    import Random

    rng = Random.Xoshiro(42)
    for _ ∈ 1:30
        R  = randn(rng, RotorF64)
        v⃗  = randn(rng, QuatVecF64)
        v_x, v_y, v_z = v⃗[2], v⃗[3], v⃗[4]   # QuatVec: [1]=w=0, [2..4]=x,y,z
        w, x, y, z = R[1], R[2], R[3], R[4]
        # Inlined formula used by compute_t′:
        vdotn = 2v_x * (w * y + x * z) + 2v_y * (y * z - w * x) + v_z * (w^2 + z^2 - x^2 - y^2)
        # Reference: apply rotor to 𝐤, then take dot product component-wise
        n̂ = R(𝐤)
        vdotn_ref = v⃗[2] * n̂[2] + v⃗[3] * n̂[3] + v⃗[4] * n̂[4]
        @test vdotn ≈ vdotn_ref atol = 4eps(Float64)
    end
end

@testitem "compute_t′: output preserves length and is strictly monotone" tags = [:unit, :fast] begin
    import Quaternionic: RotorF64, QuatVec
    import Random

    rng = Random.Xoshiro(7)
    for _ ∈ 1:10
        Nm = 20
        β  = 0.1 + 0.6 * rand(rng)
        Rₚ = randn(rng, RotorF64, Nm)
        αₚ = 0.05 .* randn(rng, Nm)   # small enough that range doesn't collapse
        v⃗  = QuatVec(0.0, 0.0, β)
        # Non-uniform time grid centered near zero
        dt  = abs.(randn(rng, 100)) .+ 0.01
        t   = cumsum(dt) .- sum(dt) / 2
        (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
        @test length(t′) == length(t)
        @test all(diff(t′) .> 0)
    end
end

@testitem "compute_t′: identity (β=0, α=0) gives t′=t" tags = [:unit, :fast] begin
    import Quaternionic: RotorF64, QuatVec
    import Random

    rng = Random.Xoshiro(12)
    for _ ∈ 1:5
        Nm  = 10
        Rₚ  = randn(rng, RotorF64, Nm)
        αₚ  = zeros(Nm)
        v⃗   = QuatVec(0.0, 0.0, 0.0)
        t   = collect(range(-5.0, 5.0; length = 51))
        (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
        # With β=0 and α=0, k⁻¹=1 for every pixel, so t′ₘᵢₙ=tₘᵢₙ,
        # t′ₘₐₓ=tₘₐₓ, scale=1, and t′ = t exactly.
        @test t′ ≈ t atol = 4eps(maximum(abs, t))
    end
end

@testitem "compute_t′: tᵪ is a fixed point of the affine map t↦t′(t)" tags = [
    :unit, :fast, :validation,
] begin
    import Quaternionic: RotorF64, QuatVec
    import Random

    rng = Random.Xoshiro(77)
    for _ ∈ 1:10
        Nm  = 15
        β   = 0.1 + 0.5 * rand(rng)
        Rₚ  = randn(rng, RotorF64, Nm)
        αₚ  = 0.05 .* randn(rng, Nm)
        v⃗   = QuatVec(0.0, 0.0, β)
        t   = collect(range(-10.0, 10.0; length = 101))
        (t′, tᵪ) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
        # Reconstruct the affine map from returned endpoints and verify t′(tᵪ) = tᵪ.
        tₘᵢₙ, tₘₐₓ     = t[begin], t[end]
        t′ₘᵢₙ, t′ₘₐₓ   = t′[begin], t′[end]
        scale             = (t′ₘₐₓ - t′ₘᵢₙ) / (tₘₐₓ - tₘᵢₙ)
        t′_at_tᵪ          = t′ₘᵢₙ + scale * (tᵪ - tₘᵢₙ)
        @test t′_at_tᵪ ≈ tᵪ atol = 4eps(max(abs(tᵪ), 1.0))
    end
end

@testitem "compute_t′: uniformly-spaced t gives uniformly-spaced t′" tags = [:unit, :fast] begin
    import Quaternionic: RotorF64, QuatVec
    import Random

    rng = Random.Xoshiro(55)
    for _ ∈ 1:10
        Nm  = 12
        β   = 0.2 + 0.4 * rand(rng)
        Rₚ  = randn(rng, RotorF64, Nm)
        αₚ  = 0.05 .* randn(rng, Nm)
        v⃗   = QuatVec(0.0, 0.0, β)
        t   = collect(range(-8.0, 8.0; length = 81))
        (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
        Δt′ = diff(t′)
        # The map t↦t′ is affine, so equal steps in t produce equal steps in t′.
        # Tolerance is relative to the full range, not the step: the forced
        # t′[end] = t′ₘₐₓ assignment can shift the last step by up to eps(range).
        @test all(isapprox.(Δt′, Δt′[begin]; atol = 4eps(t′[end] - t′[begin])))
    end
end

@testitem "compute_t′: t grid entirely positive works correctly" tags = [:unit, :fast] begin
    import Quaternionic: Rotor, QuatVec

    β     = 0.3
    γ     = 1 / √(1 - β^2)
    k_inv = γ * (1 - β)
    # Single pixel at +z, no supertranslation, boost along +z.
    Rₚ = [Rotor(1.0, 0.0, 0.0, 0.0)]
    αₚ = [0.0]
    v⃗  = QuatVec(0.0, 0.0, β)
    t  = collect(range(5.0, 20.0; length = 101))
    (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
    @test length(t′) == length(t)
    @test t′[begin] ≈ t[begin] / k_inv atol = 4eps(t′[begin])
    @test t′[end]   ≈ t[end]   / k_inv atol = 4eps(t′[end])
end

@testitem "compute_t′: t grid entirely negative works correctly" tags = [:unit, :fast] begin
    import Quaternionic: Rotor, QuatVec

    β     = 0.3
    γ     = 1 / √(1 - β^2)
    k_inv = γ * (1 - β)
    # Single pixel at +z, no supertranslation, boost along +z.
    Rₚ = [Rotor(1.0, 0.0, 0.0, 0.0)]
    αₚ = [0.0]
    v⃗  = QuatVec(0.0, 0.0, β)
    t  = collect(range(-20.0, -5.0; length = 101))
    (t′, _) = Scri.compute_t′(t, αₚ, Rₚ, v⃗)
    @test length(t′) == length(t)
    @test t′[begin] ≈ t[begin] / k_inv atol = 4eps(t′[begin])
    @test t′[end]   ≈ t[end]   / k_inv atol = 4eps(t′[end])
end

@testitem "compute_t′: collapsed range throws an error" tags = [:unit, :fast] begin
    import Quaternionic: Rotor, QuatVec

    # With β=0, k⁻¹=1 for all pixels, so:
    #   t′ₘᵢₙ = tₘᵢₙ − min(αₚ)
    #   t′ₘₐₓ = tₘₐₓ − max(αₚ)
    # Range collapses when max(αₚ)−min(αₚ) > tₘₐₓ−tₘᵢₙ.
    # Here the supertranslation spread (10) exceeds the time range (2).
    Rₚ = [Rotor(1.0, 0.0, 0.0, 0.0), Rotor(1.0, 0.0, 0.0, 0.0)]
    αₚ = [-5.0, 5.0]
    v⃗  = QuatVec(0.0, 0.0, 0.0)
    t  = collect(range(-1.0, 1.0; length = 51))
    @test_throws ErrorException Scri.compute_t′(t, αₚ, Rₚ, v⃗)
end
