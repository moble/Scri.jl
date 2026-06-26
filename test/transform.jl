# Tests for `transform!` (src/transform.jl), including the ℐ⁻ (εᴵ = -1) conformal-factor
# convention and the BMS-element bridge method.

@testitem "transform!: pure-boost spin-0 conformal factor (ℐ⁺ and ℐ⁻)" tags = [
    :validation, :integration, :fast
] begin
    using Quaternionic: QuatVec, Rotor, 𝐤, vec, absvec
    using SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
    using LinearAlgebra: lu, dot

    # For a pure boost (no rotation, no supertranslation) acting on a spin-0 ψ₂ field with
    # the rest of its peeling tower set to zero, there is no component mixing, so the law is
    # simply ψ₂′ = κ⁻³ ψ₂ with 1/κ = γ(1 - εᴵ v⃗⋅n̂).  We reconstruct that pointwise from the
    # *output* modes and compare against an independent evaluation — pinning the conformal
    # factor (transform.jl) and the past-vs-future direction map (`aberration`'s `emitted`)
    # at both null infinities.  At εᴵ = +1 this is also a regression guard for ℐ⁺.
    ℓ = 8
    N = (ℓ + 1)^2
    Rs = golden_ratio_spiral_rotors(0, ℓ, Float64)   # the uniform B-frame grid transform! uses
    # A smooth, band-limited, real spin-0 field.
    f = [1.0 + 0.3 * vec(R(𝐤))[3] + 0.2 * vec(R(𝐤))[1]^2 for R ∈ Rs]
    α0 = lu(ₛ𝐘(0, ℓ, Float64, Rs)) \ complex(f)

    v⃗ = QuatVec(0.1, -0.15, 0.35)
    β = absvec(v⃗)
    γ = 1 / √(1 - β^2)
    t = [0.0, 1.0, 2.0, 3.0]

    for (εᴵ, comps) ∈ ((+1, (:ψ₂, :ψ₃, :ψ₄)), (-1, (:ψ₂, :ψ₁, :ψ₀)))
        dc = Scri.DataComponents(comps...; εᴵ)
        data = zeros(ComplexF64, N, 4, 3)
        for j ∈ 1:4
            data[:, j, 1] .= α0   # ψ₂ is the first component in both `comps`; the rest are 0
        end
        Scri.transform!(data, copy(t), v⃗, one(Rotor{Float64}), zeros(ComplexF64, 1), dc)

        # Replicate transform!'s rest-frame grid and the expected pixel values.
        emitted = (εᴵ == 1)
        Rₚ = [Scri.aberration(R′ₚ, v⃗; emitted) for R′ₚ ∈ Rs]   # R = 1, so R*R′ₚ = R′ₚ
        ψ₂_rest = ₛ𝐘(0, ℓ, Float64, Rₚ) * α0                    # input ψ₂ at rest directions
        κ⁻¹ = [γ * (1 - εᴵ * dot(vec(v⃗), vec(Rₚ[p](𝐤)))) for p ∈ eachindex(Rₚ)]
        expected = @. κ⁻¹^3 * ψ₂_rest                            # ψ₂′ = κ⁻³ ψ₂ (no mixing)

        # The square s=0 analysis is exactly invertible, so synthesizing the output modes on
        # the B-frame grid recovers the transformed pixel values.
        out_pixels = ₛ𝐘(0, ℓ, Float64, Rs) * data[:, 1, 1]
        @test maximum(abs, out_pixels .- expected) < 1e-12

        # ψ₂′ is constant in time (constant input, time-independent κ⁻³), so every slice
        # agrees — a check that the time interpolation is exact here.
        @test maximum(abs, data[:, 1, 1] .- data[:, 3, 1]) < 1e-12
    end
end

@testitem "transform!: BMS-element bridge" tags = [:validation, :fast] begin
    using Quaternionic: QuatVec, rotor
    using SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
    using LinearAlgebra: lu

    ℓ = 6
    N = (ℓ + 1)^2
    Rs = golden_ratio_spiral_rotors(0, ℓ, Float64)
    α0 = lu(ₛ𝐘(0, ℓ, Float64, Rs)) \ complex([cospi(0.3) + 0.2 * R[1] for R ∈ Rs])
    t = [0.0, 1.0, 2.0, 3.0]
    mkdata() = (d=zeros(ComplexF64, N, 4, 3); foreach(j -> (d[:, j, 1] .= α0), 1:4); d)

    v⃗ = QuatVec(0.05, 0.1, 0.2)
    R = rotor(0.4, 0.0, 1.0, 0.0)

    # ℐ⁺: the bridge must reproduce the main method called with the unpacked parts.
    g⁺ = Scri.BMS{Float64}(; boost_velocity=v⃗, frame_rotation=R)
    dc⁺ = Scri.DataComponents(:ψ₂, :ψ₃, :ψ₄)
    d_bridge = mkdata()
    Scri.transform!(d_bridge, t, g⁺, dc⁺)
    d_core = mkdata()
    Scri.transform!(
        d_core,
        t,
        Scri.boost_velocity(g⁺),
        Scri.frame_rotation(g⁺),
        Scri.supertranslation(g⁺),
        dc⁺,
    )
    @test d_bridge == d_core

    # ℐ⁻: the same element (BMS is null-infinity-agnostic) applied with ℐ⁻ data components,
    # which carry εᴵ = -1 and the reversed peeling tower (ψ₂ requires ψ₁, ψ₀).
    dc⁻ = Scri.DataComponents(:ψ₂, :ψ₁, :ψ₀; εᴵ=-1)
    d_bridge⁻ = mkdata()
    Scri.transform!(d_bridge⁻, t, g⁺, dc⁻)
    d_core⁻ = mkdata()
    Scri.transform!(
        d_core⁻,
        t,
        Scri.boost_velocity(g⁺),
        Scri.frame_rotation(g⁺),
        Scri.supertranslation(g⁺),
        dc⁻,
    )
    @test d_bridge⁻ == d_core⁻

    # The same element transforms differently on the two null infinities (εᴵ comes from dc).
    @test d_bridge⁻ != d_bridge
end

@testitem "transform!: pure supertranslation pins mixing sign and shear ½" tags = [
    :validation, :integration, :fast
] begin
    import Random
    using Quaternionic: QuatVec, Rotor
    using SphericalFunctions: ₛ𝐘, ð, golden_ratio_spiral_rotors
    using LinearAlgebra: lu

    # A pure supertranslation (no boost, no rotation) with time-independent data isolates the
    # null-rotation mixing and the shear shift at κ = 1.  We check transform! against the
    # analytic laws — ψₙ' = Σₖ C(4−n,k) bᵏ ψₙ₊ₖ with b = −ðα/2, and σ' = σ + ½ð²α — evaluated
    # independently on the same grid.  This pins the *signs* (e.g. that b = −ðα/2, not +ðα/2)
    # and the factor of ½ on σ, which the component-level tests (which take b and ð²α as given)
    # cannot.
    ℓ = 8
    N = (ℓ + 1)^2
    Rs = golden_ratio_spiral_rotors(0, ℓ, Float64)
    rng = Random.Xoshiro(20260625)

    # A real supertranslation α restricted to ℓ ≤ 2 (so the mixing products b·ψ stay within
    # ℓₘₐₓ and the spin-weighted analysis round-trips exactly), with its eth-derivatives on the
    # grid — computed exactly as transform! does internally (same ð, same grid).
    αmodes = zeros(ComplexF64, N)
    αmodes[1:9] .= randn(rng, ComplexF64, 9)   # ℓ ≤ 2 only
    α = Scri.impose_reality(αmodes, ℓ, 1)
    ðα_g = ₛ𝐘(1, ℓ, Float64, Rs) * (ð(0, 0, ℓ, Float64) * α)[2:end]
    ð²α_g = ₛ𝐘(2, ℓ, Float64, Rs) * (ð(1, 0, ℓ, Float64) * ð(0, 0, ℓ, Float64) * α)[5:end]
    b_g = -ðα_g ./ 2

    # Time-independent inputs for ψ₂(s=0), ψ₃(s=−1), ψ₄(s=−2), σ(s=+2), band-limited to ℓ ≤ 3.
    dc = Scri.DataComponents(:ψ₂, :ψ₃, :ψ₄, :σ)
    spins = (0, -1, -2, 2)
    Nᵗ = 4
    modes = [randn(rng, ComplexF64, N) for _ ∈ 1:4]
    for k ∈ 1:4
        modes[k][1:(spins[k] ^ 2)] .= 0  # ℓ < |s| (ignored)
        modes[k][17:end] .= 0          # ℓ > 3, so b·ψ stays within ℓₘₐₓ = 8 (no aliasing)
    end
    in_pix = [ₛ𝐘(spins[k], ℓ, Float64, Rs) * modes[k][(spins[k] ^ 2 + 1):end] for k ∈ 1:4]

    data = zeros(ComplexF64, N, Nᵗ, 4)
    for k ∈ 1:4, j ∈ 1:Nᵗ
        data[:, j, k] .= modes[k]
    end
    Scri.transform!(
        data, collect(0.0:(Nᵗ - 1)), QuatVec(0.0, 0.0, 0.0), one(Rotor{Float64}), α, dc
    )
    out_pix = [
        ₛ𝐘(spins[k], ℓ, Float64, Rs) * data[(spins[k] ^ 2 + 1):end, 1, k] for k ∈ 1:4
    ]

    ψ₂ᵢ, ψ₃ᵢ, ψ₄ᵢ, σᵢ = in_pix
    ψ₂ₒ, ψ₃ₒ, ψ₄ₒ, σₒ = out_pix
    tol = 1e-9 * maximum(abs, vcat(in_pix...))
    @test maximum(abs, ψ₄ₒ .- ψ₄ᵢ) < tol                                    # ψ₄ unmixed
    @test maximum(abs, ψ₃ₒ .- (ψ₃ᵢ .+ b_g .* ψ₄ᵢ)) < tol                    # +b, b = −ðα/2
    @test maximum(abs, ψ₂ₒ .- (ψ₂ᵢ .+ 2 .* b_g .* ψ₃ᵢ .+ b_g .^ 2 .* ψ₄ᵢ)) < tol
    @test maximum(abs, σₒ .- (σᵢ .+ ð²α_g ./ 2)) < tol                      # +½ð²α
end

@testitem "transform!: pure supertranslation at ℐ⁻ uses the conjugate parameter" tags = [
    :validation, :integration, :fast
] begin
    import Random
    using Quaternionic: QuatVec, Rotor
    using SphericalFunctions: ₛ𝐘, ð, golden_ratio_spiral_rotors

    # The ℐ⁻ mirror of the previous test.  The generator is l̃, so the peeling tower must use
    # the *conjugate* parameter b̄ = ð̄v'/2κ (spin −1), running ψ₀(s=+2) → ψ₄(s=−2); the shear
    # shift flips sign, σ' = σ − ½ð²α.  Pins the conjugation in the Eᴵ=−1 branch.
    ℓ = 8
    N = (ℓ + 1)^2
    Rs = golden_ratio_spiral_rotors(0, ℓ, Float64)
    rng = Random.Xoshiro(20260626)

    αmodes = zeros(ComplexF64, N)
    αmodes[1:9] .= randn(rng, ComplexF64, 9)   # ℓ ≤ 2 only
    α = Scri.impose_reality(αmodes, ℓ, 1)
    ðα_g = ₛ𝐘(1, ℓ, Float64, Rs) * (ð(0, 0, ℓ, Float64) * α)[2:end]
    ð²α_g = ₛ𝐘(2, ℓ, Float64, Rs) * (ð(1, 0, ℓ, Float64) * ð(0, 0, ℓ, Float64) * α)[5:end]
    b̄_g = conj.(-ðα_g ./ 2)   # b̄ = conj(b), b = −ðα/2

    # ℐ⁻ tower order: ψ₀ unmixed, builds ψ₁, ψ₂ from it.  ψ₀(s=2), ψ₁(s=1), ψ₂(s=0), σ(s=2).
    dc = Scri.DataComponents(:ψ₂, :ψ₁, :ψ₀, :σ; εᴵ=-1)
    spins = (0, 1, 2, 2)
    Nᵗ = 4
    modes = [randn(rng, ComplexF64, N) for _ ∈ 1:4]
    for k ∈ 1:4
        modes[k][1:(spins[k] ^ 2)] .= 0
        modes[k][17:end] .= 0
    end
    in_pix = [ₛ𝐘(spins[k], ℓ, Float64, Rs) * modes[k][(spins[k] ^ 2 + 1):end] for k ∈ 1:4]

    data = zeros(ComplexF64, N, Nᵗ, 4)
    for k ∈ 1:4, j ∈ 1:Nᵗ
        data[:, j, k] .= modes[k]
    end
    Scri.transform!(
        data, collect(0.0:(Nᵗ - 1)), QuatVec(0.0, 0.0, 0.0), one(Rotor{Float64}), α, dc
    )
    out_pix = [
        ₛ𝐘(spins[k], ℓ, Float64, Rs) * data[(spins[k] ^ 2 + 1):end, 1, k] for k ∈ 1:4
    ]

    ψ₂ᵢ, ψ₁ᵢ, ψ₀ᵢ, σᵢ = in_pix
    ψ₂ₒ, ψ₁ₒ, ψ₀ₒ, σₒ = out_pix
    tol = 1e-9 * maximum(abs, vcat(in_pix...))
    @test maximum(abs, ψ₀ₒ .- ψ₀ᵢ) < tol                                    # ψ₀ unmixed
    @test maximum(abs, ψ₁ₒ .- (ψ₁ᵢ .+ b̄_g .* ψ₀ᵢ)) < tol                    # conjugate parameter
    @test maximum(abs, ψ₂ₒ .- (ψ₂ᵢ .+ 2 .* b̄_g .* ψ₁ᵢ .+ b̄_g .^ 2 .* ψ₀ᵢ)) < tol
    @test maximum(abs, σₒ .- (σᵢ .- ð²α_g ./ 2)) < tol                      # −½ð²α at ℐ⁻
    # The un-conjugated parameter would be wrong (spin-weight mismatch); confirm it differs.
    @test maximum(abs, ψ₁ₒ .- (ψ₁ᵢ .+ (-ðα_g ./ 2) .* ψ₀ᵢ)) > tol
end

@testitem "compute_ðt′╱2κ: λ-formula and the boost×supertranslation cross-term sign" tags = [
    :unit, :validation, :fast
] begin
    import Scri: compute_ðt′╱2κ
    using Quaternionic: QuatVec, Rotor, vec
    using SphericalFunctions: golden_ratio_spiral_rotors
    import Random

    # Direct, alias-free unit test of the helper extracted from `transform!`.  This is the
    # coverage that was missing: the cross term `ðt′╱2κ[2,i]·αₚ[i]` in row 1 is nonzero only
    # when BOTH the boost (`v⃗≠0`, so `[2,i]≠0`) AND a supertranslation (`αₚ≠0`) are present —
    # every full-`transform!` test zeroes one of the two, so a sign error there is invisible.
    rng = Random.Xoshiro(42)
    ℓ = 5
    Rₚ = golden_ratio_spiral_rotors(0, ℓ, Float64)   # arbitrary grid; only the rotors matter
    Nₚ = length(Rₚ)
    v⃗ = QuatVec(0.1, -0.2, 0.3)

    for Eᴵ ∈ (1, -1)
        # Independent `[2,i] = (λˣ+iλʸ)/2(λᶻ−εᴵ)`, with `λ = R̃ v⃗ R` formed by Quaternionic
        # rotor conjugation — a different code path than the hand-expanded polynomial inside
        # `compute_ðt′╱2κ`, so this checks that polynomial.
        b₁ = map(Rₚ) do R
            λ = vec(conj(R) * v⃗ * R)
            return (λ[1] + im * λ[2]) / 2(λ[3] - Eᴵ)
        end

        # Case A — pure boost (no supertranslation): row 2 = b₁, row 1 = 0.
        M = compute_ðt′╱2κ(Rₚ, v⃗, zeros(Nₚ), zeros(ComplexF64, Nₚ), Eᴵ)
        @test M[2, :] ≈ b₁
        @test all(iszero, M[1, :])

        # Case B — boost + a *constant* supertranslation α = δt (so ðα = 0).  Physically
        # `ðt′/2κ|₀ = ð(κ·δt)/2κ = δt·(ðκ/2κ) = δt·b₁`, so `[1,i] = −δt·b₁`.  This pins the
        # cross-term SIGN (minus) from a pure Leibniz identity, independent of the helper.
        δt = 0.7
        Mc = compute_ðt′╱2κ(Rₚ, v⃗, fill(δt, Nₚ), zeros(ComplexF64, Nₚ), Eᴵ)
        @test Mc[1, :] ≈ -δt .* b₁
        @test !isapprox(Mc[1, :], +δt .* b₁; rtol=1e-6)   # the wrong (+) sign would fail here

        # Case C — generic αₚ, ðαₚ: full law `[1,i] = −(b₁·αₚ + ðαₚ/2)`, combining the
        # independently-verified b₁ with the −½ ðα term (checks both signs at once).
        αₚ = randn(rng, Nₚ)
        ðαₚ = randn(rng, ComplexF64, Nₚ)
        Mr = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, Eᴵ)
        @test Mr[1, :] ≈ -(b₁ .* αₚ .+ ðαₚ ./ 2)
    end
end

@testitem "compute_ðt′╱2κ: explicit SphericalFunctions cross-check (small boost, high ℓₘₐₓ)" tags = [
    :validation, :integration, :fast
] begin
    import Scri: compute_ðt′╱2κ
    using Quaternionic: QuatVec, Rotor, 𝐤, vec
    using SphericalFunctions: ₛ𝐘, ð, golden_ratio_spiral_rotors
    using LinearAlgebra: dot
    using DoubleFloats: Double64
    import Random

    # Fully independent, end-to-end numerical check of `compute_ðt′╱2κ`.  We build the actual
    # field `t′(n̂) = κ(n̂)·(t − α(n̂))` on the sphere, take its `ð` with SphericalFunctions,
    # divide by `2κ`, and compare to `ðt′╱2κ[1,:] + t·ðt′╱2κ[2,:]`.  Nothing here touches the
    # hand-expanded `λ = R̃v⃗R` polynomial or the derived `−(b·αₚ + ðα/2)` form: `n̂ = R𝐤R̃` comes
    # from Quaternionic rotor action and `ð` from SphericalFunctions, so this independently
    # verifies BOTH the closed form AND that SphericalFunctions' `ð` matches the convention the
    # closed form was derived in.
    #
    # `κ = 1/(γ(1 − εᴵv⃗·n̂))` is not band-limited, so `t′` is not either — but a tiny boost
    # (β ≈ 1e-3) with a large `ℓₘₐₓ` shrinks the out-of-band tail to `~β^ℓₘₐₓ`, far below
    # roundoff, so the spin-0 analysis (square `ₛ𝐘`) and the `ð` are exact and there is no
    # aliasing.  The cross term `ðt′╱2κ[2,:]·αₚ` (the part that hid two sign errors) is `~1e-3`
    # here — many orders above the agreement, so its sign is genuinely pinned.  Running at both
    # `Float64` (≈1e-11) and `Double64` (≈4e-20) shows the residual tracks precision exactly,
    # confirming it is roundoff/conditioning rather than a real discrepancy.
    rng = Random.Xoshiro(2718)
    for (T, rtol) ∈ ((Float64, 1e-8), (Double64, 1e-18))
        v⃗ = QuatVec{T}(7e-4, -5e-4, 9e-4)        # β ≈ 1.2e-3
        γ = 1 / √(1 - sum(abs2, vec(v⃗)))
        for ℓₘₐₓ ∈ (12, 16), Eᴵ ∈ (1, -1)
            N = (ℓₘₐₓ + 1)^2
            Rₚ = golden_ratio_spiral_rotors(0, ℓₘₐₓ, T)

            # A low-ℓ real supertranslation α (ℓ ≤ 2), so κ·α stays well within the band.
            αmodes = zeros(Complex{T}, N)
            αmodes[1:9] .= [Complex{T}(randn(rng), randn(rng)) for _ ∈ 1:9]
            α = Scri.impose_reality(αmodes, ℓₘₐₓ, 1)

            Y0 = ₛ𝐘(0, ℓₘₐₓ, T, Rₚ)              # N×N, square ⟹ exact spin-0 analysis
            Y1 = ₛ𝐘(1, ℓₘₐₓ, T, Rₚ)              # N×(N−1) spin-1 synthesis
            ð0 = ð(0, 0, ℓₘₐₓ, T)                # spin 0 → 1
            αₚ = real.(Y0 * α)                    # α on the grid
            ðαₚ = Y1 * (ð0 * α)[2:end]            # ðα on the grid (drops the ℓ=0 zero)

            # κ on the grid from n̂ = R𝐤R̃ (Quaternionic rotor action; not the λ polynomial).
            n̂ = [vec(R * 𝐤 * conj(R)) for R ∈ Rₚ]
            κ = [1 / (γ * (1 - Eᴵ * dot(vec(v⃗), n))) for n ∈ n̂]

            M = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, Eᴵ)

            for t ∈ (T(0), T(1), T(10))
                t′ = κ .* (t .- αₚ)                # the field whose ð we want
                ðt′ = Y1 * (ð0 * (Y0 \ t′))[2:end] # ð via analyze → ð → synthesize
                oracle = ðt′ ./ (2 .* κ)
                code = M[1, :] .+ t .* M[2, :]
                @test maximum(abs, code .- oracle) < rtol * (1 + maximum(abs, oracle))
            end

            # Sensitivity: flipping the cross-term sign (the old bug) must break the t=0 match.
            code_wrong = (M[1, :] .+ 2 .* M[2, :] .* αₚ)  # = +b·αₚ − ðα/2 instead of −b·αₚ − ðα/2
            oracle0 = (Y1 * (ð0 * (Y0 \ (κ .* (T(0) .- αₚ))))[2:end]) ./ (2 .* κ)
            @test maximum(abs, code_wrong .- oracle0) > 1e-6
        end
    end
end
