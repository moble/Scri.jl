# Tests for `transform!` (src/transform.jl), including the ℐ⁻ (ℐ = -1) conformal-factor
# convention and the BMS-element bridge method.

@testitem "transform!: pure-boost spin-0 conformal factor (ℐ⁺ and ℐ⁻)" tags = [
    :validation, :integration, :fast
] begin
    using Quaternionic: QuatVec, Rotor, 𝐤, vec, absvec
    using SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
    using LinearAlgebra: lu, dot

    # For a pure boost (no rotation, no supertranslation) acting on a spin-0 ψ₂ field with
    # the rest of its peeling tower set to zero, there is no component mixing, so the law is
    # simply ψ₂′ = κ⁻³ ψ₂ with 1/κ = γ(1 - ℐ v⃗⋅k̂).  We reconstruct that pointwise from the
    # *output* modes and compare against an independent evaluation — pinning the conformal
    # factor (transform.jl) and the past-vs-future direction map (`aberration`'s `ℐ`)
    # at both null infinities.  At ℐ = +1 this is also a regression guard for ℐ⁺.
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

    for (ℐ, comps) ∈ ((+1, (:ψ₂, :ψ₃, :ψ₄)), (-1, (:ψ₂, :ψ₁, :ψ₀)))
        dc = Scri.DataComponents(comps...; ℐ)
        data = zeros(ComplexF64, N, 4, 3)
        for j ∈ 1:4
            data[:, j, 1] .= α0   # ψ₂ is the first component in both `comps`; the rest are 0
        end
        Scri.transform!(data, copy(t), v⃗, one(Rotor{Float64}), zeros(ComplexF64, 1), dc)

        # Replicate transform!'s rest-frame grid and the expected pixel values.
        Rₚ = [Scri.aberration(R′ₚ, v⃗, ℐ) for R′ₚ ∈ Rs]   # R = 1, so R*R′ₚ = R′ₚ
        ψ₂_rest = ₛ𝐘(0, ℓ, Float64, Rₚ) * α0                    # input ψ₂ at rest directions
        κ⁻¹ = [γ * (1 - ℐ * dot(vec(v⃗), vec(Rₚ[p](𝐤)))) for p ∈ eachindex(Rₚ)]
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

    # ℐ⁻: the element applied with ℐ⁻ data components, which is given by ℐ = -1 and the
    # reversed peeling tower (ψ₂ requires ψ₁, ψ₀).  This element has no supertranslation, so
    # its re-representation to match dc is trivial and the raw parts agree.
    dc⁻ = Scri.DataComponents(:ψ₂, :ψ₁, :ψ₀; ℐ=-1)
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

    # The same element transforms differently on the two null infinities (ℐ comes from dc).
    @test d_bridge⁻ != d_bridge

    # With a supertranslation including odd-ℓ content, the bridge must first re-represent
    # the ℐ⁺ element in dc's ℐ⁻ labeling (the antipodal mode flip via the `BMS` conversion
    # constructor), so it agrees with the core method called on the *converted* modes...
    # The data must be time-dependent here: a supertranslation only shifts retarded time, so
    # it is invisible on constant-in-time data.
    function mktdata()
        d = zeros(ComplexF64, N, 4, 3)
        for j ∈ 1:4
            d[:, j, 1] .= (1 + t[j]) .* α0
        end
        return d
    end
    parts = (;
        boost_velocity=v⃗,
        frame_rotation=R,
        time_translation=0.5,
        space_translation=[0.3, -0.2, 0.1],
    )
    gˢ = Scri.BMS{Float64}(; parts...)
    d_bridgeˢ = mktdata()
    Scri.transform!(d_bridgeˢ, t, gˢ, dc⁻)
    gᴵ = Scri.BMS(gˢ; ℐ=-1)
    d_coreˢ = mktdata()
    Scri.transform!(
        d_coreˢ,
        t,
        Scri.boost_velocity(gᴵ),
        Scri.frame_rotation(gᴵ),
        Scri.supertranslation(gᴵ),
        dc⁻,
    )
    @test d_bridgeˢ == d_coreˢ

    # ...whereas an element already in the ℐ⁻ representation passes its modes through
    # unchanged — and the two differ, because the same raw modes mean different functions in
    # the two labelings.
    g⁻ = Scri.BMS{Float64}(; parts..., ℐ=-1)
    d_bridge⁻ˢ = mktdata()
    Scri.transform!(d_bridge⁻ˢ, t, g⁻, dc⁻)
    d_core⁻ˢ = mktdata()
    Scri.transform!(
        d_core⁻ˢ,
        t,
        Scri.boost_velocity(g⁻),
        Scri.frame_rotation(g⁻),
        Scri.supertranslation(g⁻),
        dc⁻,
    )
    @test d_bridge⁻ˢ == d_core⁻ˢ
    @test d_bridgeˢ != d_bridge⁻ˢ
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
    # shift flips sign, σ' = σ − ½ð²α.  Pins the conjugation in the I=−1 branch.
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

    # ℐ⁻ tower order: ψ₀ unmixed, builds ψ₁, ψ₂ from it.  ψ₀(s=2), ψ₁(s=1), ψ₂(s=0).
    # The ℐ⁻ radiative shear is λ (NP's n-congruence shear), spin weight −2.
    dc = Scri.DataComponents(:ψ₂, :ψ₁, :ψ₀, :λ; ℐ=-1)
    spins = (0, 1, 2, -2)
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

    ψ₂ᵢ, ψ₁ᵢ, ψ₀ᵢ, λᵢ = in_pix
    ψ₂ₒ, ψ₁ₒ, ψ₀ₒ, λₒ = out_pix
    tol = 1e-9 * maximum(abs, vcat(in_pix...))
    @test maximum(abs, ψ₀ₒ .- ψ₀ᵢ) < tol                                    # ψ₀ unmixed
    @test maximum(abs, ψ₁ₒ .- (ψ₁ᵢ .+ b̄_g .* ψ₀ᵢ)) < tol                    # conjugate parameter
    @test maximum(abs, ψ₂ₒ .- (ψ₂ᵢ .+ 2 .* b̄_g .* ψ₁ᵢ .+ b̄_g .^ 2 .* ψ₀ᵢ)) < tol
    @test maximum(abs, λₒ .- (λᵢ .+ conj.(ð²α_g) ./ 2)) < tol               # +½ð̄²α at ℐ⁻
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

    for I ∈ (1, -1)
        # Independent `[2,i] = (λˣ+iλʸ)/2(λᶻ−ℐ)`, with `λ = R̃ v⃗ R` formed by Quaternionic
        # rotor conjugation — a different code path than the hand-expanded polynomial inside
        # `compute_ðt′╱2κ`, so this checks that polynomial.
        b₁ = map(Rₚ) do R
            λ = vec(conj(R) * v⃗ * R)
            return (λ[1] + im * λ[2]) / 2(λ[3] - I)
        end

        # Case A — pure boost (no supertranslation): row 2 = b₁, row 1 = 0.
        M = compute_ðt′╱2κ(Rₚ, v⃗, zeros(Nₚ), zeros(ComplexF64, Nₚ), I)
        @test M[2, :] ≈ b₁
        @test all(iszero, M[1, :])

        # Case B — boost + a *constant* supertranslation α = δt (so ðα = 0).  Physically
        # `ðt′/2κ|₀ = ð(κ·δt)/2κ = δt·(ðκ/2κ) = δt·b₁`, so `[1,i] = −δt·b₁`.  This pins the
        # cross-term SIGN (minus) from a pure Leibniz identity, independent of the helper.
        δt = 0.7
        Mc = compute_ðt′╱2κ(Rₚ, v⃗, fill(δt, Nₚ), zeros(ComplexF64, Nₚ), I)
        @test Mc[1, :] ≈ -δt .* b₁
        @test !isapprox(Mc[1, :], +δt .* b₁; rtol=1e-6)   # the wrong (+) sign would fail here

        # Case C — generic αₚ, ðαₚ: full law `[1,i] = −(b₁·αₚ + ðαₚ/2)`, combining the
        # independently-verified b₁ with the −½ ðα term (checks both signs at once).
        αₚ = randn(rng, Nₚ)
        ðαₚ = randn(rng, ComplexF64, Nₚ)
        Mr = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)
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

    # Fully independent, end-to-end numerical check of `compute_ðt′╱2κ`.  We build the
    # actual field `t′(k̂) = κ(k̂)·(t − α(k̂))` on the sphere, take its `ð` with
    # SphericalFunctions, divide by `2κ`, and compare to `ðt′╱2κ[1,:] + t·ðt′╱2κ[2,:]`.
    # Nothing here touches the hand-expanded `λ = R̃v⃗R` polynomial or the derived `−(b·αₚ +
    # ðα/2)` form: `k̂ = R𝐤R̃` comes from Quaternionic rotor action and `ð` from
    # SphericalFunctions, so this independently verifies BOTH the closed form AND that
    # SphericalFunctions' `ð` matches the convention the closed form was derived in.
    #
    # `κ = 1/(γ(1 − ℐv⃗·k̂))` is not band-limited, so `t′` is not either — but a tiny boost
    # (β ≈ 1e-3) with a large `ℓₘₐₓ` shrinks the out-of-band tail to `~β^ℓₘₐₓ`, far below
    # roundoff, so the spin-0 analysis (square `ₛ𝐘`) and the `ð` are exact and there is no
    # aliasing.  The cross term `ðt′╱2κ[2,:]·αₚ` (the part that hid two sign errors) is
    # `~1e-3` here — many orders above the agreement, so its sign is truly pinned down.
    # Running at both `Float64` (≈1e-11) and `Double64` (≈4e-20) shows the residual tracks
    # precision exactly, confirming it is roundoff/conditioning rather than a real
    # discrepancy.
    rng = Random.Xoshiro(2718)
    for (T, rtol) ∈ ((Float64, 1e-8), (Double64, 1e-18))
        v⃗ = QuatVec{T}(7e-4, -5e-4, 9e-4)        # β ≈ 1.2e-3
        γ = 1 / √(1 - sum(abs2, vec(v⃗)))
        for ℓₘₐₓ ∈ (12, 16), I ∈ (1, -1)
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

            # κ on the grid from k̂ = R𝐤R̃ (Quaternionic rotor action; not the λ polynomial).
            k̂ = [vec(R * 𝐤 * conj(R)) for R ∈ Rₚ]
            κ = [1 / (γ * (1 - I * dot(vec(v⃗), n))) for n ∈ k̂]

            M = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)

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

@testitem "transform!: convention commuting square (represent! ∘ transform! commutes)" tags = [
    :validation, :integration, :fast
] begin
    using Quaternionic: QuatVec, rotor
    using Random

    # THE load-bearing conventions test: for data in a generic convention X, transforming
    # natively in X must equal converting to the native SXS convention, transforming there,
    # and converting back.  Because the conversion factors Fₙ differ *between* components,
    # this pins the convention-rescaled mixing parameter — b → (c_l c_m) b on ℐ⁺ and
    # b̄ → b̄/(c_l c_m) on ℐ⁻ — against the factor ratios Fₙ/Fₙ₊ₖ, jointly across the Weyl
    # and Faraday towers, and exercises the F_σ/F_h shift plumbing and t′ agreement.
    rng = Xoshiro(37)
    ℓ, Nᵗ = 6, 12
    N = (ℓ + 1)^2
    # The radiative-shear slot is σ on ℐ⁺ and λ on ℐ⁻.
    comps⁺ = (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News, :φ₀, :φ₁, :φ₂)
    comps⁻ = (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :λ, :h, :News, :φ₀, :φ₁, :φ₂)
    # A maximally generic convention: every sign flipped, dyad rescaled and rotated in
    # phase, strain rescaled.  (c_ð is set too, but must not affect any transform.)
    X = Scri.Conventions(;
        c_s=-1,
        c_R=-1,
        c_ψ=-1,
        c_σ=-1,
        c_λ=-1,
        c_φ=-1,
        c_l=(-√2),
        c_m=cis(π / 4),
        c_h=2,
        c_ð=1 / √2,
    )
    v⃗ = QuatVec(0.02, -0.03, 0.05)
    R = rotor(QuatVec(0.1, 0.2, -0.3))
    α = 0.03 * randn(rng, ComplexF64, 9)
    α[1] = 0.1
    t = collect(range(-10.0, 10.0; length=Nᵗ))
    mkdata() = begin
        Random.seed!(rng, 4)
        d = 1e-2 * randn(rng, ComplexF64, N, Nᵗ, length(comps⁺))
        for j ∈ 2:Nᵗ
            d[:, j, :] .= d[:, 1, :] .* (1 + 0.01j)  # smooth in time
        end
        d
    end

    for ℐ ∈ (+1, -1)
        comps = ℐ == 1 ? comps⁺ : comps⁻
        dcS = Scri.DataComponents(comps...; ℐ)
        dcX = Scri.DataComponents(comps...; ℐ, conventions=X)

        # Path A: convert SXS → X, then transform natively in X.
        dA = mkdata()
        dA, _ = Scri.represent!(dA, dcS, X)
        dA, t′A = Scri.transform!(dA, copy(t), v⃗, R, copy(α), dcX)

        # Path B: transform natively in SXS, then convert → X.
        dB = mkdata()
        dB, t′B = Scri.transform!(dB, copy(t), v⃗, R, copy(α), dcS)
        dB, _ = Scri.represent!(dB, dcS, X)

        @test maximum(abs, dA .- dB) < 1e-12 * maximum(abs, dB)
        @test t′A == t′B
    end
end

@testitem "transform!: c_α is the time-law sign (and the BMS bridge honors it)" tags = [
    :validation, :fast
] begin
    using Quaternionic: QuatVec, rotor
    using Random

    # The sign formerly passed as c_α is now stored in Conventions.c_α: transforming with
    # c_α = −1 must equal transforming with the default conventions and −α, exactly.
    rng = Xoshiro(21)
    ℓ, Nᵗ = 4, 10
    N = (ℓ + 1)^2
    comps = (:ψ₄, :ψ₃, :σ, :h)
    v⃗ = QuatVec(0.01, 0.02, -0.03)
    R = rotor(QuatVec(0.2, -0.1, 0.15))
    α = 0.05 * randn(rng, ComplexF64, 4)
    α[1] = 0.2
    t = collect(range(-8.0, 8.0; length=Nᵗ))
    mkdata() = begin
        Random.seed!(rng, 5)
        d = 1e-2 * randn(rng, ComplexF64, N, Nᵗ, length(comps))
        for j ∈ 2:Nᵗ
            d[:, j, :] .= d[:, 1, :] .* (1 + 0.02j)
        end
        d
    end

    dc₋ = Scri.DataComponents(comps...; conventions=Scri.Conventions(; c_α=-1))
    dc₊ = Scri.DataComponents(comps...)
    d1 = mkdata()
    d1, t′1 = Scri.transform!(d1, copy(t), v⃗, R, copy(α), dc₋)
    d2 = mkdata()
    d2, t′2 = Scri.transform!(d2, copy(t), v⃗, R, -copy(α), dc₊)
    @test d1 == d2
    @test t′1 == t′2

    # The keyword form threads conventions through to the DataComponents it builds.
    d3 = mkdata()
    d3, _ = Scri.transform!(
        d3,
        copy(t),
        v⃗,
        R,
        copy(α);
        data_components=comps,
        conventions=Scri.Conventions(; c_α=-1),
    )
    @test d3 == d1

    # The BMS bridge re-represents the element with the data conventions' supertranslation
    # sign, so it must agree exactly with the manual unpacking at that sign.
    g = Scri.BMS{Float64}(; boost_velocity=v⃗, frame_rotation=R, supertranslation=copy(α))
    d4 = mkdata()
    Scri.transform!(d4, copy(t), g, dc₋)
    gᴵ = Scri.BMS(g; c_α=-1, ℐ=+1)
    d5 = mkdata()
    Scri.transform!(
        d5,
        copy(t),
        Scri.boost_velocity(gᴵ),
        Scri.frame_rotation(gᴵ),
        Scri.supertranslation(gᴵ),
        dc₋,
    )
    @test d4 == d5
end

@testitem "transform!: plain-vector convenience method dispatches" tags = [:unit, :fast] begin
    using Quaternionic: QuatVec, Rotor

    # Regression test: this method's signature was previously written with abstract
    # element types (`Array{Complex}`, `Vector{Complex}`), which no concrete user data
    # could ever match, so the method was dead code.  It must accept plain `Vector`
    # inputs for `v⃗` and `R` and agree exactly with the `QuatVec`/`Rotor` form.
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 8
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    data = randn(ComplexF64, N, Nᵗ, 2)
    α = randn(ComplexF64, N)
    v⃗ = [0.0, 1e-3, 2e-3]
    R = [1.0, 0.0, 0.0, 0.0]
    dc = Scri.DataComponents(:h, :ψ₄)

    d1 = copy(data)
    d1, t′1 = Scri.transform!(d1, copy(t), v⃗, R, copy(α); data_components=(:h, :ψ₄))
    d2 = copy(data)
    d2, t′2 = Scri.transform!(d2, copy(t), QuatVec(v⃗), Rotor(R), copy(α), dc)
    @test d1 == d2
    @test t′1 == t′2
end

@testitem "transform!: user-input validation throws ArgumentError" tags = [:unit, :fast] begin
    using Quaternionic: QuatVec, Rotor

    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 8
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    data = randn(ComplexF64, N, Nᵗ, 2)
    α = zeros(ComplexF64, N)
    v⃗ = QuatVec(0.0, 0.0, 1e-3)
    R = one(Rotor{Float64})
    dc = Scri.DataComponents(:h, :ψ₄)

    # Explicit `ArgumentError`s (not `@assert`s) so validation survives disabled asserts.
    @test_throws ArgumentError Scri.transform!(
        copy(data), t, QuatVec(0.0, 0.0, 1.5), R, copy(α), dc
    )  # |v⃗| ≥ 1
    @test_throws ArgumentError Scri.transform!(copy(data), t[1:3], v⃗, R, copy(α), dc)  # too few time samples
    @test_throws ArgumentError Scri.transform!(
        copy(data)[1:(N - 1), :, :], t, v⃗, R, copy(α), dc
    )  # non-square mode count
    @test_throws ArgumentError Scri.transform!(
        copy(data), t, v⃗, R, copy(α), Scri.DataComponents(:h)
    )  # component-count mismatch
    @test_throws ArgumentError Scri.transform!(
        copy(data), t, v⃗, R, zeros(ComplexF64, N + 11), dc
    )  # α too long
    # Mixed precision: data must be at least as wide as the other inputs.
    @test_throws ArgumentError Scri.transform!(
        randn(ComplexF32, N, Nᵗ, 2), t, v⃗, R, copy(α), dc
    )
end

@testitem "transform!: pure rotation rotates modes and preserves the time grid" tags = [
    :validation, :fast
] begin
    using Quaternionic: QuatVec, Rotor, randn
    using SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
    import Random

    # For a pure rotation (no boost, no supertranslation) κ ≡ 1 and the mixing parameter
    # vanishes, so the transformation is a rigid rotation of each component on the sphere
    # (including the spin phase) with the time array unchanged.  The output pixel R′ₚ
    # samples the input field at the rotor R·R′ₚ — for β = 0 the aberration K factor *is*
    # R·R′ₚ — so we can check the output against a direct rotation of the input, pointwise,
    # including nonzero spin weights.
    rng = Random.Xoshiro(6262)
    ℓ = 6
    N = (ℓ + 1)^2
    Nᵗ = 5
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    R = randn(rng, Rotor{Float64})
    v⃗ = QuatVec(0.0, 0.0, 0.0)
    dc = Scri.DataComponents(:h, :ψ₄)

    data = zeros(ComplexF64, N, Nᵗ, 2)
    for k ∈ 1:2
        modes = Random.randn(rng, ComplexF64, N)
        modes[1:4] .= 0  # ℓ < |s| modes must be zero for s = ∓2 fields
        for j ∈ 1:Nᵗ
            data[:, j, k] .= modes
        end
    end
    data₀ = copy(data)

    data′, t′ = Scri.transform!(data, copy(t), v⃗, R, zeros(ComplexF64, 1), dc)

    # The time array is unchanged (up to roundoff in the affine regrid).
    @test maximum(abs, t′ .- t) ≤ 4 * eps() * maximum(abs, t)

    # Pointwise check on the golden-ratio grid, per component and spin weight.
    Rs = golden_ratio_spiral_rotors(0, ℓ, Float64)
    Rs_rot = [R * R′ₚ for R′ₚ ∈ Rs]
    for (k, s) ∈ ((1, -2), (2, -2))
        Y = ₛ𝐘(s, ℓ, Float64, Rs)
        Y_rot = ₛ𝐘(s, ℓ, Float64, Rs_rot)
        valid = (s ^ 2 + 1):N
        out = Y * data′[valid, 3, k]
        expected = Y_rot * data₀[valid, 3, k]
        @test maximum(abs, out .- expected) < 1e-11
    end
end

@testitem "transform!: pure time translation is exact at the shifted knots" tags = [
    :validation, :fast
] begin
    using Quaternionic: QuatVec, Rotor
    import Random

    # A pure time translation δt has κ ≡ 1, ðα ≡ 0, and ð²α ≡ 0: no conformal scaling and
    # no mixing.  The output grid is t′ = t − δt, and the data at t′[j] are the input data
    # at t′[j] + δt = t[j] — exactly a knot, where the spline is exact.  So the output mode
    # array equals the input, to roundoff.
    rng = Random.Xoshiro(6363)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 30
    t = collect(LinRange(-30.0, 30.0, Nᵗ))
    δt = 3.75
    α = zeros(ComplexF64, 1)
    α[1] = 2 * √π * δt  # ℓ=0 mode of the constant function δt
    dc = Scri.DataComponents(:ψ₄)

    data = randn(rng, ComplexF64, N, Nᵗ, 1)
    data[1:4, :, :] .= 0  # ℓ < |s| modes must be zero on input for the s = -2 field
    data₀ = copy(data)
    data′, t′ = Scri.transform!(
        data, copy(t), QuatVec(0.0, 0.0, 0.0), one(Rotor{Float64}), α, dc
    )

    @test maximum(abs, (t .- δt) .- t′) ≤ 32 * eps() * (abs(δt) + maximum(abs, t))
    @test maximum(abs, data′ .- data₀) < 1e-11 * maximum(abs, data₀)
end

@testitem "transform!: pure time translation off the knots matches direct interpolation" tags = [
    :validation, :fast
] begin
    using Quaternionic: QuatVec, Rotor
    import Random

    # Shift by δt and request an output grid t′ that does NOT land on the input knots
    # (via the t′ keyword).  Each output sample is then the natural cubic spline of the
    # input mode series evaluated at t′[j] + δt, which we replicate directly from the
    # analytic time dependence of the constructed data — smooth and slowly varying, so
    # the spline error is far below the test tolerance.
    rng = Random.Xoshiro(6464)
    ℓ = 2
    N = (ℓ + 1)^2
    Nᵗ = 201
    T = 500.0
    t = collect(LinRange(-T, T, Nᵗ))
    f(m, t) = exp(-(t / 300)^2) * (1 + m / 10) * cis(t / 150 + m)
    δt = 2.6  # incommensurate with the grid spacing of 5
    α = zeros(ComplexF64, 1)
    α[1] = 2 * √π * δt
    dc = Scri.DataComponents(:ψ₂, :ψ₃, :ψ₄)

    data = zeros(ComplexF64, N, Nᵗ, 3)
    for k ∈ 1:3, (i, m) ∈ enumerate((-ℓ):ℓ), j ∈ 1:Nᵗ
        data[ℓ ^ 2 + i, j, k] = f(m + k, t[j])
    end
    # Stay ≳10 knots away from the ends: the natural-spline boundary condition (d̈ = 0)
    # creates an O(f″δt²) error layer at the endpoints that decays geometrically inward.
    t′ = collect(LinRange(-T + 60, T - 60, Nᵗ))
    data′, _ = Scri.transform!(
        data, copy(t), QuatVec(0.0, 0.0, 0.0), one(Rotor{Float64}), α, dc; t′
    )
    maxerr = maximum(
        abs(data′[ℓ ^ 2 + i, j, k] - f(m + k, t′[j] + δt)) for
        k ∈ 1:3, (i, m) ∈ enumerate((-ℓ):ℓ), j ∈ 1:Nᵗ
    )
    # The measured interior spline floor for this waveform and spacing is ≈ 2e-8.
    @test maxerr < 1e-7
end

@testitem "transform!: round trip through a BMS element and its inverse" tags = [
    :validation, :integration
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # Transform by g, then by inv(g), and compare against the *analytic* original data.
    # Each transformation shrinks the span of times with full-sphere coverage, so the
    # round-trip result is returned on a "squeezed" grid; comparing against the closed-form
    # time dependence of the input evaluates the original exactly at those final times, with
    # no interpolation oracle needed.  The tolerance is set by the spline floor and the
    # band-limit truncation of the mild boost.
    rng = Random.Xoshiro(6565)
    ℓ_data = 2   # band limit of the physical content
    ℓ = 8        # padded band limit, with headroom for the transformation
    N = (ℓ + 1)^2
    Nᵗ = 201
    T = 500.0
    t = collect(LinRange(-T, T, Nᵗ))
    f(m, k, t) = exp(-(t / 300)^2) * (1 + m / 10 + k / 7) * cis(t / 150 + m)
    dc = Scri.DataComponents(:h, :ψ₄)

    data = zeros(ComplexF64, N, Nᵗ, 2)
    for k ∈ 1:2, (i, m) ∈ enumerate((-ℓ_data):ℓ_data), j ∈ 1:Nᵗ
        data[ℓ_data ^ 2 + i, j, k] = f(m, k, t[j])
    end

    v⃗ = QuatVec(1e-4, -2e-4, 3e-4)
    R = randn(rng, Rotor{Float64})
    α = Scri.impose_reality(1e-2 * Random.randn(rng, ComplexF64, 9), 2, 1)
    g = Scri.BMS{Float64}(; boost_velocity=v⃗, frame_rotation=R, supertranslation=α)
    g⁻¹ = inv(g; ℓₘₐₓ=6, ℓʷ=16)

    d1, t1 = Scri.transform!(copy(data), copy(t), g, dc)
    d2, t2 = Scri.transform!(d1, t1, g⁻¹, dc)

    scale = maximum(abs, data)
    maxerr = maximum(
        abs(d2[ℓ_data ^ 2 + i, j, k] - f(m, k, t2[j])) for
        k ∈ 1:2, (i, m) ∈ enumerate((-ℓ_data):ℓ_data), j ∈ 1:Nᵗ
    )
    # The floor here is the band-limit truncation of the *inverse element's*
    # supertranslation sector (its α is exact only up to the ℓₘₐₓ/ℓʷ cutoffs of `inv`),
    # measured at ≈ 5e-6 of the data scale for these parameters.
    @test maxerr / scale < 2e-5
    # The modes above the physical band limit must return to (near) zero.
    @test maximum(abs, d2[((ℓ_data + 1) ^ 2 + 5):end, :, :]) / scale < 2e-5
end

@testitem "transform!: sequential transforms match the composed element" tags = [
    :validation, :integration
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # Applying g₁ then g₂ must agree with applying h = g₂∘g₁ once.  The two pipelines
    # construct different default time grids, so we force the sequential pipeline's final
    # step onto the composed run's grid via the t′ keyword, making the outputs directly
    # comparable.  Tolerances are set by the spline floor and band-limit truncation.
    rng = Random.Xoshiro(6666)
    ℓ_data = 2
    ℓ = 8
    N = (ℓ + 1)^2
    Nᵗ = 201
    T = 500.0
    t = collect(LinRange(-T, T, Nᵗ))
    f(m, k, t) = exp(-(t / 300)^2) * (1 + m / 10 + k / 7) * cis(t / 150 + m)
    dc = Scri.DataComponents(:h, :ψ₄)

    data = zeros(ComplexF64, N, Nᵗ, 2)
    for k ∈ 1:2, (i, m) ∈ enumerate((-ℓ_data):ℓ_data), j ∈ 1:Nᵗ
        data[ℓ_data ^ 2 + i, j, k] = f(m, k, t[j])
    end

    g₁ = Scri.BMS{Float64}(;
        boost_velocity=QuatVec(1e-4, 0.0, -2e-4),
        frame_rotation=randn(rng, Rotor{Float64}),
        supertranslation=Scri.impose_reality(1e-2 * Random.randn(rng, ComplexF64, 9), 2, 1),
    )
    g₂ = Scri.BMS{Float64}(;
        boost_velocity=QuatVec(0.0, 2e-4, 1e-4),
        frame_rotation=randn(rng, Rotor{Float64}),
        supertranslation=Scri.impose_reality(1e-2 * Random.randn(rng, ComplexF64, 9), 2, 1),
    )
    h = Scri.compose(g₂, g₁; ℓₘₐₓ=6, ℓʷ=16)

    # The two pipelines have slightly different valid t′ ranges, so shrink the composed
    # run's default grid a little and force BOTH runs onto that common grid.
    _, th_default = Scri.transform!(copy(data), copy(t), h, dc)
    tc = (th_default[begin] + th_default[end]) / 2
    t_common = tc .+ 0.98 .* (th_default .- tc)

    dh, th = Scri.transform!(copy(data), copy(t), h, dc; t′=copy(t_common))

    d1, t1 = Scri.transform!(copy(data), copy(t), g₁, dc)
    d2, t2 = Scri.transform!(d1, t1, g₂, dc; t′=copy(t_common))
    @test t2 == th == t_common
    # As in the round-trip test, the floor is the band-limit truncation of the composed
    # element's supertranslation sector, plus the doubled spline/regrid floor of the
    # sequential pipeline.
    @test maximum(abs, d2 .- dh) / maximum(abs, data) < 2e-5
end

@testitem "transform!: supplied t′ grid reproduces the default and validates its input" tags = [
    :unit, :fast
] begin
    using Quaternionic: QuatVec, Rotor
    import Random

    rng = Random.Xoshiro(6767)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 50
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    v⃗ = QuatVec(0.0, 0.0, 1e-2)
    R = one(Rotor{Float64})
    α = zeros(ComplexF64, N)
    dc = Scri.DataComponents(:h, :ψ₄)
    data = randn(rng, ComplexF64, N, Nᵗ, 2)

    d1, t′1 = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    d2, t′2 = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc; t′=copy(t′1))
    @test d1 == d2
    @test t′1 == t′2

    # Out-of-range, wrong-length, and non-monotone grids are rejected.
    @test_throws ArgumentError Scri.transform!(
        copy(data), copy(t), v⃗, R, copy(α), dc; t′=collect(LinRange(-100.0, 100.0, Nᵗ))
    )
    @test_throws ArgumentError Scri.transform!(
        copy(data), copy(t), v⃗, R, copy(α), dc; t′=t′1[1:(Nᵗ - 1)]
    )
    @test_throws ArgumentError Scri.transform!(
        copy(data), copy(t), v⃗, R, copy(α), dc; t′=reverse(t′1)
    )
end

@testitem "transform!: keyword form defaults to h-first components with a warning" tags = [
    :unit, :fast
] begin
    using Quaternionic: QuatVec, Rotor
    import Random

    rng = Random.Xoshiro(6868)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 8
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    v⃗ = QuatVec(0.0, 0.0, 1e-3)
    R = one(Rotor{Float64})
    α = zeros(ComplexF64, N)
    data = randn(rng, ComplexF64, N, Nᵗ, 2)

    d_default = copy(data)
    @test_logs (:warn, r"Defaulting to data components \(:h, :ψ₄\)") Scri.transform!(
        d_default, copy(t), v⃗, R, copy(α)
    )
    d_explicit, _ = Scri.transform!(
        copy(data), copy(t), v⃗, R, copy(α), Scri.DataComponents(:h, :ψ₄)
    )
    @test d_default == d_explicit
end

@testitem "transform!: declared input band limit ℓₘₐₓ₀ is exact and validated" tags = [
    :unit, :fast
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # Band-limited data zero-padded to a larger array must transform identically whether
    # or not the input band limit is declared — the declaration only skips provably-zero
    # work in the synthesis stage.  A false declaration must throw.
    rng = Random.Xoshiro(6969)
    ℓ₀ = 2
    ℓ = 6
    N = (ℓ + 1)^2
    Nᵗ = 20
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    v⃗ = QuatVec(1e-3, -2e-3, 0.0)
    R = randn(rng, Rotor{Float64})
    α = Scri.impose_reality(1e-2 * Random.randn(rng, ComplexF64, 4), 1, 1)
    dc = Scri.DataComponents(:h, :ψ₄)

    data = zeros(ComplexF64, N, Nᵗ, 2)
    data[5:((ℓ₀ + 1) ^ 2), :, :] .= Random.randn(rng, ComplexF64, (ℓ₀ + 1)^2 - 4, Nᵗ, 2)

    d1, t1 = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    d2, t2 = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc; ℓₘₐₓ₀=ℓ₀)
    @test d1 == d2
    @test t1 == t2

    # Out-of-range and violated declarations throw.
    @test_throws ArgumentError Scri.transform!(
        copy(data), copy(t), v⃗, R, copy(α), dc; ℓₘₐₓ₀=ℓ + 1
    )
    bad = copy(data)
    bad[(ℓ₀ + 1) ^ 2 + 3, 1, 1] = 1.0
    @test_throws ArgumentError Scri.transform!(bad, copy(t), v⃗, R, copy(α), dc; ℓₘₐₓ₀=ℓ₀)
end

@testitem "transform!: ForwardDiff derivatives with respect to BMS parameters" tags = [
    :unit, :validation
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import ForwardDiff
    import Random

    # `transform!` is differentiable end-to-end with ForwardDiff, provided (a) the output
    # grid is held fixed via the `t′` keyword (otherwise the default grid moves with the
    # parameters, and the derivative mixes in grid motion), and (b) ForwardDiff is loaded,
    # which activates the ScriForwardDiffExt extension: it peels dual types off the
    # parameter-independent pixel grid and analysis factorizations (`primal_float`), whose
    # `qr` would otherwise turn identically-zero perturbations into NaN partials.
    rng = Random.Xoshiro(7070)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 20
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    R₀ = randn(rng, Rotor{Float64})
    dc = Scri.DataComponents(:h, :ψ₄)
    data₀ = Random.randn(rng, ComplexF64, N, Nᵗ, 2)
    data₀[1:4, :, :] .= 0

    # A fixed output grid, safely inside the valid range for all parameter values tested.
    _, t′₀ = Scri.transform!(
        copy(data₀), copy(t), QuatVec(0.0, 0.0, 1e-3), R₀, zeros(ComplexF64, 4), dc
    )
    tc = (t′₀[begin] + t′₀[end]) / 2
    t′fix = tc .+ 0.9 .* (t′₀ .- tc)

    # Boost speed, time translation, and rotation angle, each as a scalar objective.
    function by_boost(β::T) where {T}
        data = Complex{T}.(copy(data₀))
        d, _ = Scri.transform!(
            data,
            copy(t),
            QuatVec(zero(β), zero(β), β),
            R₀,
            zeros(Complex{T}, 4),
            dc;
            t′=t′fix,
        )
        return sum(abs2, d)
    end
    function by_δt(δt::T) where {T}
        data = Complex{T}.(copy(data₀))
        α = zeros(Complex{T}, 4)
        α[1] = 2 * √T(π) * δt
        d, _ = Scri.transform!(data, copy(t), QuatVec(0.0, 0.0, 1e-3), R₀, α, dc; t′=t′fix)
        return sum(abs2, d)
    end
    function by_rotation(θ::T) where {T}
        data = Complex{T}.(copy(data₀))
        Rθ = Rotor{T}(R₀) * Rotor(cos(θ / 2), sin(θ / 2), zero(θ), zero(θ))
        d, _ = Scri.transform!(
            data,
            copy(t),
            QuatVec(zero(θ), zero(θ), zero(θ) + 1e-3),
            Rθ,
            zeros(Complex{T}, 4),
            dc;
            t′=t′fix,
        )
        return sum(abs2, d)
    end

    h = 1e-6
    for (f, x₀) ∈ ((by_boost, 1e-3), (by_δt, 0.5), (by_rotation, 0.1))
        d = ForwardDiff.derivative(f, x₀)
        fd = (f(x₀ + h) - f(x₀ - h)) / 2h
        @test isfinite(d)
        @test abs(d - fd) < 1e-6 * max(1.0, abs(fd))
    end
end

@testmodule PixelReconstruction begin
    import Scri
    import LinearAlgebra
    import SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
    using LinearAlgebra: qr

    # Reconstruct the stage-2 pixel values from `transform!` output, exactly.  The
    # augmented LU of the analysis stage solves [Q⊥ ₛY]·[ξ; modes] = pixels, so the
    # pixels are recovered as ₛY·modes + Q⊥·ξ (plain ₛY·modes for s = 0), with Q⊥
    # rebuilt exactly as in `transform!`.  The synthesis matrices are built at the
    # primal type, so dual-typed mode data reconstruct through primal matrices.
    function reconstruct_pixels(d′::AbstractArray{Complex{T},3}, dc) where {T}
        Nᵐ, Nᵗ, Nᵈ = size(d′)
        ℓₘₐₓ = isqrt(Nᵐ) - 1
        T′ = Scri.primal_float(T)
        R′ₚ = golden_ratio_spiral_rotors(0, ℓₘₐₓ, T′)
        components = typeof(dc).parameters[1]
        pixels = Array{Complex{T},3}(undef, Nᵐ, Nᵗ, Nᵈ)
        for k ∈ 1:Nᵈ
            s = Scri.spin_weight(Val(components[k]))
            ₛY = ₛ𝐘(s, ℓₘₐₓ, T′, R′ₚ)
            if s == 0
                pixels[:, :, k] = ₛY * d′[:, :, k]
            else
                F = qr(ₛY)
                Q = F.Q * Matrix{Bool}(LinearAlgebra.I, Nᵐ, Nᵐ)
                Q⊥ = Q[:, (Nᵐ - s ^ 2 + 1):end]
                pixels[:, :, k] = ₛY * d′[(s ^ 2 + 1):end, :, k] + Q⊥ * d′[1:(s ^ 2), :, k]
            end
        end
        return pixels
    end
end

@testitem "transform_objective: matches exact pixel reconstruction of transform!" tags = [
    :unit, :validation
] setup = [PixelReconstruction] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # The objective must equal Σ|pixels − target|², where the pixels are reconstructed
    # exactly from `transform!` output — same transformation, two very different code
    # paths (fused per-pixel dots vs GEMM synthesis + augmented-LU analysis roundtrip).
    rng = Random.Xoshiro(4242)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 24
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    v⃗ = QuatVec(1e-2, -2e-2, 3e-2)
    R = randn(rng, Rotor{Float64})
    α = 0.1 * Random.randn(rng, ComplexF64, 9)
    dc = Scri.DataComponents(:h, :ψ₄, :ψ₃, :ψ₂)  # includes s = ∓2, ∓1, and 0
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 4)

    d′, t′ = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    pixels = PixelReconstruction.reconstruct_pixels(d′, dc)
    target = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 4), t, dc, t′)
    expected = sum(abs2, pixels .- target)

    obj = Scri.transform_objective(data, t, v⃗, R, α, dc, target; t′=t′)
    @test obj ≈ expected rtol = 1e-11

    # `data` is read-only for the objective.
    data₀ = copy(data)
    Scri.transform_objective(data, t, v⃗, R, α, dc, target; t′=t′)
    @test data == data₀
end

@testitem "transform_objective: ForwardDiff derivatives match finite differences" tags = [
    :unit, :validation
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import ForwardDiff
    import Random

    # Clones of the `transform!` AD closures, but with no dual-typed copy of `data`
    # inside — the whole point of `transform_objective` is that `data` stays primal.
    rng = Random.Xoshiro(7171)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 20
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    R₀ = randn(rng, Rotor{Float64})
    dc = Scri.DataComponents(:h, :ψ₄)
    data₀ = Random.randn(rng, ComplexF64, N, Nᵗ, 2)
    data₀[1:4, :, :] .= 0

    _, t′₀ = Scri.transform!(
        copy(data₀), copy(t), QuatVec(0.0, 0.0, 1e-3), R₀, zeros(ComplexF64, 4), dc
    )
    tc = (t′₀[begin] + t′₀[end]) / 2
    t′fix = tc .+ 0.9 .* (t′₀ .- tc)
    target = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′fix)

    function by_boost(β::T) where {T}
        return Scri.transform_objective(
            data₀,
            t,
            QuatVec(zero(β), zero(β), β),
            R₀,
            zeros(Complex{T}, 4),
            dc,
            target;
            t′=t′fix,
        )
    end
    function by_δt(δt::T) where {T}
        α = zeros(Complex{T}, 4)
        α[1] = 2 * √T(π) * δt
        return Scri.transform_objective(
            data₀, t, QuatVec(0.0, 0.0, 1e-3), R₀, α, dc, target; t′=t′fix
        )
    end
    function by_rotation(θ::T) where {T}
        Rθ = Rotor{T}(R₀) * Rotor(cos(θ / 2), sin(θ / 2), zero(θ), zero(θ))
        return Scri.transform_objective(
            data₀,
            t,
            QuatVec(zero(θ), zero(θ), zero(θ) + 1e-3),
            Rθ,
            zeros(Complex{T}, 4),
            dc,
            target;
            t′=t′fix,
        )
    end

    h = 1e-6
    for (f, x₀) ∈ ((by_boost, 1e-3), (by_δt, 0.5), (by_rotation, 0.1))
        d = ForwardDiff.derivative(f, x₀)
        fd = (f(x₀ + h) - f(x₀ - h)) / 2h
        @test isfinite(d)
        @test abs(d - fd) < 1e-6 * max(1.0, abs(fd))
    end
end

@testitem "transform_objective: derivatives match the memory-heavy dual reference" tags = [
    :unit, :validation
] setup = [PixelReconstruction] begin
    using Quaternionic: QuatVec, Rotor, randn
    import ForwardDiff
    import Random

    # The reference computes the same objective the memory-heavy way: dual-typed copy of
    # `data` through `transform!`, exact pixel reconstruction (primal matrices × dual
    # modes), then the same L² sum.  Same function, two computations — values and
    # derivatives should agree far beyond finite-difference accuracy.
    rng = Random.Xoshiro(7272)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 20
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    R₀ = randn(rng, Rotor{Float64})
    dc = Scri.DataComponents(:h, :ψ₄)
    data₀ = Random.randn(rng, ComplexF64, N, Nᵗ, 2)
    data₀[1:4, :, :] .= 0

    _, t′₀ = Scri.transform!(
        copy(data₀), copy(t), QuatVec(0.0, 0.0, 1e-3), R₀, zeros(ComplexF64, 4), dc
    )
    tc = (t′₀[begin] + t′₀[end]) / 2
    t′fix = tc .+ 0.9 .* (t′₀ .- tc)
    target = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′fix)

    function heavy(β::T) where {T}
        data = Complex{T}.(data₀)
        d, _ = Scri.transform!(
            data,
            copy(t),
            QuatVec(zero(β), zero(β), β),
            R₀,
            zeros(Complex{T}, 4),
            dc;
            t′=t′fix,
        )
        pixels = PixelReconstruction.reconstruct_pixels(d, dc)
        return sum(abs2, pixels .- target)
    end
    function lean(β::T) where {T}
        return Scri.transform_objective(
            data₀,
            t,
            QuatVec(zero(β), zero(β), β),
            R₀,
            zeros(Complex{T}, 4),
            dc,
            target;
            t′=t′fix,
        )
    end

    β₀ = 1e-3
    @test lean(β₀) ≈ heavy(β₀) rtol = 1e-11
    @test ForwardDiff.derivative(lean, β₀) ≈ ForwardDiff.derivative(heavy, β₀) rtol = 1e-9
end

@testitem "sYlm_values!: ForwardDiff duals through the Wigner-H recursion" tags = [
    :unit, :fast
] begin
    using Quaternionic: Rotor, randn
    using SphericalFunctions: sYlm_prep, sYlm_values!
    import ForwardDiff
    import Random

    # `transform_objective` sends dual-typed rotors through `sYlm_values!`; the Wigner-H
    # recursion inside is generic (its `@fastmath` sites fall back to ordinary ops for
    # duals).  Check d/dθ of the ₛYₗₘ row under R(θ) = exp(θ𝐢/2)·R₀ against central
    # finite differences, for every spin weight used by the data components, through
    # random linear functionals of the row (covering real and imaginary parts).
    rng = Random.Xoshiro(7373)
    ℓₘₐₓ = 5
    R₀ = randn(rng, Rotor{Float64})
    c = Random.randn(rng, ComplexF64, (ℓₘₐₓ + 1)^2)
    h = 1e-6
    for s ∈ -2:2
        function Yrow(θ::T) where {T}
            storage = sYlm_prep(ℓₘₐₓ, 2, T)
            Rθ = Rotor(cos(θ / 2), sin(θ / 2), zero(θ), zero(θ)) * R₀
            return copy(sYlm_values!(storage, Rθ, s))
        end
        for f ∈ (θ -> real(sum(c .* Yrow(θ))), θ -> imag(sum(c .* Yrow(θ))))
            d = ForwardDiff.derivative(f, 0.1)
            fd = (f(0.1 + h) - f(0.1 - h)) / 2h
            @test isfinite(d)
            @test abs(d - fd) < 1e-6 * max(1.0, abs(fd))
        end
    end
end

@testitem "transform_objective: declared input band limit ℓₘₐₓ₀ is exact and validated" tags = [
    :unit, :fast
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # Mirrors the `transform!` ℓₘₐₓ₀ testitem: on zero-padded data the declaration only
    # skips provably-zero terms in the per-pixel dots, and a false declaration throws.
    rng = Random.Xoshiro(7474)
    ℓ₀ = 2
    ℓ = 6
    N = (ℓ + 1)^2
    Nᵗ = 20
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    v⃗ = QuatVec(1e-3, -2e-3, 0.0)
    R = randn(rng, Rotor{Float64})
    α = 1e-2 * Random.randn(rng, ComplexF64, 4)
    dc = Scri.DataComponents(:h, :ψ₄)

    data = zeros(ComplexF64, N, Nᵗ, 2)
    data[5:((ℓ₀ + 1) ^ 2), :, :] .= Random.randn(rng, ComplexF64, (ℓ₀ + 1)^2 - 4, Nᵗ, 2)

    _, t′ = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    target = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′)

    obj = Scri.transform_objective(data, t, v⃗, R, α, dc, target; t′=t′)
    obj₀ = Scri.transform_objective(data, t, v⃗, R, α, dc, target; t′=t′, ℓₘₐₓ₀=ℓ₀)
    @test obj₀ ≈ obj rtol = 1e-13

    # Out-of-range and violated declarations throw.
    @test_throws ArgumentError Scri.transform_objective(
        data, t, v⃗, R, α, dc, target; t′=t′, ℓₘₐₓ₀=ℓ + 1
    )
    bad = copy(data)
    bad[(ℓ₀ + 1) ^ 2 + 3, 1, 1] = 1.0
    @test_throws ArgumentError Scri.transform_objective(
        bad, t, v⃗, R, α, dc, target; t′=t′, ℓₘₐₓ₀=ℓ₀
    )
end

@testitem "transform_objective: determinism, weights, and input validation" tags = [
    :unit, :fast
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    rng = Random.Xoshiro(7575)
    ℓ = 2
    N = (ℓ + 1)^2
    Nᵗ = 12
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    v⃗ = QuatVec(1e-2, 0.0, -1e-2)
    R = randn(rng, Rotor{Float64})
    α = 1e-2 * Random.randn(rng, ComplexF64, 4)
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 2)

    _, t′ = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    target = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′)
    function obj(; weights=nothing, t′grid=t′)
        return Scri.transform_objective(data, t, v⃗, R, α, dc, target; t′=t′grid, weights)
    end

    # Deterministic: repeated calls are bitwise identical, and unit weights match the
    # default exactly.
    @test obj() == obj()
    @test obj(; weights=ones(N)) == obj()

    # Nontrivial weights match the manually weighted sum of per-pixel objectives.
    w = rand(rng, N)
    onehot(i) = [j == i ? 1.0 : 0.0 for j ∈ 1:N]
    per_pixel = [obj(; weights=onehot(i)) for i ∈ 1:N]
    @test obj(; weights=w) ≈ sum(w .* per_pixel) rtol = 1e-13

    # Bad shapes and grids throw.
    @test_throws ArgumentError obj(; weights=ones(N - 1))
    @test_throws ArgumentError Scri.transform_objective(
        data, t, v⃗, R, α, dc, target[:, :, 1:1]; t′=t′
    )
    @test_throws ArgumentError Scri.transform_objective(
        data, t, v⃗, R, α, dc, target[:, 1:(Nᵗ - 1), :]; t′=t′
    )
    @test_throws ArgumentError obj(; t′grid=reverse(t′))
    @test_throws ArgumentError Scri.transform_objective(
        data, t, v⃗, R, α, dc, target[:, 1:0, :]; t′=t′[1:0]
    )
    @test_throws ArgumentError Scri.transform_objective(
        data, t, v⃗, R, α, dc, target; t′=t′ .+ 1000.0
    )
end

@testitem "pixel_waveform: direct synthesis and interpolation-free grids" tags = [
    :unit, :fast
] begin
    using SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors
    import Random

    rng = Random.Xoshiro(7676)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 16
    t = collect(LinRange(-8.0, 8.0, Nᵗ))
    dc = Scri.DataComponents(:h, :ψ₄, :ψ₃, :ψ₂)  # spins -2, -2, -1, 0
    R′ₚ = golden_ratio_spiral_rotors(0, ℓ, Float64)

    # Time-constant modes: splining a constant is exact, so every output slice must be
    # the direct synthesis ₛ𝐘·modes, at every t′ — including a coarse, offset grid.
    modes = Random.randn(rng, ComplexF64, N, 4)
    data = repeat(reshape(modes, N, 1, 4), 1, Nᵗ, 1)
    t′ = collect(LinRange(-7.3, 6.1, 5))
    out = Scri.pixel_waveform(data, t, dc, t′)
    @test size(out) == (N, 5, 4)
    for (k, s) ∈ enumerate((-2, -2, -1, 0))
        direct = ₛ𝐘(s, ℓ, Float64, R′ₚ) * modes[(s ^ 2 + 1):end, k]
        for j ∈ 1:5
            @test out[:, j, k] ≈ direct rtol = 1e-13
        end
    end

    # t′ = t: interpolation lands exactly on the knots, so the result is just the
    # synthesis of each time slice (roundoff only).
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 4)
    out = Scri.pixel_waveform(data, t, dc, t)
    for (k, s) ∈ enumerate((-2, -2, -1, 0))
        direct = ₛ𝐘(s, ℓ, Float64, R′ₚ) * data[(s ^ 2 + 1):end, :, k]
        @test out[:, :, k] ≈ direct rtol = 1e-13
    end

    # Grids outside the span of `t`, or non-increasing, or empty, are rejected.
    @test_throws ArgumentError Scri.pixel_waveform(data, t, dc, t .- 1.0)
    @test_throws ArgumentError Scri.pixel_waveform(data, t, dc, reverse(t))
    @test_throws ArgumentError Scri.pixel_waveform(data, t, dc, Float64[])
end

@testitem "transform_objective: coarse output grids (Nᵗ′ ≠ Nᵗ)" tags = [:unit, :validation] setup = [
    PixelReconstruction
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # The kernel evaluates each output time independently, so the objective on a coarse
    # subset of a valid grid must equal the subset-sum of the fine-grid per-slice
    # contributions — while `transform!` itself still requires length(t′) == length(t).
    rng = Random.Xoshiro(7777)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 24
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    v⃗ = QuatVec(1e-2, 2e-2, -1e-2)
    R = randn(rng, Rotor{Float64})
    α = 0.1 * Random.randn(rng, ComplexF64, 4)
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 2)

    d′, t′f = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    pixels = PixelReconstruction.reconstruct_pixels(d′, dc)
    target_f = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′f)

    sub = 1:3:Nᵗ
    obj = Scri.transform_objective(data, t, v⃗, R, α, dc, target_f[:, sub, :]; t′=t′f[sub])
    expected = sum(abs2, pixels[:, sub, :] .- target_f[:, sub, :])
    @test obj ≈ expected rtol = 1e-11

    # A single-sample grid also works.
    obj₁ = Scri.transform_objective(
        data, t, v⃗, R, α, dc, target_f[:, 12:12, :]; t′=t′f[12:12]
    )
    @test obj₁ ≈ sum(abs2, pixels[:, 12, :] .- target_f[:, 12, :]) rtol = 1e-11

    # `transform!` still rejects grids of the wrong length.
    @test_throws ArgumentError Scri.transform!(
        copy(data), copy(t), v⃗, R, copy(α), dc; t′=t′f[sub]
    )
end

@testitem "transform!/transform_objective: generic v⃗/R arguments and dispatch" tags = [
    :unit, :fast
] begin
    using Quaternionic: QuatVec, Rotor, randn
    import Random

    # Plain-vector v⃗ and R must dispatch through the generic conversion methods and give
    # bitwise-identical results — so AD users can pass slices of a parameter vector
    # without touching Quaternionic.
    rng = Random.Xoshiro(7878)
    ℓ = 2
    N = (ℓ + 1)^2
    Nᵗ = 12
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    v⃗ = QuatVec(1e-2, -1e-2, 2e-2)
    R = randn(rng, Rotor{Float64})
    v⃗vec = [v⃗.x, v⃗.y, v⃗.z]
    Rvec = [R.w, R.x, R.y, R.z]
    α = 1e-2 * Random.randn(rng, ComplexF64, 4)
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 2)

    d₁, t′ = Scri.transform!(copy(data), copy(t), v⃗, R, copy(α), dc)
    d₂, t′₂ = Scri.transform!(copy(data), copy(t), v⃗vec, Rvec, copy(α), dc)
    @test d₁ == d₂
    @test t′ == t′₂

    target = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′)
    obj = Scri.transform_objective(data, t, v⃗, R, α, dc, target; t′=t′)
    @test Scri.transform_objective(data, t, v⃗vec, Rvec, α, dc, target; t′=t′) == obj
    @test Scri.transform_objective(data, t, v⃗vec, R, α, dc, target; t′=t′) == obj

    # The BMS convenience method matches the unpacked call, exactly as for `transform!`.
    # The element's accessors describe a (slightly) different transformation than the raw
    # (v⃗, R) above, so use a worst-case grid valid for both, and a matching target.
    g = Scri.BMS{Float64}(; boost_velocity=v⃗, frame_rotation=R)
    t′w = first(Scri.compute_t′(t, 0.1, 0.5))
    target_w = Scri.pixel_waveform(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′w)
    expected = Scri.transform_objective(
        data,
        t,
        Scri.boost_velocity(g),
        Scri.frame_rotation(g),
        Scri.supertranslation(g),
        dc,
        target_w;
        t′=t′w,
    )
    @test Scri.transform_objective(data, t, g, dc, target_w; t′=t′w) == expected

    # No method ambiguities are introduced by the convenience methods.  (A blanket
    # `detect_ambiguities(Scri)` would also pick up a pre-existing, unrelated clash
    # between ForwardDiff's `convert(::Type{Dual}, x)` and the `Signs` conversions
    # whenever ForwardDiff happens to be loaded, so check just these functions.)
    for f ∈ (Scri.transform!, Scri.transform_objective, Scri.pixel_waveform)
        ms = collect(methods(f))
        for i ∈ eachindex(ms), j ∈ 1:(i - 1)
            @test !Base.isambiguous(ms[i], ms[j])
        end
    end
end
