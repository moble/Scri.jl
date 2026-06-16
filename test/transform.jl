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
    # simply ψ₂′ = k⁻³ ψ₂ with 1/k = γ(1 - εᴵ v⃗⋅n̂).  We reconstruct that pointwise from the
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
        k⁻¹ = [γ * (1 - εᴵ * dot(vec(v⃗), vec(Rₚ[p](𝐤)))) for p ∈ eachindex(Rₚ)]
        expected = @. k⁻¹^3 * ψ₂_rest                            # ψ₂′ = k⁻³ ψ₂ (no mixing)

        # The square s=0 analysis is exactly invertible, so synthesizing the output modes on
        # the B-frame grid recovers the transformed pixel values.
        out_pixels = ₛ𝐘(0, ℓ, Float64, Rs) * data[:, 1, 1]
        @test maximum(abs, out_pixels .- expected) < 1e-12

        # ψ₂′ is constant in time (constant input, time-independent k⁻³), so every slice
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
