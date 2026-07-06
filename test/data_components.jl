# Tests for DataComponents (src/data_components.jl)

# ── Construction ──────────────────────────────────────────────────────────────

@testitem "DataComponents: valid sets are accepted" tags = [:unit, :fast] begin
    import Scri: DataComponents

    # Stand-alone components (no hierarchy dependency below them).
    @test DataComponents(:ψ₄) isa DataComponents
    @test DataComponents(:σ) isa DataComponents
    @test DataComponents(:h) isa DataComponents
    @test DataComponents(:News) isa DataComponents

    # Ordered Weyl chains from the base up.
    @test DataComponents(:ψ₄, :ψ₃) isa DataComponents
    @test DataComponents(:ψ₄, :ψ₃, :ψ₂) isa DataComponents
    @test DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁) isa DataComponents
    @test DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀) isa DataComponents

    # Mixed sets.
    @test DataComponents(:ψ₄, :σ, :h, :News) isa DataComponents
    @test DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ) isa DataComponents
    @test DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News) isa DataComponents

    # Type parameter encodes the exact tuple.
    dc = DataComponents(:ψ₄, :σ)
    @test dc isa DataComponents{(:ψ₄, :σ)}
end

@testitem "DataComponents: carries its Conventions" tags = [:unit, :fast] begin
    import Scri: DataComponents, Conventions

    # Defaults to the package-native conventions, stored as compile-time singletons.
    dc = DataComponents(:ψ₄, :σ)
    @test dc.conventions === Conventions()

    # A non-default convention travels with the descriptor (and its type).
    X = Conventions(:MB)
    dcX = DataComponents(:ψ₄, :σ; conventions=X)
    @test dcX.conventions === X
    @test dcX isa DataComponents{(:ψ₄, :σ),1,typeof(X)}

    # The string-parsing constructor threads it through too.
    @test DataComponents("psi4", "sigma"; conventions=X).conventions === X
end

@testitem "DataComponents: hierarchy violations are rejected" tags = [:unit, :fast] begin
    import Scri: DataComponents

    # ψₙ (n < 4) without its required upper neighbours.
    @test_throws AssertionError DataComponents(:ψ₃)           # missing ψ₄
    @test_throws AssertionError DataComponents(:ψ₂)           # missing ψ₃, ψ₄
    @test_throws AssertionError DataComponents(:ψ₁)           # missing ψ₂..ψ₄
    @test_throws AssertionError DataComponents(:ψ₀)           # missing ψ₁..ψ₄

    # Gaps in the chain.
    @test_throws AssertionError DataComponents(:ψ₂, :ψ₄)           # missing ψ₃
    @test_throws AssertionError DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₄) # missing ψ₃
end

@testitem "DataComponents: invalid symbols are rejected" tags = [:unit, :fast] begin
    import Scri: DataComponents

    @test_throws AssertionError DataComponents(:foo)
    @test_throws AssertionError DataComponents(:Weyl₀)
    @test_throws AssertionError DataComponents(:ψ₄, :bad_component)
end

# ── Accessors ─────────────────────────────────────────────────────────────────

@testitem "component_index: returns correct 1-based position" tags = [:unit, :fast] begin
    import Scri: DataComponents

    # Weyl-only, ψ₄-first ordering.
    dc = DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀)
    @test Scri.component_index(dc, Val(:ψ₄)) == 1
    @test Scri.component_index(dc, Val(:ψ₃)) == 2
    @test Scri.component_index(dc, Val(:ψ₂)) == 3
    @test Scri.component_index(dc, Val(:ψ₁)) == 4
    @test Scri.component_index(dc, Val(:ψ₀)) == 5

    # Non-Weyl, arbitrary ordering.
    dc2 = DataComponents(:σ, :h, :News)
    @test Scri.component_index(dc2, Val(:σ)) == 1
    @test Scri.component_index(dc2, Val(:h)) == 2
    @test Scri.component_index(dc2, Val(:News)) == 3
end

@testitem "component_index: returns nothing for absent components" tags = [:unit, :fast] begin
    import Scri: DataComponents

    dc = DataComponents(:ψ₄, :ψ₃)
    @test isnothing(Scri.component_index(dc, Val(:ψ₂)))
    @test isnothing(Scri.component_index(dc, Val(:ψ₁)))
    @test isnothing(Scri.component_index(dc, Val(:ψ₀)))
    @test isnothing(Scri.component_index(dc, Val(:σ)))
    @test isnothing(Scri.component_index(dc, Val(:h)))
    @test isnothing(Scri.component_index(dc, Val(:News)))
end

@testitem "has_component: consistent with component_index" tags = [:unit, :fast] begin
    import Scri: DataComponents

    dc = DataComponents(:ψ₄, :ψ₃, :ψ₂, :σ)
    for s ∈ Scri.ValidDataComponents
        present = !isnothing(Scri.component_index(dc, Val(s)))
        @test Scri.has_component(dc, Val(s)) == present
    end
end

@testitem "ncomponents: matches constructor arity" tags = [:unit, :fast] begin
    import Scri: DataComponents

    @test Scri.ncomponents(DataComponents(:ψ₄)) == 1
    @test Scri.ncomponents(DataComponents(:ψ₄, :ψ₃)) == 2
    @test Scri.ncomponents(DataComponents(:ψ₄, :ψ₃, :ψ₂)) == 3
    @test Scri.ncomponents(DataComponents(:ψ₄, :σ, :h, :News)) == 4
    @test Scri.ncomponents(DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News)) == 8
end

# ── spin_weight and conformal_weight ──────────────────────────────────────────

@testitem "spin_weight: known values for all components" tags = [:unit, :fast] begin
    @test Scri.spin_weight(Val(:ψ₀)) == 2
    @test Scri.spin_weight(Val(:ψ₁)) == 1
    @test Scri.spin_weight(Val(:ψ₂)) == 0
    @test Scri.spin_weight(Val(:ψ₃)) == -1
    @test Scri.spin_weight(Val(:ψ₄)) == -2
    @test Scri.spin_weight(Val(:σ)) == 2
    @test Scri.spin_weight(Val(:h)) == -2
    @test Scri.spin_weight(Val(:News)) == -2
end

@testitem "conformal_weight: known values for all components" tags = [:unit, :fast] begin
    # All five Weyl components share conformal weight −3.
    for s ∈ (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄)
        @test Scri.conformal_weight(Val(s)) == -3
    end
    # Shear and strain have conformal weight −1; News has −2.
    @test Scri.conformal_weight(Val(:σ)) == -1
    @test Scri.conformal_weight(Val(:h)) == -1
    @test Scri.conformal_weight(Val(:News)) == -2
end

# ── mix_components! ───────────────────────────────────────────────────────────

@testitem "mix_components!: identity (κ⁻¹=1, ðt′╱2κ=0, ð²α=0)" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    rng = Random.Xoshiro(42)
    dc = DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News)
    for _ ∈ 1:5
        data = randn(rng, ComplexF64, 8)
        orig = copy(data)
        Scri.mix_components!(data, 1.0, 0.0 + 0im, 0.0 + 0im, dc)
        @test data == orig
    end
end

@testitem "mix_components!: pure conformal scaling (ðt′╱2κ=0, ð²α=0)" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    # When ðt′╱2κ = 0 (e.g., at u = 0 for a boost with no supertranslation, since
    # ðt′╱2κ = −(ðκ/2κ)·u), each component scales by κ^(conformal_weight).
    # Weyl: κ⁻³; σ,h: κ⁻¹; News: κ⁻².
    rng = Random.Xoshiro(7)
    dc = DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News)
    for _ ∈ 1:5
        data = randn(rng, ComplexF64, 8)
        orig = copy(data)
        κ⁻¹ = 0.4 + 0.3 * randn(rng)
        Scri.mix_components!(data, κ⁻¹, 0.0 + 0im, 0.0 + 0im, dc)
        for s ∈ (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄)
            i = Scri.component_index(dc, Val(s))
            @test data[i] ≈ κ⁻¹^3 * orig[i]
        end
        @test data[Scri.component_index(dc, Val(:σ))] ≈
            κ⁻¹ * orig[Scri.component_index(dc, Val(:σ))]
        @test data[Scri.component_index(dc, Val(:h))] ≈
            κ⁻¹ * orig[Scri.component_index(dc, Val(:h))]
        @test data[Scri.component_index(dc, Val(:News))] ≈
            κ⁻¹^2 * orig[Scri.component_index(dc, Val(:News))]
    end
end

@testitem "mix_components!: ψ₄-seed propagation (κ⁻¹=1)" tags = [:unit, :fast, :validation] begin
    import Random
    import Scri: DataComponents

    # When only ψ₄ = z is non-zero and κ⁻¹=1, the lower Weyl components receive
    # the values ψₙ' = (ðt′╱2κ)^(4−n) · z — purely from the nested polynomial.
    rng = Random.Xoshiro(11)
    dc = DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀)
    for _ ∈ 1:8
        z = randn(rng, ComplexF64)
        f = randn(rng, ComplexF64)   # ðt′╱2κ
        data = ComplexF64[z, 0, 0, 0, 0]   # ψ₄=z, ψ₃=ψ₂=ψ₁=ψ₀=0
        Scri.mix_components!(data, 1.0, f, 0.0 + 0im, dc)
        for (s, exp) ∈ ((:ψ₄, 0), (:ψ₃, 1), (:ψ₂, 2), (:ψ₁, 3), (:ψ₀, 4))
            i = Scri.component_index(dc, Val(s))
            @test data[i] ≈ f^exp * z atol = 4eps(Float64) * abs(f)^exp * abs(z)
        end
    end
end

@testitem "mix_components!: all-ones Weyl gives binomial pattern" tags = [
    :unit, :fast, :validation
] begin
    import Scri: DataComponents

    # When all five ψ inputs equal 1, each output is κ⁻³·(1+ðt′╱2κ)^(4−n).
    # This follows from the binomial expansion of (1 + ðt′╱2κ * ∂_u)^4 acting on 1.
    dc = DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀)
    κ⁻¹ = 2.0
    f = 3.0 + 2.0im
    data = ones(ComplexF64, 5)
    Scri.mix_components!(data, κ⁻¹, f, 0.0 + 0im, dc)
    for (s, exp) ∈ ((:ψ₄, 0), (:ψ₃, 1), (:ψ₂, 2), (:ψ₁, 3), (:ψ₀, 4))
        i = Scri.component_index(dc, Val(s))
        @test data[i] ≈ κ⁻¹^3 * (1 + f)^exp
    end
end

@testitem "mix_components!: σ and h shift by ð²α and its conjugate" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    # σ' = κ⁻¹·(σ + ½ð²α),  h' = κ⁻¹·(h + ½conj(ð²α))
    rng = Random.Xoshiro(99)
    dc = DataComponents(:σ, :h)
    for _ ∈ 1:8
        σ_v = randn(rng, ComplexF64)
        h_v = randn(rng, ComplexF64)
        κ⁻¹ = 0.5 + randn(rng)
        ð²α = randn(rng, ComplexF64)
        data = ComplexF64[σ_v, h_v]
        Scri.mix_components!(data, κ⁻¹, 0.0 + 0im, ð²α, dc)
        @test data[1] ≈ κ⁻¹ * (σ_v + ð²α / 2)
        @test data[2] ≈ κ⁻¹ * (h_v + conj(ð²α) / 2)
    end
end

@testitem "mix_components!: convention factors appear as in the formulas" tags = [
    :unit, :fast, :validation
] begin
    import Random
    import Scri: DataComponents, Conventions, dyad_factor, shear_factor, strain_factor

    # With a generic convention carried by `dc`, the laws read (docs, "Convention
    # dependence"): towers on c_l c_m·ðu′/2κ at ℐ⁺ and conj(ðt′/2κ)/(c_l c_m) at ℐ⁻;
    # shifts F_σ·ð²α/2 and F_h·ð̄²α/2.
    rng = Random.Xoshiro(17)
    X = Conventions(; c_s=-1, c_ψ=-1, c_σ=-1, c_l=(-√2), c_m=cis(π / 4), c_h=2)
    q = dyad_factor(X)
    for _ ∈ 1:4
        ψs = randn(rng, ComplexF64, 2)
        σ_v, h_v = randn(rng, ComplexF64), randn(rng, ComplexF64)
        κ⁻¹ = 0.5 + randn(rng)
        f = randn(rng, ComplexF64)     # ðt′╱2κ
        ð²α = randn(rng, ComplexF64)

        # ℐ⁺: ψ₃′ = κ⁻³(ψ₃ + (c_l c_m)·ðu′/2κ·ψ₄), σ/h shifts by F_σ, F_h.
        dc⁺ = DataComponents(:ψ₄, :ψ₃, :σ, :h; conventions=X)
        data = ComplexF64[ψs[1], ψs[2], σ_v, h_v]
        Scri.mix_components!(data, κ⁻¹, f, ð²α, dc⁺)
        @test data[2] ≈ κ⁻¹^3 * (ψs[2] + q * f * ψs[1])
        @test data[3] ≈ κ⁻¹ * (σ_v + shear_factor(X, +1) * ð²α / 2)
        @test data[4] ≈ κ⁻¹ * (h_v + strain_factor(X) * conj(ð²α) / 2)

        # ℐ⁻: the tower mixes downward on conj(ðt′/2κ)/(c_l c_m), and the shifts negate.
        dc⁻ = DataComponents(:ψ₀, :ψ₁, :σ, :h; ℐ=-1, conventions=X)
        data = ComplexF64[ψs[1], ψs[2], σ_v, h_v]
        Scri.mix_components!(data, κ⁻¹, f, ð²α, dc⁻)
        @test data[2] ≈ κ⁻¹^3 * (ψs[2] + (conj(f) / q) * ψs[1])
        @test data[3] ≈ κ⁻¹ * (σ_v - shear_factor(X, -1) * ð²α / 2)
        @test data[4] ≈ κ⁻¹ * (h_v - strain_factor(X) * conj(ð²α) / 2)
    end
end

# ── represent! ────────────────────────────────────────────────────────────────

@testitem "represent!: identity, factors, and round trips" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents, Conventions, One, conversion_factor, represent!

    rng = Random.Xoshiro(3)
    comps = (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News, :φ₀, :φ₁, :φ₂)
    X = Conventions(; c_s=-1, c_ψ=-1, c_σ=-1, c_φ=-1, c_l=(-√2), c_m=cis(π / 4), c_h=2)
    Y = Conventions(:MB)

    for ℐ ∈ (+1, -1)
        dcS = DataComponents(comps...; ℐ)
        data = randn(rng, ComplexF64, 9, 3, length(comps))
        orig = copy(data)

        # SXS → SXS is an exact no-op (every factor ratio is One()).
        d, dc′ = represent!(copy(data), dcS, Conventions())
        @test d == orig
        @test dc′ === dcS

        # SXS → X multiplies each component slice by its conversion factor.
        d, dcX = represent!(copy(data), dcS, X)
        @test dcX.conventions === X
        for (k, S) ∈ enumerate(comps)
            @test d[:, :, k] ≈ conversion_factor(X, Val(S), ℐ) .* orig[:, :, k]
        end

        # Round trip X → Y → X is the identity to roundoff.
        dY, dcY = represent!(copy(d), dcX, Y)
        dX, _ = represent!(dY, dcY, X)
        @test dX ≈ d rtol = 4eps(Float64)
    end

    # Component-count mismatch is caught.
    dc2 = DataComponents(:ψ₄, :h)
    @test_throws AssertionError represent!(zeros(ComplexF64, 4, 2, 3), dc2, X)
end

@testitem "mix_components!: News scales by κ⁻² with no mixing" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    rng = Random.Xoshiro(13)
    dc = DataComponents(:News)
    for _ ∈ 1:8
        news = randn(rng, ComplexF64)
        κ⁻¹ = 0.5 + randn(rng)
        data = ComplexF64[news]
        Scri.mix_components!(data, κ⁻¹, randn(rng, ComplexF64), randn(rng, ComplexF64), dc)
        @test data[1] ≈ κ⁻¹^2 * news
    end
end
