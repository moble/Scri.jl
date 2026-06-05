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

@testitem "mix_components!: identity (k⁻¹=1, ðt′╱k=0, ð²α=0)" tags = [:unit, :fast] begin
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

@testitem "mix_components!: pure conformal scaling (ðt′╱k=0, ð²α=0)" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    # When ðt′╱k = 0 (e.g., at u = 0 for a boost with no supertranslation, since
    # ðt′╱k = −(ðk/k)·u), each component scales by k^(conformal_weight).
    # Weyl: k⁻³; σ,h: k⁻¹; News: k⁻².
    rng = Random.Xoshiro(7)
    dc = DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News)
    for _ ∈ 1:5
        data = randn(rng, ComplexF64, 8)
        orig = copy(data)
        k⁻¹ = 0.4 + 0.3 * randn(rng)
        Scri.mix_components!(data, k⁻¹, 0.0 + 0im, 0.0 + 0im, dc)
        for s ∈ (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄)
            i = Scri.component_index(dc, Val(s))
            @test data[i] ≈ k⁻¹^3 * orig[i]
        end
        @test data[Scri.component_index(dc, Val(:σ))] ≈
            k⁻¹ * orig[Scri.component_index(dc, Val(:σ))]
        @test data[Scri.component_index(dc, Val(:h))] ≈
            k⁻¹ * orig[Scri.component_index(dc, Val(:h))]
        @test data[Scri.component_index(dc, Val(:News))] ≈
            k⁻¹^2 * orig[Scri.component_index(dc, Val(:News))]
    end
end

@testitem "mix_components!: ψ₄-seed propagation (k⁻¹=1)" tags = [:unit, :fast, :validation] begin
    import Random
    import Scri: DataComponents

    # When only ψ₄ = z is non-zero and k⁻¹=1, the lower Weyl components receive
    # the values ψₙ' = (−ðt′╱k)^(4−n) · z — purely from the nested polynomial.
    rng = Random.Xoshiro(11)
    dc = DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀)
    for _ ∈ 1:8
        z = randn(rng, ComplexF64)
        f = randn(rng, ComplexF64)   # ðt′╱k
        data = ComplexF64[z, 0, 0, 0, 0]   # ψ₄=z, ψ₃=ψ₂=ψ₁=ψ₀=0
        Scri.mix_components!(data, 1.0, f, 0.0 + 0im, dc)
        for (s, exp) ∈ ((:ψ₄, 0), (:ψ₃, 1), (:ψ₂, 2), (:ψ₁, 3), (:ψ₀, 4))
            i = Scri.component_index(dc, Val(s))
            @test data[i] ≈ (-f)^exp * z atol = 4eps(Float64) * abs(f)^exp * abs(z)
        end
    end
end

@testitem "mix_components!: all-ones Weyl gives binomial pattern" tags = [
    :unit, :fast, :validation
] begin
    import Scri: DataComponents

    # When all five ψ inputs equal 1, each output is k⁻³·(1−ðt′╱k)^(4−n).
    # This follows from the binomial expansion of (1 − ðt′╱k * ∂_u)^4 acting on 1.
    dc = DataComponents(:ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀)
    k⁻¹ = 2.0
    f = 3.0 + 2.0im
    data = ones(ComplexF64, 5)
    Scri.mix_components!(data, k⁻¹, f, 0.0 + 0im, dc)
    for (s, exp) ∈ ((:ψ₄, 0), (:ψ₃, 1), (:ψ₂, 2), (:ψ₁, 3), (:ψ₀, 4))
        i = Scri.component_index(dc, Val(s))
        @test data[i] ≈ k⁻¹^3 * (1 - f)^exp
    end
end

@testitem "mix_components!: σ and h shift by ð²α and its conjugate" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    # σ' = k⁻¹·(σ + ð²α),  h' = k⁻¹·(h + conj(ð²α))
    rng = Random.Xoshiro(99)
    dc = DataComponents(:σ, :h)
    for _ ∈ 1:8
        σ_v = randn(rng, ComplexF64)
        h_v = randn(rng, ComplexF64)
        k⁻¹ = 0.5 + randn(rng)
        ð²α = randn(rng, ComplexF64)
        data = ComplexF64[σ_v, h_v]
        Scri.mix_components!(data, k⁻¹, 0.0 + 0im, ð²α, dc)
        @test data[1] ≈ k⁻¹ * (σ_v + ð²α)
        @test data[2] ≈ k⁻¹ * (h_v + conj(ð²α))
    end
end

@testitem "mix_components!: News scales by k⁻² with no mixing" tags = [:unit, :fast] begin
    import Random
    import Scri: DataComponents

    rng = Random.Xoshiro(13)
    dc = DataComponents(:News)
    for _ ∈ 1:8
        news = randn(rng, ComplexF64)
        k⁻¹ = 0.5 + randn(rng)
        data = ComplexF64[news]
        Scri.mix_components!(data, k⁻¹, randn(rng, ComplexF64), randn(rng, ComplexF64), dc)
        @test data[1] ≈ k⁻¹^2 * news
    end
end
