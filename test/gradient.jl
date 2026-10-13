# Tests for the analytic gradient of the objective (src/gradient.jl)

@testitem "spline_transpose!: the transpose of the kernel's interpolation" tags = [
    :unit, :fast
] begin
    using LinearAlgebra: dot
    import Random

    # ⟨y, S x⟩ = ⟨Sᵀ y, x⟩, where S is interpolation exactly as `transform_pixel!` does
    # it (with mixing that changes nothing), on unevenly spaced knots, including the
    # smallest grid of four knots, times on knots, and a time a hair below the first knot.
    rng = Random.Xoshiro(6161)
    dc = Scri.DataComponents(:ψ₂, :ψ₃, :ψ₄)
    for Nᵗ ∈ (4, 9, 30)
        t = cumsum(0.2 .+ rand(rng, Nᵗ))
        cache = Scri.CubicSplineCache(t)
        τ = [t[1] - 1e-13, t[1], (t[2] + t[3]) / 2, t[end - 1], t[end] - 0.01, t[end]]
        function S(x)
            d = zeros(ComplexF64, 3, Nᵗ)
            d[1, :] = x
            out = zeros(ComplexF64, 3, length(τ))
            Scri.transform_pixel!(
                out, d, similar(d), cache, t, τ, 1.0, 0.0, 0.0im, 0.0im, 0.0im, dc
            )
            return out[1, :]
        end
        for _ ∈ 1:3
            x = Random.randn(rng, ComplexF64, Nᵗ)
            y = Random.randn(rng, ComplexF64, length(τ))
            Sᵀy = zeros(ComplexF64, Nᵗ)
            Scri.spline_transpose!(Sᵀy, y, cache, t, τ)
            @test dot(y, S(x)) ≈ dot(Sᵀy, x) rtol = 1e-13
        end
    end
end

@testitem "mixing_matrix!: mix_components! as a matrix and an offset" tags = [:unit, :fast] begin
    import Random

    rng = Random.Xoshiro(6262)
    for dc ∈ (
        Scri.DataComponents(:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
        Scri.DataComponents(:σ, :ψ₄, :ψ₃; conventions=Scri.Conventions(:MB)),
        Scri.DataComponents(:λ, :ψ₀, :ψ₁, :ψ₂; ℐ=-1),
        Scri.DataComponents(:News, :φ₂, :φ₁, :φ₀),
    )
        Nᵈ = Scri.ncomponents(dc)
        M = zeros(ComplexF64, Nᵈ, Nᵈ)
        κ⁻¹, ð, ð² = 1.1, 0.2 - 0.1im, 0.3 + 0.05im
        Scri.mixing_matrix!(M, κ⁻¹, ð, ð², dc)
        offset = zeros(ComplexF64, Nᵈ)
        Scri.mix_components!(offset, κ⁻¹, ð, ð², dc)
        for _ ∈ 1:3
            x = Random.randn(rng, ComplexF64, Nᵈ)
            mixed = copy(x)
            Scri.mix_components!(mixed, κ⁻¹, ð, ð², dc)
            @test mixed ≈ M * x + offset rtol = 1e-14
        end
    end
end

@testitem "bms_from_parameters: rapidity and rotation vectors, and the inverse" tags = [
    :unit, :fast
] begin
    using DoubleFloats: Double64
    using Quaternionic: QuatVec, Rotor, absvec
    import ForwardDiff
    import Random

    rng = Random.Xoshiro(6363)
    θ = [0.1, -0.2, 0.3, 0.2, 0.1, -0.15, Random.randn(rng, 9)...]
    v⃗, R, α = Scri.bms_from_parameters(θ)
    𝐛 = θ[4:6]
    @test collect(vec(v⃗)) ≈ tanh(sqrt(sum(abs2, 𝐛))) .* 𝐛 ./ sqrt(sum(abs2, 𝐛)) rtol = 1e-15
    @test R ≈ exp(QuatVec(θ[1:3]...) / 2)
    @test α ≈ Scri.impose_reality(α, 2, 1) rtol = 1e-15  # already real on the sphere
    @test Scri.parameters_from_bms(v⃗, R, α) ≈ θ rtol = 1e-14

    # Zero is the identity, and nothing is singular there or at small rapidities, where
    # the series and the closed form of tanh(√x)/√x meet.
    v₀, R₀, α₀ = Scri.bms_from_parameters(zeros(7))
    @test iszero(absvec(v₀)) && R₀ == one(Rotor{Float64}) && α₀ == [0.0im]
    x₀ = eps(Float64)^(1 / 5)
    for x ∈ (0.0, x₀, 1.0)
        f(y) = Scri.tanhc_sqrt(y)
        @test f(x) ≈ (x == 0 ? 1.0 : tanh(sqrt(x)) / sqrt(x)) rtol = 1e-15
        @test isfinite(ForwardDiff.derivative(f, x))
    end
    @test Scri.tanhc_sqrt(x₀ * (1 - 1e-12)) ≈ Scri.tanhc_sqrt(x₀ * (1 + 1e-12)) rtol = 1e-14

    # Both series give way to the closed forms while they are still exact to roundoff, so
    # each function is as accurate as its type on both sides of the switch.
    for T ∈ (Double64, BigFloat), x ∈ T.((1e-20, 1e-8, 1e-6, 1e-4, 1e-3, 1e-2))
        r = sqrt(big(x))
        @test Scri.tanhc_sqrt(x) ≈ tanh(r) / r rtol = 4eps(T)
        v⃗_T = QuatVec(zero(T), zero(T), T(tanh(r)))
        θ_T = Scri.parameters_from_bms(v⃗_T, one(Rotor{T}), [zero(Complex{T})])
        @test θ_T[6] ≈ r rtol = 10eps(T)
    end
    J = ForwardDiff.jacobian(θ -> collect(vec(Scri.bms_from_parameters(θ)[1])), zeros(7))
    @test J[:, 4:6] ≈ [1 0 0; 0 1 0; 0 0 1]
    @test_throws ArgumentError Scri.bms_from_parameters(zeros(8))
end

@testitem "objective_and_gradient!: matches ForwardDiff through transform_objective" tags = [
    :validation
] begin
    import ForwardDiff
    import Random

    # The value must be the objective's, and the gradient must be ForwardDiff's through
    # `transform_objective`, at the identity and away from it, for components that mix,
    # on both null infinities, with weights and normalization, and for the objective that
    # compares only the modes up to ℓ_fit.
    rng = Random.Xoshiro(6464)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 30
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    t′ = Scri.compute_t′(t, 0.1, 1.0)
    for (dc, kw) ∈ (
        (Scri.DataComponents(:h, :ψ₄, :ψ₃), (;)),
        (Scri.DataComponents(:λ, :ψ₀, :ψ₁; ℐ=-1), (; normalize=true)),
        (
            Scri.DataComponents(:σ, :ψ₄; conventions=Scri.Conventions(:MB)),
            (; time_weights=rand(rng, length(t′)), component_weights=[2.0, 0.5]),
        ),
        (Scri.DataComponents(:h, :ψ₄, :ψ₃), (; ℓ_fit=3, normalize=true)),
        (Scri.DataComponents(:λ, :ψ₀; ℐ=-1), (; ℓ_fit=2, component_weights=[0.5, 2.0])),
    )
        Nᵈ = Scri.ncomponents(dc)
        data = Random.randn(rng, ComplexF64, N, Nᵗ, Nᵈ)
        plan = Scri.ObjectivePlan(
            Random.randn(rng, ComplexF64, N, Nᵗ, Nᵈ), t, dc, t′; kw...
        )
        J(θ) = Scri.transform_objective(data, t, Scri.bms_from_parameters(θ)..., dc, plan)
        for θ ∈ (
            zeros(15),
            [0.1, -0.2, 0.3, 0.02, 0.01, -0.03, (0.05 .* Random.randn(rng, 9))...],
        )
            ∇ = similar(θ)
            value = Scri.objective_and_gradient!(∇, θ, data, t, dc, plan)
            @test value ≈ J(θ) rtol = 1e-13
            @test ∇ ≈ ForwardDiff.gradient(J, θ) rtol = 1e-11
        end
    end

    # And central differences agree, as a check independent of ForwardDiff.
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 2)
    plan = Scri.ObjectivePlan(Random.randn(rng, ComplexF64, N, Nᵗ, 2), t, dc, t′)
    θ = [0.05, 0.1, -0.1, 0.01, -0.02, 0.03, (0.05 .* Random.randn(rng, 4))...]
    ∇ = similar(θ)
    Scri.objective_and_gradient!(∇, θ, data, t, dc, plan)
    J₂(θ) = Scri.transform_objective(data, t, Scri.bms_from_parameters(θ)..., dc, plan)
    h = 1e-6
    for p ∈ eachindex(θ)
        e = [q == p ? h : 0.0 for q ∈ eachindex(θ)]
        @test ∇[p] ≈ (J₂(θ + e) - J₂(θ - e)) / 2h rtol = 1e-6 atol = 1e-6
    end

    # The gradient must fit θ.
    @test_throws ArgumentError Scri.objective_and_gradient!(zeros(3), θ, data, t, dc, plan)
end

@testitem "objective_and_gradient!: Double64" tags = [:validation] begin
    using DoubleFloats: Double64
    import ForwardDiff
    import Random

    rng = Random.Xoshiro(6565)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 16
    t = collect(LinRange(Double64(-10), Double64(10), Nᵗ))
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Complex{Double64}.(Random.randn(rng, ComplexF64, N, Nᵗ, 2))
    target = Complex{Double64}.(Random.randn(rng, ComplexF64, N, Nᵗ, 2))
    plan = Scri.ObjectivePlan(target, t, dc, Scri.compute_t′(t, Double64(0.1), Double64(1)))
    θ = Double64.([0.1, -0.2, 0.3, 0.02, 0.01, -0.03, 0.05, -0.04, 0.03, 0.02])
    J(θ) = Scri.transform_objective(data, t, Scri.bms_from_parameters(θ)..., dc, plan)
    ∇ = similar(θ)
    @test Scri.objective_and_gradient!(∇, θ, data, t, dc, plan) ≈ J(θ) rtol = 1e-28
    @test ∇ ≈ ForwardDiff.gradient(J, θ) rtol = 1e-25
end
