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

@testitem "synthesis_transpose_block!: the transpose of synthesize_block!" tags = [
    :unit, :fast
] begin
    using DoubleFloats: Double64
    import Random

    # For harmonics Y and any array ∂ν╱∂σ, the transpose ∂ν╱∂Y satisfies
    # Σ ∂ν╱∂Y Y = Σ ∂ν╱∂σ σ, where σ is synthesized from Y, with no conjugation.
    rng = Random.Xoshiro(3131)
    Nᵐ, Nᵗ, nb = 25, 6, 5
    spins = [-2, 0, -2]
    used = [s - minimum(spins) + 1 for s ∈ unique(spins)]
    for (T, N₀) ∈ ((Float64, Nᵐ), (Float64, 16), (Double64, Nᵐ), (Double64, 16))
        data = Random.randn(rng, Complex{T}, Nᵐ, Nᵗ, length(spins))
        Y = Random.randn(rng, Complex{T}, Nᵐ, nb, 3)
        ∂ν╱∂σ = Random.randn(rng, Complex{T}, Nᵗ, nb, length(spins))
        σ = similar(∂ν╱∂σ)
        Scri.synthesize_block!(σ, data, Y, spins, N₀)
        ∂ν╱∂Y = fill(Complex{T}(NaN), Nᵐ, nb, 3)  # entries that must be left alone
        Scri.synthesis_transpose_block!(∂ν╱∂Y, data, ∂ν╱∂σ, spins, N₀)
        lhs = sum(
            sum(∂ν╱∂Y[m, q, col] * Y[m, q, col] for m ∈ ((s ^ 2 + 1):N₀), q ∈ 1:nb) for
            (s, col) ∈ zip(unique(spins), used)
        )
        @test lhs ≈ sum(∂ν╱∂σ .* σ) rtol = 10eps(T)
        @test all(isnan, ∂ν╱∂Y[:, :, 2])  # the spin weight −1 is unused
        @test all(isnan, ∂ν╱∂Y[1:3, :, 1])  # ℓ < 2 for the spin weight −2
        N₀ < Nᵐ && @test all(isnan, ∂ν╱∂Y[(N₀ + 1):end, :, :])
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
    η⃗ = θ[4:6]
    @test collect(vec(v⃗)) ≈ tanh(sqrt(sum(abs2, η⃗))) .* η⃗ ./ sqrt(sum(abs2, η⃗)) rtol = 1e-15
    @test R ≈ exp(QuatVec(θ[1:3]...) / 2)
    @test α ≈ Scri.impose_reality(α, 2, 1) rtol = 1e-15  # already real on the sphere
    @test Scri.parameters_from_bms(v⃗, R, α) ≈ θ rtol = 1e-14
    @test_throws ArgumentError Scri.parameters_from_bms([0.1], R, α)

    # Either sign of the rotor gives the same, shortest rotation vector, so that `-1` is
    # the identity rather than a rotation by 2π.
    @test Scri.parameters_from_bms(v⃗, -R, α) ≈ θ rtol = 1e-14
    @test Scri.parameters_from_bms(zeros(3), -one(Rotor{Float64}), [0.0im])[1:3] == zeros(3)

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

    # Six parameters are a Lorentz transformation with no supertranslation.
    v₆, R₆, α₆ = Scri.bms_from_parameters(θ[1:6])
    @test v₆ == v⃗ && R₆ == R && α₆ == [0.0im]
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

    # Six parameters fit a Lorentz transformation alone, exactly as seven with a zero time
    # translation do.
    ∇₆, ∇₇ = zeros(6), zeros(7)
    @test Scri.objective_and_gradient!(∇₆, θ[1:6], data, t, dc, plan) ==
        Scri.objective_and_gradient!(∇₇, [θ[1:6]; 0.0], data, t, dc, plan)
    @test ∇₆ == ∇₇[1:6]

    # A boost beyond the speed `t′` was built for is an error, and its message reports the
    # valid range in plain numbers, not the dual numbers the gradient computes with.
    θ_far = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, zeros(4)...]
    err = try
        Scri.objective_and_gradient!(∇, θ_far, data, t, dc, plan)
    catch e
        e
    end
    @test err isa ArgumentError
    @test occursin("exceeds the valid range", err.msg)
    @test !occursin("Dual", err.msg)
end

@testitem "objective_and_gradient!: any number of parameters, in chunks of any width" tags = [
    :validation
] begin
    import ForwardDiff
    import Random

    # The geometry's derivatives are computed for one chunk of the parameters at a time.
    # With 42 parameters there are six chunks of eight, the last with slots to spare, and
    # the gradient must still be ForwardDiff's.  It must also not depend on the width of
    # the chunks: one parameter at a time, a width that does not divide their number, and
    # all of them at once.
    rng = Random.Xoshiro(6666)
    ℓ = 6
    N = (ℓ + 1)^2
    Nᵗ = 12
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    t′ = Scri.compute_t′(t, 0.1, 1.0)
    dc = Scri.DataComponents(:h, :ψ₄, :ψ₃)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 3)
    θ = [0.1, -0.2, 0.3, 0.02, 0.01, -0.03, (0.02 .* Random.randn(rng, 36))...]
    P = length(θ)
    for kw ∈ ((;), (; ℓ_fit=4, normalize=true))
        plan = Scri.ObjectivePlan(Random.randn(rng, ComplexF64, N, Nᵗ, 3), t, dc, t′; kw...)
        J(θ) = Scri.transform_objective(data, t, Scri.bms_from_parameters(θ)..., dc, plan)
        ∇ = similar(θ)
        Φ = Scri.objective_and_gradient!(∇, θ, data, t, dc, plan)
        @test Φ ≈ J(θ) rtol = 1e-13
        @test ∇ ≈ ForwardDiff.gradient(J, θ) rtol = 1e-11
        for W ∈ (1, 5, P)
            Φᵂ, ∇ᵂ = Scri.gradient_core(data, t, θ, dc, plan, nothing, Val(W))
            @test Φᵂ ≈ Φ rtol = 1e-14
            @test ∇ᵂ ≈ ∇ rtol = 1e-13
        end
    end
end

@testitem "objective_and_gradient!: the same with several blocks of pixels" tags = [
    :validation
] begin
    import Random

    # The pixels are taken in blocks, and with ℓ_fit the main pass takes them in reverse
    # order and reuses the harmonics and series that the first pass left from its last
    # block.  A tiny memory budget forces several blocks of the smallest size, four pixels
    # for each thread, with a shorter last block, on a grid fine enough for at least three
    # of them.  The results must be those of the default budget, which makes one block.
    rng = Random.Xoshiro(6868)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 30
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    t′ = Scri.compute_t′(t, 0.1, 1.0)
    dc = Scri.DataComponents(:h, :ψ₄, :ψ₃)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 3)
    target = Random.randn(rng, ComplexF64, N, Nᵗ, 3)
    θ = [0.1, -0.2, 0.3, 0.02, 0.01, -0.03, (0.05 .* Random.randn(rng, 9))...]
    Bₘᵢₙ = 4 * Threads.nthreads()
    Nᵖ(L) = (L + 1) * (2L + 1)  # the pixels of `sphere_quadrature(L)`
    L = first(L for L ∈ ℓ:100 if Nᵖ(L) ≥ 3Bₘᵢₙ && Nᵖ(L) % Bₘᵢₙ ≠ 0)
    tiny = 1
    for kw ∈ ((;), (; ℓ_fit=3, normalize=true))
        plan = Scri.ObjectivePlan(target, t, dc, t′; grid=Scri.sphere_quadrature(L), kw...)
        @test length(plan.R′ₚ) == Nᵖ(L)
        # Any budget too small for the smallest blocks gives the same blocks.
        blocks = Scri.synthesis_blocks(data, plan.R′ₚ, plan.spins; budget=tiny)
        @test length(blocks) ≥ 3
        @test length(last(blocks)) < length(first(blocks))
        @test length(Scri.synthesis_blocks(data, plan.R′ₚ, plan.spins)) == 1

        Φ, ∇ = Scri.gradient_core(data, t, θ, dc, plan, nothing, Val(8))
        Φᵇ, ∇ᵇ = Scri.gradient_core(data, t, θ, dc, plan, nothing, Val(8), tiny)
        @test Φᵇ ≈ Φ rtol = 1e-13
        @test ∇ᵇ ≈ ∇ rtol = 1e-13
        J = Scri.transform_objective(data, t, Scri.bms_from_parameters(θ)..., dc, plan)
        @test Φᵇ ≈ J rtol = 1e-13

        # The objective's projection, in the tiny blocks, is the one in a single block.
        if plan.ℓ_fit ≥ 0
            cache = Scri.CubicSplineCache(t)
            geometry = Scri.pixel_geometry(
                t, t′, Scri.bms_from_parameters(θ)..., plan.R′ₚ, dc
            )
            args = (data, t, t′, cache, geometry, dc, plan.spins, N, plan.Ȳ)
            μ̂ = Scri.transformed_projection(args..., Float64)
            workspace = Scri.block_workspace(
                blocks, data, t′, geometry.Rₚ, plan.spins, Float64
            )
            μᵇ = Array{ComplexF64}(undef, length(t′), length(first(blocks)), 3)
            μ̂ᵇ = Scri.transformed_projection!(similar.(μ̂), μᵇ, blocks, args..., workspace)
            @test μ̂ᵇ ≈ μ̂ rtol = 1e-13
        end
    end
end

@testitem "objective_and_gradient!: Double64" tags = [:validation] begin
    using DoubleFloats: Double64
    import ForwardDiff
    import Random

    # In Double64 the products are generic rather than BLAS, also with ℓ_fit, where the
    # projected residual is evaluated back at the pixels by a product with its adjoint.
    rng = Random.Xoshiro(6565)
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 16
    t = collect(LinRange(Double64(-10), Double64(10), Nᵗ))
    t′ = Scri.compute_t′(t, Double64(0.1), Double64(1))
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Complex{Double64}.(Random.randn(rng, ComplexF64, N, Nᵗ, 2))
    target = Complex{Double64}.(Random.randn(rng, ComplexF64, N, Nᵗ, 2))
    θ = Double64.([0.1, -0.2, 0.3, 0.02, 0.01, -0.03, 0.05, -0.04, 0.03, 0.02])
    for kw ∈ ((;), (; ℓ_fit=2))
        plan = Scri.ObjectivePlan(target, t, dc, t′; kw...)
        J(θ) = Scri.transform_objective(data, t, Scri.bms_from_parameters(θ)..., dc, plan)
        ∇ = similar(θ)
        @test Scri.objective_and_gradient!(∇, θ, data, t, dc, plan) ≈ J(θ) rtol = 1e-28
        @test ∇ ≈ ForwardDiff.gradient(J, θ) rtol = 1e-25
    end
end
