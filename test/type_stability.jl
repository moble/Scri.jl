# Tests that the main entry points compile to code with no dynamic dispatch in Scri itself,
# as reported by JET's optimization analysis.
#
# `target_modules=(Scri,)` keeps only the reports whose dynamic dispatch happens in Scri's
# own code.  Reports from other packages, and from Base, are left out on purpose: they are
# in error-message printing and task scheduling, which run once per call or once per block
# of pixels, not once per pixel.
#
# The keyword forms of `transform!`, whose `data_components` are given as symbols, are not
# tested: the type of the `DataComponents` they build is not known until the call runs, so
# they dispatch once to the positional form by design.

@testitem "JET: no dynamic dispatch in the main entry points" tags = [:type_stability] begin
    using JET
    using Quaternionic: QuatVec, Rotor, vec, components
    using DoubleFloats: Double64
    using Random: Xoshiro

    rng = Xoshiro(1234)
    function inputs(T, ℓₘₐₓ, Nᵗ, Nᵈ)
        Nᵐ = (ℓₘₐₓ + 1)^2
        t = collect(LinRange(T(-50), T(50), Nᵗ))
        v⃗ = T(1e-3) * randn(rng, QuatVec{T})
        R = randn(rng, Rotor{T})
        αᵢₙ = T(1e-3) * randn(rng, Complex{T}, Nᵐ)
        data = randn(rng, Complex{T}, Nᵐ, Nᵗ, Nᵈ)
        return data, t, v⃗, R, αᵢₙ
    end

    for dc ∈ (
        DataComponents(:h, :ψ₄),
        DataComponents(:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
        DataComponents(:λ, :ψ₀, :ψ₁; ℐ=-1),
        DataComponents(:ψ₄, :h; conventions=Conventions(:MB)),
    )
        data, t, v⃗, R, αᵢₙ = inputs(Float64, 4, 41, Scri.ncomponents(dc))
        plan = TransformPlan(data, dc)
        t′ = compute_t′(t, 1e-2, 1e-2)
        @test_opt target_modules = (Scri,) transform!(data, t, v⃗, R, αᵢₙ, dc)
        @test_opt target_modules = (Scri,) transform!(
            data, t, v⃗, R, αᵢₙ, dc; t′, plan, ℓₘₐₓ₀=2
        )
        @test_opt target_modules = (Scri,) transform!(
            data, t, collect(vec(v⃗)), collect(components(R)), αᵢₙ, dc
        )
        @test_opt target_modules = (Scri,) TransformPlan(data, dc)
        @test_opt target_modules = (Scri,) pixel_waveform(data, t, dc, t′)
        objective_plan = ObjectivePlan(data, t, dc, t′)
        θ = zeros(6 + 9)
        @test_opt target_modules = (Scri,) ObjectivePlan(data, t, dc, t′)
        @test_opt target_modules = (Scri,) transform_objective(
            data, t, v⃗, R, αᵢₙ, dc, objective_plan
        )
        # On Julia 1.12, inference gives up on part of the gradient's call graph and leaves
        # it unoptimized, which JET reports as "failed to optimize due to recursion."  The
        # loop in `supertranslation_on_grid` runs under OhMyThreads' serial fallback, which
        # is a `mapreduce`, and normalizing the dual-number rotor in `aberration` calls
        # `hypot`, which 1.12 implements with another `mapreduce`.  Inference takes the
        # second call for unbounded recursion and stops optimizing every frame between the
        # two.  Julia 1.13 does not do this, and the code runs correctly either way, so on
        # 1.12 we accept those reports but still reject any runtime dispatch.
        if VERSION ≥ v"1.13"
            @test_opt target_modules = (Scri,) objective_and_gradient!(
                similar(θ), θ, data, t, dc, objective_plan
            )
        else
            result = @report_opt target_modules = (Scri,) objective_and_gradient!(
                similar(θ), θ, data, t, dc, objective_plan
            )
            @test all(r -> r isa JET.OptimizationFailureReport, JET.get_reports(result))
        end
    end

    # Other precisions take the same code paths, with different element types.
    let dc = DataComponents(:h, :ψ₄)
        for T ∈ (Float32, Double64)
            data, t, v⃗, R, αᵢₙ = inputs(T, 4, 41, 2)
            @test_opt target_modules = (Scri,) transform!(data, t, v⃗, R, αᵢₙ, dc)
        end
    end
end

@testitem "JET: no dynamic dispatch in the work transform! runs in tasks" tags = [
    :type_stability
] begin
    # JET cannot follow a call into the body of a task, so the functions that `transform!`
    # runs with `OhMyThreads.@spawn` are not covered by the test above.  They are tested
    # here directly, with the arguments `transform!` would give them.
    using JET
    using Quaternionic: QuatVec, Rotor, Boost
    using Random: Xoshiro
    using Scri:
        transform_grid,
        supertranslation_on_grid,
        augmented_factorization,
        CubicSplineCache,
        compute_ðt′╱2κ,
        compute_t′

    rng = Xoshiro(1234)
    ℓₘₐₓ = 4
    t = collect(LinRange(-50.0, 50.0, 41))
    v⃗ = 1e-3 * randn(rng, QuatVec{Float64})
    R = randn(rng, Rotor{Float64})
    αᵢₙ = 1e-3 * randn(rng, ComplexF64, (ℓₘₐₓ + 1)^2)
    for ℐ ∈ (1, -1)
        R′ₚ = transform_grid(ℓₘₐₓ, Float64)
        Λℐ = Boost(ℐ * v⃗) * R
        Rₚ, αₚ, ðαₚ, _ = supertranslation_on_grid(R′ₚ, Λℐ, αᵢₙ, 1)
        @test_opt target_modules = (Scri,) supertranslation_on_grid(R′ₚ, Λℐ, αᵢₙ, 1)
        @test_opt target_modules = (Scri,) augmented_factorization(R′ₚ, ℓₘₐₓ, -2)
        @test_opt target_modules = (Scri,) CubicSplineCache(t)
        @test_opt target_modules = (Scri,) compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, ℐ)
        @test_opt target_modules = (Scri,) compute_t′(t, αₚ, Rₚ, v⃗, ℐ)
    end
end
