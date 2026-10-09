# Tests for ObjectivePlan (src/objective_plan.jl)

@testitem "ObjectivePlan: default grid, weights, and the target's norm" tags = [
    :unit, :fast, :validation
] begin
    import Random

    rng = Random.Xoshiro(9292)
    ℓ = 4
    N = (ℓ + 1)^2
    Nᵗ = 12
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    dc = Scri.DataComponents(:h, :ψ₄)
    target = Random.randn(rng, ComplexF64, N, Nᵗ, 2)

    # On the knots themselves, interpolation is exact, so the plan's norm is the
    # trapezoid sum over time of the sum over modes with ℓ ≥ |s| (Parseval, through the
    # exact sphere quadrature).
    plan = Scri.ObjectivePlan(target, t, dc, t)
    R, w = Scri.sphere_quadrature(ℓ)
    @test plan.R′ₚ == R && plan.w == w
    @test plan.target ≈ Scri.pixel_waveform(target, t, dc, t, R) rtol = 1e-14
    v = Scri.trapezoid_weights(t)
    @test plan.v == v && sum(v) ≈ t[end] - t[begin]
    @test plan.c == [1.0, 1.0]
    expected = sum(v[j] * sum(abs2, target[5:end, j, k]) for j ∈ 1:Nᵗ, k ∈ 1:2)
    @test plan.norm² ≈ expected rtol = 1e-12
    @test !plan.normalize
    @test Scri.trapezoid_weights([3.0]) == [1.0]

    # Weights multiply as documented.
    pw = rand(rng, length(R))
    p2 = Scri.ObjectivePlan(
        target, t, dc, t; weights=pw, time_weights=ones(Nᵗ), component_weights=[2.0, 0.5]
    )
    @test p2.w ≈ w .* pw && p2.v == ones(Nᵗ) && p2.c == [2.0, 0.5]

    # A grid may be given, here the golden-ratio pixels with unit weights.
    golden = Scri.golden_ratio_spiral_rotors(0, ℓ, Float64)
    p3 = Scri.ObjectivePlan(target, t, dc, t[2:5]; grid=(golden, ones(N)))
    @test length(p3.R′ₚ) == N && p3.target ≈ Scri.pixel_waveform(target, t, dc, t[2:5])

    # Bad weights and shapes are rejected.
    @test_throws ArgumentError Scri.ObjectivePlan(target, t, dc, t; weights=ones(3))
    @test_throws ArgumentError Scri.ObjectivePlan(target, t, dc, t; weights=(-pw))
    @test_throws ArgumentError Scri.ObjectivePlan(target, t, dc, t; time_weights=ones(2))
    @test_throws ArgumentError Scri.ObjectivePlan(target, t, dc, t; component_weights=[1.0])
    @test_throws ArgumentError Scri.ObjectivePlan(target, t, dc, t; grid=(golden, ones(3)))
    @test_throws ArgumentError Scri.ObjectivePlan(0 * target, t, dc, t; normalize=true)
    @test_throws ArgumentError Scri.ObjectivePlan(target, t, dc, t .+ 100)
end

@testitem "transform_objective: the default plan's minimum is the true L² minimum" tags = [
    :validation
] begin
    using Quaternionic: QuatVec, Rotor
    import ForwardDiff

    # The target is the data shifted in time by δ₂ in its ℓ = 2 modes and by δ₃ in its
    # ℓ = 3 modes, so that no single time translation matches it, and the best one depends
    # on how the norm weighs the two.  For a pure time translation the transformed field
    # stays band-limited, so with the exact sphere quadrature the objective must equal the
    # same weighted difference computed over modes, at every translation, and so must
    # have the same minimum.  The old objective, an unweighted sum over the pixels of
    # `transform!`, weighs the two degrees by amounts that depend on the field, and misses.
    ℓ = 3
    N = (ℓ + 1)^2
    Nᵗ = 41
    t = collect(LinRange(-20.0, 20.0, Nᵗ))
    dc = Scri.DataComponents(:h, :ψ₄)
    f(i, k, τ) = exp(-(τ / 12)^2) * cis(0.4τ + i + k) * (1 + i / 7 + k / 3)
    δ₂, δ₃ = 0.3, -0.25
    data = zeros(ComplexF64, N, Nᵗ, 2)
    target = zeros(ComplexF64, N, Nᵗ, 2)
    for k ∈ 1:2, i ∈ 5:N, j ∈ 1:Nᵗ
        data[i, j, k] = f(i, k, t[j])
        target[i, j, k] = f(i, k, t[j] + (i ≤ 9 ? δ₂ : δ₃))
    end
    t′ = t[6:(end - 5)]  # knots, so the target needs no interpolation
    plan = Scri.ObjectivePlan(target, t, dc, t′)
    golden = Scri.golden_ratio_spiral_rotors(0, ℓ, Float64)
    old = Scri.ObjectivePlan(target, t, dc, t′; grid=(golden, fill(4π / N, N)))
    # A time translation by δ is the supertranslation α = δ, whose only mode is ℓ = 0.
    J(p, δ) = Scri.transform_objective(
        data, t, QuatVec(0.0, 0.0, 0.0), one(Rotor{Float64}), [complex(2√π * δ)], dc, p
    )

    # The same quantity over modes, with each mode's time series interpolated by a natural
    # cubic spline from a dense solve, independent of the package's kernel.
    h = diff(t)
    A = zeros(Nᵗ, Nᵗ)
    A[1, 1] = A[Nᵗ, Nᵗ] = 1
    for j ∈ 2:(Nᵗ - 1)
        A[j, j - 1], A[j, j], A[j, j + 1] = h[j - 1], 2(h[j - 1] + h[j]), h[j]
    end
    M = Dict{Tuple{Int,Int},Vector{ComplexF64}}()
    for i ∈ 5:N, k ∈ 1:2
        y = data[i, :, k]
        b = zeros(ComplexF64, Nᵗ)
        for j ∈ 2:(Nᵗ - 1)
            b[j] = 6 * ((y[j + 1] - y[j]) / h[j] - (y[j] - y[j - 1]) / h[j - 1])
        end
        M[(i, k)] = A \ b
    end
    function spline(i, k, τ)
        j = clamp(searchsortedlast(t, τ), 1, Nᵗ - 1)
        a, c = (t[j + 1] - τ) / h[j], (τ - t[j]) / h[j]
        y, m = view(data, :, :, k), M[(i, k)]
        return a * y[i, j] +
               c * y[i, j + 1] +
               ((a^3 - a) * m[j] + (c^3 - c) * m[j + 1]) * h[j]^2 / 6
    end
    v = Scri.trapezoid_weights(t′)
    Jmodes(δ) = sum(
        v[j] * abs2(spline(i, k, t′[j] + δ) - target[i, j + 5, k]) for
        j ∈ eachindex(t′), k ∈ 1:2, i ∈ 5:N
    )

    for δ ∈ (0.0, 0.1, 0.27)
        @test J(plan, δ) ≈ Jmodes(δ) rtol = 1e-13
    end

    # Minimize each over δ by the secant method on its ForwardDiff derivative.
    function minimizer(J)
        g(δ) = ForwardDiff.derivative(J, δ)
        a, b = -0.1, 0.1
        ga, gb = g(a), g(b)
        for _ ∈ 1:40
            a, ga, b = b, gb, b - gb * (b - a) / (gb - ga)
            gb = g(b)
            abs(b - a) < 1e-13 && break
        end
        return b
    end
    δᵐ = minimizer(Jmodes)
    @test minimizer(δ -> J(plan, δ)) ≈ δᵐ atol = 1e-10
    # Measured: the old objective's minimum is off by 0.012, about 2% of δ₂ − δ₃.
    @test abs(minimizer(δ -> J(old, δ)) - δᵐ) > 5e-3
end

@testitem "ObjectivePlan: ℓ_fit compares only the modes up to ℓ_fit" tags = [
    :unit, :validation
] begin
    using Quaternionic: QuatVec, Rotor
    import ForwardDiff
    import Random

    rng = Random.Xoshiro(9393)
    ℓ = 5
    N = (ℓ + 1)^2
    Nᵗ = 16
    t = collect(LinRange(-10.0, 10.0, Nᵗ))
    dc = Scri.DataComponents(:h, :ψ₄)
    data = Random.randn(rng, ComplexF64, N, Nᵗ, 2)
    modes = Random.randn(rng, ComplexF64, N, Nᵗ, 2)
    data[1:4, :, :] .= 0
    modes[1:4, :, :] .= 0
    t′ = t[4:(end - 3)]  # knots, so nothing is interpolated at the identity
    v = Scri.trapezoid_weights(t′)
    identity(p) = Scri.transform_objective(
        data,
        t,
        QuatVec(0.0, 0.0, 0.0),
        one(Rotor{Float64}),
        zeros(ComplexF64, 1),
        dc,
        p,
    )

    # At the identity, the projected difference is the difference of the modes themselves.
    for L ∈ (2, 3, 5)
        p = Scri.ObjectivePlan(modes, t, dc, t′; ℓ_fit=L)
        n = 5:((L + 1) ^ 2)
        expected = sum(
            v[j] * sum(abs2, data[n, j + 3, k] .- modes[n, j + 3, k]) for
            j ∈ eachindex(t′), k ∈ 1:2
        )
        @test identity(p) ≈ expected rtol = 1e-12
    end
    @test identity(Scri.ObjectivePlan(modes, t, dc, t′; ℓ_fit=ℓ)) ≈
        identity(Scri.ObjectivePlan(modes, t, dc, t′)) rtol = 1e-12

    # Normalization uses the norm of the target's projection.
    pn = Scri.ObjectivePlan(modes, t, dc, t′; ℓ_fit=3, normalize=true)
    @test pn.norm² ≈
        sum(v[j] * sum(abs2, modes[5:16, j + 3, k]) for j ∈ eachindex(t′), k ∈ 1:2) rtol =
        1e-12

    # The derivative through the projection, here with respect to a boost, matches central
    # differences.
    p₃ = Scri.ObjectivePlan(modes, t, dc, t′; ℓ_fit=3)
    f(β) = Scri.transform_objective(
        data,
        t,
        QuatVec(zero(β), zero(β), β),
        one(Rotor{Float64}),
        zeros(ComplexF64, 1),
        dc,
        p₃,
    )
    h = 1e-6
    d = ForwardDiff.derivative(f, 1e-3)
    fd = (f(1e-3 + h) - f(1e-3 - h)) / 2h
    @test isfinite(d) && abs(d - fd) < 1e-6 * max(1.0, abs(fd))

    # ℓ_fit must cover the spin weights, and stay within the default grid's band limit.
    @test_throws ArgumentError Scri.ObjectivePlan(modes, t, dc, t′; ℓ_fit=1)
    @test_throws ArgumentError Scri.ObjectivePlan(modes, t, dc, t′; ℓ_fit=ℓ + 1)
end
