# Tests for CubicSplineCache (src/cubic_spline.jl)

@testmodule CubicSplineSetup begin
    using Scri
    using DoubleFloats: Double64
    const FloatTypes = (Float32, Float64, Double64, BigFloat)

    # Complete natural cubic spline interpolation using the cache.
    # Used only in tests; mirrors the forward/backward sweep in transform.jl.
    function cubic_spline_interp(
        cache::Scri.CubicSplineCache{T}, t, d::AbstractVector, t_query
    ) where {T}
        N = length(d)
        d̈ = zeros(T, N)
        # Forward sweep: z[j] stored in d̈[j]
        d̈[2] =
            T(6) *
            (cache.h⁻¹[2] * (d[3] - d[2]) - cache.h⁻¹[1] * (d[2] - d[1])) *
            cache.u⁻¹[1]
        for j ∈ 3:(N - 1)
            r =
                T(6) *
                (cache.h⁻¹[j] * (d[j + 1] - d[j]) - cache.h⁻¹[j - 1] * (d[j] - d[j - 1]))
            d̈[j] = (r - cache.h[j - 1] * d̈[j - 1]) * cache.u⁻¹[j - 1]
        end
        # Backward sweep
        for j ∈ (N - 2):-1:2
            d̈[j] -= cache.l[j - 1] * d̈[j + 1]
        end
        # Locate interval and evaluate
        j = clamp(searchsortedlast(t, t_query), 1, N - 1)
        τ = t_query - t[j]
        return Scri.spline_eval(d[j], d[j + 1], d̈[j], d̈[j + 1], cache.h[j], cache.h⁻¹[j], τ)
    end
end

@testitem "CubicSplineCache: rejects too few knots" tags = [:unit, :fast] begin
    @test_throws ArgumentError Scri.CubicSplineCache(collect(Float64, 1:3))
    @test_throws ArgumentError Scri.CubicSplineCache(collect(Float64, 1:2))
    @test_throws ArgumentError Scri.CubicSplineCache(collect(Float64, 1:1))
end

@testitem "CubicSplineCache: rejects non-increasing knots" tags = [:unit, :fast] begin
    @test_throws ArgumentError Scri.CubicSplineCache([0.0, 1.0, 1.0, 2.0])
    @test_throws ArgumentError Scri.CubicSplineCache([0.0, 2.0, 1.0, 3.0])
end

@testitem "CubicSplineCache: array dimensions" tags = [:unit, :fast] begin
    for N ∈ [4, 5, 10, 20]
        t = collect(Float64, range(0.0, 1.0; length=N))
        c = Scri.CubicSplineCache(t)
        @test length(c.h) == N - 1
        @test length(c.h⁻¹) == N - 1
        @test length(c.l) == N - 3
        @test length(c.u⁻¹) == N - 2
        @test all(c.h .≈ diff(t))
        @test all(c.h⁻¹ .≈ inv.(diff(t)))
    end
end

@testitem "CubicSplineCache: h and h⁻¹ are consistent" tags = [:unit, :fast] setup = [
    CubicSplineSetup
] begin
    import Random
    using .CubicSplineSetup: FloatTypes

    rng = Random.Xoshiro(42)
    for T ∈ FloatTypes
        dt = T.(abs.(randn(rng, 19))) .+ T(0.01)
        t = T.(cumsum(dt))
        c = Scri.CubicSplineCache(t)
        @test all(c.h .* c.h⁻¹ .≈ ones(T, length(c.h)))
    end
end

@testitem "CubicSplineCache: interpolates exactly at knots" tags = [:unit, :fast] setup = [
    CubicSplineSetup
] begin
    import Random
    using .CubicSplineSetup: FloatTypes, cubic_spline_interp

    rng = Random.Xoshiro(7)
    for T ∈ FloatTypes
        N = 12
        dt = T.(abs.(randn(rng, N - 1))) .+ T(0.1)
        t = [zero(T); T.(cumsum(dt))]   # N elements: prepend t=0
        d = T.(randn(rng, N))
        c = Scri.CubicSplineCache(t)
        for j ∈ 1:N
            val = cubic_spline_interp(c, t, d, t[j])
            @test val ≈ d[j] atol = 100 * eps(T) * maximum(abs, d)
        end
    end
end

@testitem "CubicSplineCache: exact for linear polynomials" tags = [
    :unit, :fast, :validation
] setup = [CubicSplineSetup] begin
    import Random
    using .CubicSplineSetup: FloatTypes, cubic_spline_interp

    # A natural cubic spline must reproduce any linear polynomial exactly
    # (up to floating-point rounding): linear data has zero second divided
    # differences, so the Thomas algorithm produces d̈ = 0 everywhere and
    # spline_eval reduces to linear interpolation.
    rng = Random.Xoshiro(13)
    for T ∈ FloatTypes
        N = 15
        dt = T.(abs.(randn(rng, N - 1))) .+ T(0.1)
        t = T.(cumsum(dt))
        # random linear: p(x) = a + b·x
        a, b = T.(randn(rng, 2))
        poly(x) = a + b*x
        d = poly.(t)
        cache = Scri.CubicSplineCache(t)
        # query at 30 interior points
        for _ ∈ 1:30
            xq = t[1] + rand(rng) * (t[end] - t[1])
            s = cubic_spline_interp(cache, t, d, T(xq))
            @test s ≈ poly(T(xq)) atol = 100 * eps(T) * (abs(a) + abs(b)*t[end] + 1)
        end
    end
end

@testitem "CubicSplineCache: O(h⁴) convergence on smooth function" tags = [
    :unit, :validation
] setup = [CubicSplineSetup] begin
    using .CubicSplineSetup: cubic_spline_interp

    # Verify that the max error decreases as h⁴ under uniform grid refinement.
    # Test on f(x) = sin(x) over [0, 2π].
    f = sin
    x_test = range(0.1, 2π - 0.1; length=40)

    function max_err(N)
        t = collect(Float64, range(0.0, 2π; length=N))
        d = f.(t)
        c = Scri.CubicSplineCache(t)
        maximum(x_test) do xq
            return abs(cubic_spline_interp(c, t, d, xq) - f(xq))
        end
    end

    err16 = max_err(16)
    err32 = max_err(32)
    err64 = max_err(64)
    # Doubling N ≈ halving h; error ratio should be ≈ 2⁴ = 16.
    ratio1 = err16 / err32
    ratio2 = err32 / err64
    @test ratio1 > 10   # O(h⁴): ratio ≈ 16; allow slack for small N
    @test ratio2 > 10
end

@testitem "CubicSplineCache: non-uniform grid" tags = [:unit, :fast, :validation] setup = [
    CubicSplineSetup
] begin
    import Random
    using .CubicSplineSetup: FloatTypes, cubic_spline_interp

    # On a non-uniform grid, the spline must still reproduce linear polynomials exactly.
    # (Natural BC d̈[endpoints]=0 is satisfied by linear data, so the spline IS the line.)
    rng = Random.Xoshiro(99)
    for T ∈ FloatTypes
        N = 10
        dt = T.(sort(rand(rng, N - 1)) .* T(4) .+ T(0.05))
        t = T.(cumsum(dt))
        a, b = T.(randn(rng, 2))
        poly(x) = a + b*x
        d = poly.(t)
        cache = Scri.CubicSplineCache(t)
        for _ ∈ 1:20
            xq = t[1] + rand(rng) * (t[end] - t[1])
            s = cubic_spline_interp(cache, t, d, T(xq))
            @test s ≈ poly(T(xq)) atol = 1000 * eps(T) * (abs(a) + abs(b)*t[end] + 1)
        end
    end
end

@testitem "CubicSplineCache: natural BC gives zero curvature at endpoints" tags = [
    :unit, :fast, :validation
] setup = [CubicSplineSetup] begin
    import Random
    using .CubicSplineSetup: FloatTypes

    # The natural spline has d̈[1]=d̈[N]=0.  We verify this by reconstructing d̈
    # from the cache and checking the endpoint values, and that the interior
    # tridiagonal system A·d̈ = r is satisfied to machine precision.
    rng = Random.Xoshiro(55)
    for T ∈ (Float64,)
        N = 12
        dt = T.(abs.(randn(rng, N - 1))) .+ T(0.1)
        t = [zero(T); T.(cumsum(dt))]   # N elements: prepend t=0
        d = T.(randn(rng, N))
        c = Scri.CubicSplineCache(t)
        # Forward sweep
        d̈ = zeros(T, N)
        d̈[2] = T(6) * (c.h⁻¹[2]*(d[3]-d[2]) - c.h⁻¹[1]*(d[2]-d[1])) * c.u⁻¹[1]
        for j ∈ 3:(N - 1)
            r = T(6) * (c.h⁻¹[j]*(d[j + 1]-d[j]) - c.h⁻¹[j - 1]*(d[j]-d[j - 1]))
            d̈[j] = (r - c.h[j - 1] * d̈[j - 1]) * c.u⁻¹[j - 1]
        end
        # Backward sweep
        for j ∈ (N - 2):-1:2
            d̈[j] -= c.l[j - 1] * d̈[j + 1]
        end
        @test d̈[1] == 0
        @test d̈[end] == 0
        # The interior system A·d̈[2..N-1] = r should be satisfied to machine precision.
        for j ∈ 2:(N - 1)
            r = T(6) * (c.h⁻¹[j]*(d[j + 1]-d[j]) - c.h⁻¹[j - 1]*(d[j]-d[j - 1]))
            lhs =
                (j > 2 ? c.h[j - 1]*d̈[j - 1] : 0) +
                2*(c.h[j - 1]+c.h[j])*d̈[j] +
                (j < N-1 ? c.h[j]*d̈[j + 1] : 0)
            @test lhs ≈ r atol = 1e4 * eps(T) * (abs(r) + 1)
        end
    end
end

@testitem "CubicSplineCache: spline_eval endpoints match knot values" tags = [:unit, :fast] begin
    import Random

    rng = Random.Xoshiro(88)
    N = 8
    t = sort(rand(rng, N)) .* 10.0
    d = randn(rng, N)
    c = Scri.CubicSplineCache(t)
    for j ∈ 1:(N - 1)
        # τ=0 should give d[j], τ=h[j] should give d[j+1]
        @test Scri.spline_eval(d[j], d[j + 1], 0.0, 0.0, c.h[j], c.h⁻¹[j], 0.0) == d[j]
        @test Scri.spline_eval(d[j], d[j + 1], 0.0, 0.0, c.h[j], c.h⁻¹[j], c.h[j]) ≈
            d[j + 1]
    end
end

@testitem "transform_pixel!: matches a dense natural-spline solve" tags = [
    :unit, :fast, :validation
] begin
    import Random

    # The other spline tests check a copy of the Thomas sweep; this one checks the kernel
    # that `transform!` and `transform_objective` actually call.  The oracle solves the full
    # natural-spline system for the second derivatives densely and evaluates the textbook
    # piecewise cubic, then applies the same mixing at the same transformed times.
    function dense_natural_spline(t, y, x)
        N = length(t)
        h = diff(t)
        A = zeros(N, N)
        b = zeros(eltype(y), N)
        A[1, 1] = A[N, N] = 1
        for j ∈ 2:(N - 1)
            A[j, j - 1], A[j, j], A[j, j + 1] = h[j - 1], 2(h[j - 1] + h[j]), h[j]
            b[j] = 6 * ((y[j + 1] - y[j]) / h[j] - (y[j] - y[j - 1]) / h[j - 1])
        end
        M = A \ b
        j = clamp(searchsortedlast(t, x), 1, N - 1)
        a = (t[j + 1] - x) / h[j]
        c = (x - t[j]) / h[j]
        return a * y[j] +
               c * y[j + 1] +
               ((a^3 - a) * M[j] + (c^3 - c) * M[j + 1]) * h[j]^2 / 6
    end

    rng = Random.Xoshiro(8484)
    dc = Scri.DataComponents(:h, :ψ₄, :ψ₃)
    κ⁻¹, αₚ, ð₀, ð₁, ð²α = 1.1, 0.05, 0.02 - 0.01im, -0.03 + 0.02im, 0.04 + 0.03im
    # The smallest allowed input (four knots) with a single output time, then a longer
    # non-uniform grid with outputs in several intervals.
    for (t, t′) ∈ (
        ([-1.0, -0.3, 0.4, 2.0], [0.2]),
        ([-2.0, -1.6, -0.7, 0.0, 0.3, 1.5, 2.2], [-1.6, -1.0, -0.2, 0.6, 1.4]),
    )
        Nᵗ, Nᵗ′ = length(t), length(t′)
        dᵢ = Random.randn(rng, ComplexF64, 3, Nᵗ)
        d′ᵢ = Matrix{ComplexF64}(undef, 3, Nᵗ′)
        Scri.transform_pixel!(
            d′ᵢ,
            copy(dᵢ),
            similar(dᵢ),
            Scri.CubicSplineCache(t),
            t,
            t′,
            κ⁻¹,
            αₚ,
            ð₀,
            ð₁,
            ð²α,
            dc,
        )
        for j′ ∈ 1:Nᵗ′
            τ = t′[j′] * κ⁻¹ + αₚ
            expected = [dense_natural_spline(t, dᵢ[k, :], τ) for k ∈ 1:3]
            Scri.mix_components!(expected, κ⁻¹, ð₀ + τ * ð₁, ð²α, dc)
            @test d′ᵢ[:, j′] ≈ expected rtol = 1e-13
        end
    end
end
