# Tests for the quadrature rules (src/quadrature.jl)

@testitem "gauss_legendre: exact for polynomials of degree 2n-1, at any precision" tags = [
    :unit, :fast
] begin
    using DoubleFloats: Double64

    for T ∈ (Float64, Double64, BigFloat), n ∈ (1, 2, 5, 17, 33)
        x, w = Scri.gauss_legendre(n, T)
        @test eltype(x) === T && eltype(w) === T
        @test issorted(x) && all(-1 .< x .< 1) && all(w .> 0)
        @test x ≈ -reverse(x) atol = 4eps(T)
        for k ∈ 0:(2n - 1)
            exact = isodd(k) ? zero(T) : T(2) / (k + 1)
            @test sum(w .* x .^ k) ≈ exact atol = 50n * eps(T)
        end
    end
    # The two-point rule in closed form.
    x, w = Scri.gauss_legendre(2, BigFloat)
    @test x ≈ [-1, 1] / sqrt(BigFloat(3)) atol = 10eps(BigFloat)
    @test w ≈ [1, 1] atol = 10eps(BigFloat)
    @test_throws ArgumentError Scri.gauss_legendre(0)
end

@testitem "sphere_quadrature: Parseval for band-limited fields of every spin" tags = [
    :unit, :fast, :validation
] begin
    using SphericalFunctions: sYlm_matrix
    using DoubleFloats: Double64
    import Random

    # The weighted sum over the grid of f times the conjugate of g is the sum over modes of
    # the same, for any two fields of one spin weight band-limited to ℓₘₐₓ.
    rng = Random.Xoshiro(9191)
    for (T, ℓₘₐₓs) ∈ ((Float64, (0, 2, 8, 16)), (Double64, (2, 6)))
        for ℓₘₐₓ ∈ ℓₘₐₓs
            R, w = Scri.sphere_quadrature(ℓₘₐₓ, T)
            @test length(R) == length(w) == (ℓₘₐₓ + 1) * (2ℓₘₐₓ + 1)
            @test sum(w) ≈ 4T(π) rtol = 20eps(T)
            for s ∈ (-2):2
                abs(s) > ℓₘₐₓ && continue
                N = (ℓₘₐₓ + 1)^2 - s^2
                f = Complex{T}.(Random.randn(rng, ComplexF64, N))
                g = Complex{T}.(Random.randn(rng, ComplexF64, N))
                Y = sYlm_matrix(R, ℓₘₐₓ, s)
                @test sum(w .* abs2.(Y * f)) ≈ sum(abs2, f) rtol = 100eps(T)
                @test sum(w .* (Y * f) .* conj.(Y * g)) ≈ sum(f .* conj.(g)) atol =
                    100eps(T) * sqrt(sum(abs2, f) * sum(abs2, g))
            end
        end
    end
    @test_throws ArgumentError Scri.sphere_quadrature(-1)
end
