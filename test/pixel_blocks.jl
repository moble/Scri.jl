# Tests for the blocked synthesis of pixel time series (src/pixel_blocks.jl)

@testitem "pixel_blocks: consecutive blocks covering every pixel" tags = [:unit, :fast] begin
    budget = 2^25
    for Nᵖ ∈ (1, 7, 100, 1001), bytes ∈ (1, 10^4, 10^6, 10^9)
        blocks = Scri.pixel_blocks(Nᵖ, bytes, budget)
        @test reduce(vcat, collect.(blocks)) == 1:Nᵖ
        B = length(first(blocks))
        @test all(length(b) == B for b ∈ blocks[1:(end - 1)])
        @test B == Nᵖ || B ≥ 4 * Threads.nthreads()
        @test B == Nᵖ || B * bytes ≤ budget || B == 4 * Threads.nthreads()
    end
end

@testitem "synthesize_block!: matches pixel_dot, with and without BLAS" tags = [
    :unit, :fast
] begin
    using DoubleFloats: Double64
    using Quaternionic: Rotor, components, randn
    using SphericalFunctions: sYlmCalculator
    import ForwardDiff
    import Random

    rng = Random.Xoshiro(2027)
    ℓₘₐₓ, Nᵗ, Nᵖ = 4, 7, 11
    Nᵐ = (ℓₘₐₓ + 1)^2
    # Two components share a spin weight, and the spin weight −1 between them is unused.
    spins = [-2, 0, -2]
    spin_range = minimum(spins):maximum(spins)
    # Uneven blocks, the last one shorter than the space allocated for it.
    blocks = [1:4, 5:8, 9:11]
    B = 4

    # The real and imaginary parts of a complex number, with any partial derivatives.
    parts(x::Real) = (x,)
    parts(x::ForwardDiff.Dual) = (ForwardDiff.value(x), ForwardDiff.partials(x)...)
    parts(z::Complex) = (parts(real(z))..., parts(imag(z))...)
    flat(A) = collect(Iterators.flatten(parts.(A)))

    D = ForwardDiff.Dual{Nothing,Float64,2}
    dual(x) = D(x, ForwardDiff.Partials((Random.randn(rng), Random.randn(rng))))
    rotors = randn(rng, Rotor{Float64}, Nᵖ)
    cases = (
        # element type of the rotors, of the data, the band limit N₀, and exact agreement
        (Float64, Float64, Nᵐ, false),  # BLAS
        (Float64, Float64, 16, false),  # BLAS, with a declared band limit below ℓₘₐₓ
        (Double64, Double64, Nᵐ, true),  # pixel_dot
        (Float64, Float32, Nᵐ, true),  # pixel_dot, since the types differ
        (D, Float64, 16, false),  # dual harmonics, primal data
    )
    for (Tᴿ, Tᵈ, N₀, exact) ∈ cases
        Rₚ = if Tᴿ === D
            [Rotor{D}(dual.(components(R))...) for R ∈ rotors]
        else
            [Rotor{Tᴿ}(Tᴿ.(components(R))...) for R ∈ rotors]
        end
        data = Random.randn(rng, Complex{Tᵈ}, Nᵐ, Nᵗ, length(spins))
        calc = sYlmCalculator(one(eltype(Rₚ)), ℓₘₐₓ, spin_range)
        Yᵢ = Matrix{Complex{Tᴿ}}(undef, Nᵐ, length(spin_range))
        Tˣ = promote_type(Complex{Tᴿ}, Complex{Tᵈ})
        expected = Array{Tˣ}(undef, Nᵗ, Nᵖ, length(spins))
        for i ∈ 1:Nᵖ
            Scri.sYlm_rows!(Yᵢ, calc, Rₚ[i])
            for k ∈ eachindex(spins), h ∈ 1:Nᵗ
                s = spins[k]
                Yˢ = view(Yᵢ, :, s - first(spin_range) + 1)
                expected[h, i, k] = Scri.pixel_dot(Yˢ, view(data, :, h, k), s^2 + 1, N₀)
            end
        end

        calcs = [sYlmCalculator(one(eltype(Rₚ)), ℓₘₐₓ, spin_range) for _ ∈ 1:3]
        Y = Array{Complex{Tᴿ}}(undef, Nᵐ, B, length(spin_range))
        X = Array{Tˣ}(undef, Nᵗ, B, length(spins))
        for block ∈ blocks
            nb = length(block)
            Yb, Xb = view(Y, :, 1:nb, :), view(X, :, 1:nb, :)
            Scri.harmonics_block!(Yb, calcs, Rₚ, block)
            Scri.synthesize_block!(Xb, data, Yb, spins, N₀)
            if exact
                @test Xb == expected[:, block, :]
            else
                @test flat(Xb) ≈ flat(expected[:, block, :]) rtol = 1e-13
            end
        end
    end
end

@testitem "stack_partials: one BLAS product for primal times dual matrices" tags = [
    :unit, :fast
] begin
    import ForwardDiff
    import LinearAlgebra
    import Random

    rng = Random.Xoshiro(4141)
    D = ForwardDiff.Dual{Nothing,Float64,3}
    dual() = D(Random.randn(rng), ForwardDiff.Partials(Tuple(Random.randn(rng, 3))))
    B = [complex(dual(), dual()) for _ ∈ 1:5, _ ∈ 1:4]
    A = Random.randn(rng, ComplexF64, 6, 5)
    C = [complex(dual(), dual()) for _ ∈ 1:6, _ ∈ 1:4]

    parts(x::ForwardDiff.Dual) = (ForwardDiff.value(x), ForwardDiff.partials(x)...)
    parts(z::Complex) = (parts(real(z))..., parts(imag(z))...)
    flat(M) = collect(Iterators.flatten(parts.(M)))

    # Stacking and unstacking are inverses.
    @test Scri.unstack_partials!(similar(B), Scri.stack_partials(B), false) == B

    # The stacked BLAS product agrees with the generic product of dual numbers.
    expected = LinearAlgebra.mul!(copy(C), A, B, true, true)
    @test flat(Scri.add_product!(copy(C), A, B)) ≈ flat(expected) rtol = 1e-14
    @test flat(Scri.add_product!(copy(C), view(A, :, 2:4), view(B, 2:4, :))) ≈
        flat(LinearAlgebra.mul!(copy(C), view(A, :, 2:4), view(B, 2:4, :), true, true)) rtol =
        1e-14
end
