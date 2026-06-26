# Tests for the Conventions struct (src/conventions.jl)

@testitem "Conventions: defaults are SpEC + Newman–Penrose ð" tags = [:unit, :fast] begin
    import Scri: Conventions

    c = Conventions()
    @test c.c_s == 1
    @test c.c_R == 1
    @test c.c_Ψ == 1
    @test c.c_σ == 1
    @test c.c_φ == 1
    @test c.c_l == 1
    @test c.c_m == 0
    @test c.c_h == 1
    @test c.c_ð == 1
    @test c isa Conventions{Float64}
    # The default coincides with the :SpEC preset.
    @test c == Conventions(:SpEC)
end

@testitem "Conventions: named presets match the appendix table" tags = [:unit, :fast] begin
    import Scri: Conventions

    # (c_s, c_Ψ, c_σ, c_R, c_l, c_m, c_h) ↔ (s₀, s₁, s₂, s₃, λ, Θ, ζ); N/A → default 1.
    rt2 = sqrt(2.0)
    sp = Conventions(:SpEC)
    @test (sp.c_s, sp.c_Ψ, sp.c_σ, sp.c_R, sp.c_l, sp.c_m, sp.c_h) == (1, 1, 1, 1, 1, 0, 1)

    mb = Conventions(:MB)
    @test (mb.c_s, mb.c_Ψ, mb.c_σ, mb.c_R, mb.c_m) == (-1, 1, 1, 1, 0)
    @test mb.c_l ≈ -√2 && mb.c_h == 2

    np = Conventions(:NP)
    @test (np.c_s, np.c_Ψ, np.c_σ, np.c_R) == (-1, -1, 1, 1)
    @test np.c_l ≈ -√2 && np.c_m ≈ π

    adlk = Conventions(:ADLK)
    @test (adlk.c_s, adlk.c_Ψ, adlk.c_σ, adlk.c_R, adlk.c_m) == (1, 1, -1, 1, 0)
    @test adlk.c_l ≈ -√2

    br = Conventions(:BR)
    @test (br.c_s, br.c_Ψ, br.c_R, br.c_l, br.c_m) == (1, -1, 1, 1, 0)

    cc = Conventions(:C)
    @test (cc.c_s, cc.c_Ψ, cc.c_σ, cc.c_R) == (-1, -1, 1, 1)
    @test cc.c_l ≈ -√2

    # c_φ and c_ð are not in the GW table; presets leave them at the NP default.
    for name ∈ (:SpEC, :MB, :NP, :ADLK, :BR, :C)
        @test Conventions(name).c_φ == 1
        @test Conventions(name).c_ð == 1
    end
end

@testitem "Conventions: construction is validated" tags = [:unit, :fast] begin
    import Scri: Conventions

    @test_throws ArgumentError Conventions(; c_s=2)        # sign must be ±1
    @test_throws ArgumentError Conventions(; c_Ψ=0)
    @test_throws ArgumentError Conventions(; c_l=0)        # scales must be nonzero
    @test_throws ArgumentError Conventions(; c_h=0)
    @test_throws ArgumentError Conventions(; c_ð=0)
    @test_throws ArgumentError Conventions(:NotAConvention)
end

@testitem "Conventions: precision is carried by the type parameter" tags = [:unit, :fast] begin
    import Scri: Conventions
    using DoubleFloats: Double64

    np = Conventions(:NP; T=Double64)
    @test np isa Conventions{Double64}
    @test np.c_l ≈ -sqrt(Double64(2)) atol = 10eps(Double64)
    @test abs(np.c_l + sqrt(Double64(2))) < 1e-30   # full Double64 precision, not Float64
    @test np.c_m == Double64(π)
end

@testitem "Conventions: conversion factors and round-trips" tags = [:unit, :fast] begin
    import Scri: Conventions, weyl_factor, strain_factor, convert_weyl, convert_strain

    sp = Conventions(:SpEC)
    # SpEC is the reference, so all its factors are unity.
    for n ∈ 0:4
        @test weyl_factor(sp, n) ≈ 1
    end
    @test strain_factor(sp) ≈ 1

    # NP vs SpEC: ψₙ picks up c_s c_Ψ c_R (c_l e^{i c_m})^{2-n} = (√2)^{2-n}
    # (since c_s c_Ψ = 1, c_l e^{iπ} = +√2).
    np = Conventions(:NP)
    for n ∈ 0:4
        @test weyl_factor(np, n) ≈ (sqrt(2.0))^(2 - n)
    end

    # Round-trips: converting A→B→A is the identity, for any pair.
    import Random
    rng = Random.Xoshiro(8)
    for (A, B) ∈ ((sp, np), (Conventions(:MB), Conventions(:ADLK)), (np, Conventions(:C)))
        for n ∈ 0:4
            ψ = randn(rng, ComplexF64)
            @test convert_weyl(convert_weyl(ψ, n, A, B), n, B, A) ≈ ψ
        end
        h = randn(rng, ComplexF64)
        @test convert_strain(convert_strain(h, A, B), B, A) ≈ h
    end
end
