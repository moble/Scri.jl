# Tests for the Conventions struct (src/conventions.jl)

@testitem "Conventions: defaults are SpEC + Newman–Penrose ð" tags = [:unit, :fast] begin
    import Scri: Conventions, One

    c = Conventions()
    @test c.c_s == 1
    @test c.c_R == 1
    @test c.c_ψ == 1
    @test c.c_σ == 1
    @test c.c_φ == 1
    @test c.c_α == 1
    @test c.c_l == 1
    @test c.c_m == 1          # the spin-phase *factor* e^{i·0} = 1, not the angle
    @test c.c_h == 1
    @test c.c_ð == 1
    # Every default field is the `One` singleton → maximal compile-time elision.
    for f ∈ fieldnames(Conventions)
        @test getfield(c, f) === One()
    end
end

@testitem "Conventions: construction is validated" tags = [:unit, :fast] begin
    import Scri: Conventions

    @test_throws ArgumentError Conventions(; c_s=2)        # sign must be ±1
    @test_throws ArgumentError Conventions(; c_ψ=0)
    @test_throws ArgumentError Conventions(; c_α=2)        # c_α is a sign too
    @test_throws ArgumentError Conventions(; c_l=0)        # factors must be nonzero
    @test_throws ArgumentError Conventions(; c_m=0)        # c_m is a (nonzero) factor now
    @test_throws ArgumentError Conventions(; c_h=0)
    @test_throws ArgumentError Conventions(; c_ð=0)

    # ℓ is a real null vector, so c_l must be real (any nonzero real scale is fine).
    @test_throws ArgumentError Conventions(; c_l=1im)
    @test_throws ArgumentError Conventions(; c_l=1 + 1im)
    @test Conventions(; c_l=(-√2)).c_l === -√2

    # c_m is the spin-phase factor e^(iΘ), so it must have unit modulus.
    @test_throws ArgumentError Conventions(; c_m=2)
    @test_throws ArgumentError Conventions(; c_m=0.5im)
    @test Conventions(; c_m=cis(π / 4)).c_m === cis(π / 4)
    @test Conventions(; c_m=1im).c_m === 1im
end

@testitem "Conventions: named presets match the appendix table" tags = [:unit, :fast] begin
    import Scri: Conventions, One, MinusOne

    # :SXS (alias :SpEC) is the package default.
    @test Conventions(:SXS) === Conventions()
    @test Conventions(:SpEC) === Conventions()

    # :MB — Boyle/Lehner.
    c = Conventions(:MB)
    @test c.c_s === MinusOne() && c.c_l === -√2 && c.c_h === 2
    @test c.c_R === One() && c.c_ψ === One() && c.c_σ === One() && c.c_m === One()

    # :NP — Newman–Penrose (1968); the old angle Θ = π is the factor c_m = −1.
    c = Conventions(:NP)
    @test c.c_s === MinusOne() && c.c_ψ === MinusOne()
    @test c.c_l === -√2 && c.c_m === MinusOne()

    # :ADLK — Ashtekar et al.
    c = Conventions(:ADLK)
    @test c.c_σ === MinusOne() && c.c_l === -√2 && c.c_s === One()

    # :BR — Bishop et al.
    c = Conventions(:BR)
    @test c.c_ψ === MinusOne() && c.c_l === One()

    # :C — Chandrasekhar.
    c = Conventions(:C)
    @test c.c_s === MinusOne() && c.c_ψ === MinusOne() && c.c_l === -√2

    # The precision keyword controls the type of inexact entries.
    using DoubleFloats: Double64
    @test Conventions(:MB; T=Double64).c_l === -sqrt(Double64(2))

    # Unknown names are rejected.
    @test_throws ArgumentError Conventions(:nope)
end

@testitem "Conventions: conversion factors (values and One-elision)" tags = [:unit, :fast] begin
    import Scri:
        Conventions,
        One,
        MinusOne,
        dyad_factor,
        weyl_factor,
        faraday_factor,
        shear_factor,
        strain_factor,
        news_factor,
        conversion_factor

    # At the default conventions every factor is the One() singleton, so it costs nothing.
    c₀ = Conventions()
    @test dyad_factor(c₀) === One()
    for n ∈ 0:4
        @test weyl_factor(c₀, n) === One()
    end
    for n ∈ 0:2
        @test faraday_factor(c₀, n) === One()
    end
    for ℐ ∈ (+1, -1)
        @test shear_factor(c₀, ℐ) === One()
    end
    @test strain_factor(c₀) === One()
    @test news_factor(c₀) === One()

    # A generic convention: check every factor against its defining formula, written
    # in export form q^[X] = F q^[SXS].
    c_l, c_m, c_h = -√2, cis(π / 4), 3 + 4im
    X = Conventions(; c_s=-1, c_R=-1, c_ψ=-1, c_σ=-1, c_φ=-1, c_l, c_m, c_h)
    q = c_l * c_m
    @test dyad_factor(X) ≈ q
    for n ∈ 0:4
        @test weyl_factor(X, n) ≈ (-1)^3 * q^(2 - n)  # c_s c_ψ c_R = (−1)³
    end
    for n ∈ 0:2
        @test faraday_factor(X, n) ≈ -q^(1 - n)       # c_φ = −1
    end
    # F_σ = c_σ c_l c_m² on ℐ⁺ — no c_s (two raised m legs give c_s² = 1, and the
    # defining relation fixes the one-form l_b directly).  On ℐ⁻ the radiative shear is
    # built from the n leg, so c_l inverts.
    @test shear_factor(X, +1) ≈ -c_l * c_m^2
    @test shear_factor(X, -1) ≈ -c_m^2 / c_l
    # F_h = c_s c_h⁻¹ c_m⁻²; the News inherits it exactly (N = ∂ᵤh, shared coordinates).
    @test strain_factor(X) ≈ -1 / (c_h * c_m^2)
    @test news_factor(X) == strain_factor(X)

    # conversion_factor dispatches to the right helper for every component.
    for (S, expected) ∈ (
        (:ψ₀, weyl_factor(X, 0)),
        (:ψ₃, weyl_factor(X, 3)),
        (:φ₂, faraday_factor(X, 2)),
        (:h, strain_factor(X)),
        (:News, news_factor(X)),
    )
        @test conversion_factor(X, Val(S), 1) == expected
    end
    @test conversion_factor(X, Val(:σ), -1) == shear_factor(X, -1)

    # Out-of-range component indices are rejected.
    @test_throws ArgumentError weyl_factor(X, 5)
    @test_throws ArgumentError faraday_factor(X, 3)
end

@testitem "Conventions: non-±1 factors keep their input type (no floating)" tags = [
    :unit, :fast
] begin
    import Scri: Conventions
    using DoubleFloats: Double64

    # Integers stay integers, extended precision stays exact — nothing is silently widened.
    @test Conventions(; c_l=2).c_l === 2
    @test Conventions(; c_h=3 + 4im).c_h === 3 + 4im
    d = sqrt(Double64(2))
    @test Conventions(; c_l=d).c_l === d
end
