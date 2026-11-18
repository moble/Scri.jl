# Tests for the Conventions struct (src/conventions.jl)

@testitem "Conventions: defaults are SXS" tags = [:unit, :fast] begin
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
    @test c.c_N == 1
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
    @test_throws ArgumentError Conventions(; c_λ=0)
    @test_throws ArgumentError Conventions(; c_h=0)
    @test_throws ArgumentError Conventions(; c_N=0)
    @test_throws ArgumentError Conventions(; c_ð=0)

    # The code always uses c_ð = 1, so the keyword constructor warns about any other value,
    # but the presets, which record a non-unit c_ð on purpose, do not.
    @test_logs (:warn, r"c_ð") Conventions(; c_ð=1 / √2)
    @test_logs Conventions(; c_ð=1)
    @test_logs Conventions(:MB)
    @test_logs Conventions(:scri2024abd)
    @test_logs Conventions(:scri2024wm)

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

@testitem "Conventions: named presets have the documented values" tags = [:unit, :fast] begin
    import Scri: Conventions, One, MinusOne

    # :SXS is the package default.
    @test Conventions(:SXS) === Conventions()

    # :MB — Moreschi/Boyle.  See the "Tetrad" page footnote: lowering the Boyle (2015)
    # cotetrad gives c_s = -1 and c_l = +√2 exactly, with the m-leg sign (relative to
    # the polar spin frame) booked as c_m = -1.
    c = Conventions(:MB)
    @test c.c_s === MinusOne() && c.c_l === √2 && c.c_m === MinusOne()
    @test c.c_R === One() && c.c_ψ === MinusOne() && c.c_σ === MinusOne()
    @test c.c_h === -0.5 && c.c_ð === 1 / √2 && c.c_λ === One() && c.c_φ === One()
    @test c.c_N === MinusOne()

    # The precision keyword controls the type of every entry other than ±1.
    using DoubleFloats: Double64
    c64 = Conventions(:MB; T=Double64)
    @test c64.c_l === sqrt(Double64(2)) && c64.c_h === -inv(Double64(2))
    @test c64.c_ð === inv(sqrt(Double64(2)))

    # :scri2024abd — scri's AsymptoticBondiData, which is :MB but for the strain and News.
    c = Conventions(:scri2024abd)
    @test c.c_s === MinusOne() && c.c_l === √2 && c.c_m === MinusOne()
    @test c.c_R === One() && c.c_ψ === MinusOne() && c.c_σ === MinusOne()
    @test c.c_h === MinusOne() && c.c_N === One()
    @test c.c_ð === 1 / √2 && c.c_λ === One() && c.c_φ === One()
    c64 = Conventions(:scri2024abd; T=Double64)
    @test c64.c_l === sqrt(Double64(2)) && c64.c_h === MinusOne()

    # :scri2024wm — scri's WaveformModes, which is SXS but for the shear.
    c = Conventions(:scri2024wm)
    @test c.c_σ === √2
    for f ∈ fieldnames(Conventions)
        f === :c_σ || @test getfield(c, f) === One()
    end
    @test Conventions(:scri2024wm; T=Double64).c_σ === sqrt(Double64(2))

    # Unknown names are rejected.
    @test_throws ArgumentError Conventions(:nope)
end

@testitem "Conventions: the :MB preset follows Boyle (2015)" tags = [:unit, :fast] begin
    import Scri: Conventions, weyl_factor, shear_factor, strain_factor, news_factor

    # Boyle (2015) defines σ in spinor form (Eq. 18), which is m m ∇l in its +--- metric,
    # so that σ = h̄/2 in terms of the SXS strain.
    # The SXS shear is σ = h̄/2√2, so the MB shear is √2 times the SXS one.  The paper's
    # ψ₄ = -∂²ᵤσ̄ (Eq. 20) then gives ψₙ factors 2, -√2, 1, -1/√2, and 1/2, which are also
    # the factors that the Python package scri applies in `create_abd_from_h5`.
    c = Conventions(:MB)
    for (n, factor) ∈ zip(0:4, (2, -√2, 1, -1 / √2, 1 / 2))
        @test weyl_factor(c, n) ≈ factor
    end
    @test shear_factor(c) ≈ √2

    # The paper defines the strain as h = σ̄, which is then half the SXS strain, and the
    # news as n = -∂ᵤσ̄ (Eq. 22), which is minus the time derivative of that strain.
    @test strain_factor(c) ≈ 1 / 2
    @test news_factor(c) ≈ -1 / 2

    # Hence ψ₄ = -ḧ in both conventions, as Eq. 20 has with h = σ̄, and the Schwarzschild
    # ψ₂ = -M is the same in both, as the mass aspect -Re ψ₂ of scri needs.
    @test weyl_factor(c, 4) ≈ strain_factor(c)
    @test weyl_factor(c, 2) == 1
end

@testitem "Conventions: the :scri2024abd and :scri2024wm presets match scri 2024" tags = [
    :unit, :fast
] begin
    import Scri: Conventions, One, weyl_factor, shear_factor, strain_factor, news_factor

    # The AsymptoticBondiData of scri 2024 holds σ and ψₙ as Boyle (2015) defines them, with
    # the same factors as :MB, but its strain `h` is 2σ̄, the SXS strain, and the news it
    # computes is 2σ̄̇ = ∂ᵤh, the SXS News.  Both factors are therefore exactly one.
    c = Conventions(:scri2024abd)
    mb = Conventions(:MB)
    for n ∈ 0:4
        @test weyl_factor(c, n) == weyl_factor(mb, n)
    end
    @test shear_factor(c) == shear_factor(mb) ≈ √2
    @test strain_factor(c) === One()
    @test news_factor(c) === One()

    # So ψ₄ = -ḧ/2 here, as scri's constraint ψ₄ = -σ̄̈ with h = 2σ̄ requires.
    @test weyl_factor(c, 4) ≈ strain_factor(c) / 2

    # A WaveformModes of scri 2024 holds SXS data, except that its `sigma` is the
    # AsymptoticBondiData shear h̄/2, which scri shifts by ð_GHP²α = ð²α/2.
    w = Conventions(:scri2024wm)
    for n ∈ 0:4
        @test weyl_factor(w, n) === One()
    end
    @test strain_factor(w) === One()
    @test news_factor(w) === One()
    @test shear_factor(w) == shear_factor(c)
    # The shift of σ in the transformation law is F_σ ð²α/2√2, which is ð²α/2 here.
    @test shear_factor(w) / (2√2) ≈ 1 / 2
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
        lambda_factor,
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
    @test shear_factor(c₀) === One()
    @test lambda_factor(c₀) === One()
    @test strain_factor(c₀) === One()
    @test news_factor(c₀) === One()

    # A generic convention: check every factor against its defining formula, written
    # in export form q^[X] = F q^[SXS].
    c_l, c_m, c_h, c_N, c_λ = -√2, cis(π / 4), 3 + 4im, 2 - 1im, -1
    X = Conventions(; c_s=-1, c_R=-1, c_ψ=-1, c_σ=-1, c_λ, c_φ=-1, c_l, c_m, c_h, c_N)
    q = c_l * c_m
    @test dyad_factor(X) ≈ q
    for n ∈ 0:4
        @test weyl_factor(X, n) ≈ (-1)^3 * q^(2 - n)  # c_s c_ψ c_R = (−1)³
    end
    for n ∈ 0:2
        @test faraday_factor(X, n) ≈ q^(1 - n)        # c_s c_φ = (−1)²
    end
    # The signature alone flips every φₙ, as it flips every ψₙ: both are contractions of a
    # tensor lowered by one index from a mixed tensor that all conventions share.
    Y = Conventions(; c_s=-1)
    for n ∈ 0:2
        @test faraday_factor(Y, n) == -1 == weyl_factor(Y, n)
    end
    # F_σ = c_s c_σ c_l c_m² (ℐ⁺ shear, of the l congruence): m legs are vectors (c_m²,
    # no metric), the lowered l_b = g l contributes c_s c_l.  Here c_s c_σ = (−1)² = 1.
    @test shear_factor(X) ≈ c_l * c_m^2
    # F_λ = c_s c_λ / (c_l c_m²) (ℐ⁻ shear = NP's λ, of the n congruence, n ~ 1/c_l).
    @test lambda_factor(X) ≈ (-1) * c_λ / (c_l * c_m^2)
    # F_h = c_s c_h — no tetrad leg enters h = c_h(h₊ − ih×), so no c_m; the News
    # takes the strain factor times c_N (N = c_N ∂ᵤh, shared coordinates).
    @test strain_factor(X) ≈ -c_h
    @test news_factor(X) ≈ c_N * strain_factor(X)

    # conversion_factor dispatches to the right helper for every component.
    for (S, expected) ∈ (
        (:ψ₀, weyl_factor(X, 0)),
        (:ψ₃, weyl_factor(X, 3)),
        (:φ₂, faraday_factor(X, 2)),
        (:σ, shear_factor(X)),
        (:λ, lambda_factor(X)),
        (:h, strain_factor(X)),
        (:News, news_factor(X)),
    )
        @test conversion_factor(X, Val(S)) == expected
    end

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
