# Tests for the Conventions struct (src/conventions.jl)

@testitem "Conventions: defaults are SpEC + Newman–Penrose ð" tags = [:unit, :fast] begin
    import Scri: Conventions, One

    c = Conventions()
    @test c.c_s == 1
    @test c.c_R == 1
    @test c.c_Ψ == 1
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
    @test_throws ArgumentError Conventions(; c_Ψ=0)
    @test_throws ArgumentError Conventions(; c_α=2)        # c_α is a sign too
    @test_throws ArgumentError Conventions(; c_l=0)        # factors must be nonzero
    @test_throws ArgumentError Conventions(; c_m=0)        # c_m is a (nonzero) factor now
    @test_throws ArgumentError Conventions(; c_h=0)
    @test_throws ArgumentError Conventions(; c_ð=0)
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
