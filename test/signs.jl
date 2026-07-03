# Tests for the One/MinusOne sign singletons (src/signs.jl)

@testitem "Signs: multiplication elides at the type level" tags = [:unit, :fast] begin
    import Scri: One, MinusOne

    # `One` is a true multiplicative identity and `MinusOne` a true negation: the result is the
    # very same object (or its negation), so the sign vanishes from the generated code.
    for x ∈ (3, 3.0, 3.0 + 4.0im, 3 + 4im, -2.5)
        @test One() * x === x
        @test x * One() === x
        @test MinusOne() * x === -x
        @test x * MinusOne() === -x
        @test x / One() === x
        @test x / MinusOne() === -x
    end
end

@testitem "Signs: composition stays in the singletons" tags = [:unit, :fast] begin
    import Scri: One, MinusOne

    @test One() * One() === One()
    @test One() * MinusOne() === MinusOne()
    @test MinusOne() * One() === MinusOne()
    @test MinusOne() * MinusOne() === One()
    @test -One() === MinusOne()
    @test -MinusOne() === One()
    @test inv(One()) === One()
    @test inv(MinusOne()) === MinusOne()
    # A chain of signs collapses to a single singleton.
    @test One() * MinusOne() * MinusOne() * One() === One()
end

@testitem "Signs: behave as ordinary ±1 in generic code" tags = [:unit, :fast] begin
    import Scri: One, MinusOne
    using DoubleFloats: Double64

    @test One() == 1
    @test MinusOne() == -1
    @test !(MinusOne() == 1)
    @test isone(One()) && !isone(MinusOne())
    @test !iszero(One()) && !iszero(MinusOne())
    @test isreal(One()) && isreal(MinusOne())
    @test abs(MinusOne()) === One()
    @test sign(MinusOne()) === MinusOne()
    @test isapprox(One(), 1.0) && isapprox(MinusOne(), -1.0)
    for T ∈ (Int, Float64, Double64, ComplexF64)
        @test convert(T, One()) == one(T)
        @test convert(T, MinusOne()) == -one(T)
    end
end

@testitem "Signs: no method ambiguities across the numeric types" tags = [:unit, :fast] begin
    using Scri: Scri
    using Test: detect_ambiguities

    # Only ambiguities involving our own sign singletons are this refactor's responsibility
    # (the codebase may carry other, unrelated ambiguities).
    involves_sign(m) = occursin(r"\b(One|MinusOne)\b", string(m.sig))
    sign_ambig = filter(
        a ->
            (a[1].module === Scri || a[2].module === Scri) &&
            (involves_sign(a[1]) || involves_sign(a[2])),
        detect_ambiguities(Scri),
    )
    @test isempty(sign_ambig)
end

@testitem "Signs: mixed products are type-stable and allocation-free" tags = [:unit, :fast] begin
    import Scri: One, MinusOne

    negate(x) = MinusOne() * x
    identish(x) = One() * x
    z = 3.0 + 4.0im
    @test @inferred(negate(z)) === -z
    @test @inferred(identish(z)) === z
    negate(z), identish(z)                       # warm up
    @test @allocated(negate(z)) == 0
    @test @allocated(identish(z)) == 0
end
