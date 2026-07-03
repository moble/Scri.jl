"""
    One
    MinusOne

Singleton types representing the numbers ``+1`` and ``-1`` *at the type level*.

Their whole purpose is to make a sign that is fixed by a choice of convention disappear from the
generated code.  Because the value lives in the type, multiplication is resolved at compile time:
`One() * x` returns `x` itself (the very same object) and `MinusOne() * x` returns `-x`, so a
convention factor of ``±1`` costs nothing (at most a negation) rather than a runtime multiply.
Products of these singletons stay singletons (`MinusOne() * MinusOne() === One()`), so a chain
like ``c_s c_Ψ c_R`` collapses to a single `One`/`MinusOne` before it ever touches field data.

They deliberately do **not** subtype `Number`.  A ``±1`` newtype that did would have to define
`*(::One, ::Number)` and `convert(::Type{<:Number}, ::One)`, which then clash — irreconcilably,
across dozens of ambiguities — with every other package that adds its own `Number` subtype
(quaternions, static integers, extended-precision floats, …).  Staying outside the `Number`
hierarchy means those downstream methods simply never apply to a singleton, so the identity and
negation methods below win cleanly.  The small numeric interface they *do* need (`==`, `≈`, `^`,
conversion, …) is supplied here by concrete dispatch.

These are internal helpers used to store the sign fields of [`Conventions`](@ref); values are
normalized to them by the `Conventions` constructor, so user code rarely names them directly.
"""
struct One end

@doc (@doc One) struct MinusOne end

const SignSingleton = Union{One,MinusOne}

# Underlying integer value; used only where an actual number is unavoidable (`+`, `-`, `≈`).
signval(::One) = 1
signval(::MinusOne) = -1

# Identity / negation multiply — the point of the whole exercise.  These return the other operand
# unchanged (or negated), so the singleton vanishes at the Julia type level.  `x` ranges over
# `Number`, which the singletons themselves are not, so there is no self-overlap here.
Base.:*(::One, x::Number) = x
Base.:*(x::Number, ::One) = x
Base.:*(::MinusOne, x::Number) = -x
Base.:*(x::Number, ::MinusOne) = -x

# Composition: a product of signs is again a sign, resolved at compile time.
Base.:*(::One, ::One) = One()
Base.:*(::One, ::MinusOne) = MinusOne()
Base.:*(::MinusOne, ::One) = MinusOne()
Base.:*(::MinusOne, ::MinusOne) = One()

# Division and inverse (needed by e.g. the strain factor's ``c_h^{-1}`` and ratios of Weyl
# factors); ``±1`` is its own inverse, and dividing by/into it is the same elision.
Base.:/(x::Number, ::One) = x
Base.:/(x::Number, ::MinusOne) = -x
Base.:/(::One, x::Number) = inv(x)
Base.:/(::MinusOne, x::Number) = -inv(x)
Base.:/(::One, ::One) = One()
Base.:/(::One, ::MinusOne) = MinusOne()
Base.:/(::MinusOne, ::One) = MinusOne()
Base.:/(::MinusOne, ::MinusOne) = One()
Base.inv(x::SignSingleton) = x

# Unary minus flips the sign, staying in the singleton world.
Base.:-(::One) = MinusOne()
Base.:-(::MinusOne) = One()

# Integer powers stay in the singletons too.
Base.:^(::One, ::Integer) = One()
Base.:^(::MinusOne, p::Integer) = iseven(p) ? One() : MinusOne()

# Sums/differences of two signs fall back to plain integers.
Base.:+(a::SignSingleton, b::SignSingleton) = signval(a) + signval(b)
Base.:-(a::SignSingleton, b::SignSingleton) = signval(a) - signval(b)

# Enough of a numeric interface to be well-behaved in generic code, all by concrete dispatch.
Base.one(::SignSingleton) = One()
Base.one(::Type{<:SignSingleton}) = One()
Base.oneunit(::SignSingleton) = One()
Base.oneunit(::Type{<:SignSingleton}) = One()
Base.zero(::SignSingleton) = 0
Base.zero(::Type{<:SignSingleton}) = 0
Base.isone(::One) = true
Base.isone(::MinusOne) = false
Base.iszero(::SignSingleton) = false
Base.isfinite(::SignSingleton) = true
Base.isnan(::SignSingleton) = false
Base.isinf(::SignSingleton) = false
Base.sign(x::SignSingleton) = x
Base.abs(::SignSingleton) = One()
Base.isreal(::SignSingleton) = true
Base.real(x::SignSingleton) = x
Base.imag(::SignSingleton) = 0
Base.conj(x::SignSingleton) = x

Base.:(==)(::One, y::Number) = isone(y)
Base.:(==)(y::Number, ::One) = isone(y)
Base.:(==)(::MinusOne, y::Number) = y == -1
Base.:(==)(y::Number, ::MinusOne) = y == -1
Base.isapprox(a::SignSingleton, b::Number; kw...) = isapprox(signval(a), b; kw...)
Base.isapprox(a::Number, b::SignSingleton; kw...) = isapprox(a, signval(b); kw...)

Base.convert(::Type{T}, ::One) where {T<:Number} = one(T)
Base.convert(::Type{T}, ::MinusOne) where {T<:Number} = -one(T)
Base.promote_rule(::Type{<:SignSingleton}, ::Type{T}) where {T<:Number} = T

Base.show(io::IO, ::One) = print(io, "One()")
Base.show(io::IO, ::MinusOne) = print(io, "MinusOne()")

"""
    signify(x)

Normalize a convention value to its most specialized representation: an exact ``+1`` becomes
`One()`, an exact ``-1`` becomes `MinusOne()`, and anything else is returned unchanged.  This is
how the [`Conventions`](@ref) constructor turns the common ``±1`` factors into type-level
constants while leaving genuine numbers (`-√2`, `1/√2`, a complex phase, …) as they are.
"""
function signify(x)
    if x == 1
        One()
    elseif x == -1
        MinusOne()
    else
        x
    end
end
signify(x::SignSingleton) = x
