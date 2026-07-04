"""
    Conventions{S,R,L,M,Ψ,Σ,H,Φ,Ð,A}

Sign and scale factors setting a choice of conventions for asymptotic field quantities,
following the conventions appendix (sign/scale parameters ``s₀, s₁, s₂, s₃, λ, Θ, ζ`` there,
renamed here to mnemonic ``c`` parameters).  Each field is a ``c``-parameter; the defaults
reproduce the **SXS** conventions.

Construct with keywords (validated), or by name from the table of published conventions:

    Conventions()                   # SXS + Newman–Penrose ð (the package default)
    Conventions(; c_s=-1, c_l=-√2)  # a custom convention

"""
struct Conventions{S,R,L,M,Ψ,Σ,H,Φ,Ð,A}
    c_s::S
    c_R::R
    c_l::L
    c_m::M
    c_Ψ::Ψ
    c_σ::Σ
    c_h::H
    c_φ::Φ
    c_ð::Ð
    c_α::A
end

function Conventions(; c_s=1, c_R=1, c_l=1, c_m=1, c_Ψ=1, c_σ=1, c_h=1, c_φ=1, c_ð=1, c_α=1)
    for (name, c) ∈
        (("c_s", c_s), ("c_R", c_R), ("c_Ψ", c_Ψ), ("c_σ", c_σ), ("c_φ", c_φ), ("c_α", c_α))
        c == 1 || c == -1 || throw(ArgumentError("$name is a sign and must be ±1; got $c"))
    end
    (iszero(c_l) || iszero(c_m) || iszero(c_h) || iszero(c_ð)) &&
        throw(ArgumentError("c_l, c_m, c_h, and c_ð must be nonzero"))
    # Exact ±1 → singleton (so common conventions elide); every other value is kept in its input
    # type verbatim, with no float conversion — extended-precision inputs stay exact, and it is
    # the caller's job to combine them only through precision-preserving operations (e.g. `x / y`
    # rather than `x * y^-1`, since an integer base to a negative power throws a `DomainError`).
    return Conventions(
        signify(c_s),
        signify(c_R),
        signify(c_l),
        signify(c_m),
        signify(c_Ψ),
        signify(c_σ),
        signify(c_h),
        signify(c_φ),
        signify(c_ð),
        signify(c_α),
    )
end
