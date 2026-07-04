"""
    Conventions{S,R,L,M,Ψ,Σ,H,Φ,Ð,A}

Sign and scale factors setting a choice of conventions for asymptotic field quantities,
following the conventions appendix (sign/scale parameters ``s₀, s₁, s₂, s₃, λ, Θ, ζ`` there,
renamed here to mnemonic ``c`` parameters).  Each field is a ``c``-parameter; the defaults
reproduce the **SXS** conventions.

Every parameter is defined in **export form**: it takes a quantity from the package-native
SXS convention *to* the convention being described, as in ``q^{[X]} = c_q q^{[SXS]}`` (and
per tetrad leg, ``ℓ^{[X]} = c_l ℓ^{[SXS]}``, ``m^{[X]} = c_m m^{[SXS]}``, with ``n^{[X]} =
n^{[SXS]}/c_l`` forced by the ``ℓ·n`` normalization).  Because ``ℓ`` is a real null vector,
`c_l` must be real; because ``m·m̄`` is fixed, `c_m` must be a unit-modulus phase.  Both are
validated.

Construct with keywords (validated), or by name from the table of published conventions:

    Conventions()                   # SXS + Newman–Penrose ð (the package default)
    Conventions(; c_s=-1, c_l=-√2)  # a custom convention
    Conventions(:NP)                # a named preset; see the symbol-argument docstring

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
    (iszero(c_l) || iszero(c_h) || iszero(c_ð)) &&
        throw(ArgumentError("c_l, c_h, and c_ð must be nonzero"))
    isreal(c_l) || throw(
        ArgumentError("c_l scales the real null vector ℓ, so it must be real; got $c_l")
    )
    abs2(c_m) ≈ 1 || throw(
        ArgumentError(
            "c_m is the spin-phase factor e^(iΘ) of the m leg, " *
            "so it must have unit modulus; got $c_m",
        ),
    )
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

"""
    Conventions(name::Symbol; T=Float64)

Named presets from the conventions-appendix table.  Available: `:SXS` (the package default),
`:MB` (Moreschi/Boyle), `:NP` (Newman–Penrose), `:ADLK` (Ashtekar et al.), `:BR` (Bishop et
al.), `:C` (Chandrasekhar).  `T` sets the type used for inexact values such as `-√2`; exact
values (`±1`, `2`) are stored exactly regardless.  Values that are not clear from the
references are left at the default.
"""
function Conventions(name::Symbol; T::Type{<:Real}=Float64)
    s2 = sqrt(T(2))
    return if name === :SXS || name === :SpEC
        Conventions()
    elseif name === :MB
        Conventions(; c_s=-1, c_l=(-s2), c_h=2)
    elseif name === :NP
        Conventions(; c_s=-1, c_Ψ=-1, c_l=(-s2), c_m=-1)
    elseif name === :ADLK
        Conventions(; c_σ=-1, c_l=(-s2))
    elseif name === :BR
        Conventions(; c_Ψ=-1)
    elseif name === :C
        Conventions(; c_s=-1, c_Ψ=-1, c_l=(-s2))
    else
        throw(
            ArgumentError(
                "Unknown convention preset :$name; " *
                "expected one of :SXS, :SpEC, :MB, :NP, :ADLK, :BR, :C",
            ),
        )
    end
end

"""
    dyad_factor(c::Conventions)

The dyad scaling ``c_l c_m`` — the only convention combination that enters the homogeneous
(peeling-tower) transformation laws.  A transform uses the mixing parameter ``b^{[c]} = (c_l
c_m) b^{[SXS]}`` on ``ℐ⁺``, and ``b̄^{[c]} = b̄^{[SXS]} / (c_l c_m)`` on ``ℐ⁻`` (the ``ℐ⁻``
tower mixes *downward*, so the factor ratio ``Fₙ/Fₙ₋ₖ = (c_l c_m)^{-k}`` inverts the
scaling).  See the "Convention dependence" section of the "BMS action on fields"
documentation page.
"""
dyad_factor(c::Conventions) = c.c_l * c.c_m

"""
    weyl_factor(c::Conventions, n::Integer)

Factor taking the Weyl component ``ψₙ`` from the SXS convention to convention `c`:
``ψₙ^{[c]} = c_s c_Ψ c_R (c_l c_m)^{2-n} ψₙ^{[SXS]}``.
"""
function weyl_factor(c::Conventions, n::Integer)
    0 ≤ n ≤ 4 || throw(ArgumentError("Weyl component index n must be in 0:4; got $n"))
    q = dyad_factor(c)
    F = c.c_s * c.c_Ψ * c.c_R
    # Written without negative integer powers, which would throw for exact integer scales.
    return if n == 0
        F * q^2
    elseif n == 1
        F * q
    elseif n == 2
        F
    elseif n == 3
        F / q
    else
        F / q^2
    end
end

"""
    faraday_factor(c::Conventions, n::Integer)

Factor taking the Faraday component ``φₙ`` from the SXS convention to convention `c`:
``φₙ^{[c]} = c_φ (c_l c_m)^{1-n} φₙ^{[SXS]}`` (two tetrad legs instead of the Weyl four).
"""
function faraday_factor(c::Conventions, n::Integer)
    0 ≤ n ≤ 2 || throw(ArgumentError("Faraday component index n must be in 0:2; got $n"))
    q = dyad_factor(c)
    return if n == 0
        c.c_φ * q
    elseif n == 1
        c.c_φ
    else
        c.c_φ / q
    end
end

"""
    shear_factor(c::Conventions, ℐ=1)

Factor taking the asymptotic shear ``σ`` from the SXS convention to convention `c`.  On
``ℐ⁺`` the shear is built from the ``ℓ`` leg, ``σ = -c_σ mᵃmᵇ∇ₐl_b``, giving ``F_σ = c_σ c_l
c_m²`` — with **no** ``c_s``, because the two raised ``m`` indices contribute ``c_s² = 1``
and the defining relation fixes the one-form ``l_b`` directly.  On ``ℐ⁻`` the radiative
shear is instead built from the ``n`` leg (the ``λ̄``-type coefficient, still spin weight
+2), so ``c_l`` inverts: ``F_σ = c_σ c_m² / c_l``.
"""
function shear_factor(c::Conventions, ℐ::Integer=1)
    F = c.c_σ * c.c_m^2
    return ℐ == 1 ? F * c.c_l : F / c.c_l
end

"""
    strain_factor(c::Conventions)

Factor taking the strain ``h`` from the SXS convention to convention `c`: ``h^{[c]} = c_s
c_h^{-1} c_m^{-2} h^{[SXS]}``.  The strain comes from the metric perturbation, with no
``ℓ`` or ``n`` leg, so the factor is the same on both null infinities.
"""
strain_factor(c::Conventions) = c.c_s / (c.c_h * c.c_m^2)

"""
    news_factor(c::Conventions)

Factor taking the News from the SXS convention to convention `c`.  The News is ``∂ᵤ`` of the
strain, and the coordinates (including ``u``) are shared by all conventions — only tetrad
and field *definitions* differ — so the News inherits the strain factor exactly.
"""
news_factor(c::Conventions) = strain_factor(c)

"""
    conversion_factor(c::Conventions, ::Val{S}, ℐ=1)

Factor taking data component `S` from the SXS convention to convention `c`, as in
``q^{[c]} = \\texttt{conversion\\_factor}(c, \\mathrm{Val}(S), ℐ)\\, q^{[SXS]}``.  The
general convention-`X`-to-convention-`Y` factor is the ratio ``F_Y / F_X`` of two of these.
Only the shear factor depends on `ℐ`; see [`shear_factor`](@ref).
"""
conversion_factor(c::Conventions, ::Val{:ψ₀}, ℐ::Integer=1) = weyl_factor(c, 0)
conversion_factor(c::Conventions, ::Val{:ψ₁}, ℐ::Integer=1) = weyl_factor(c, 1)
conversion_factor(c::Conventions, ::Val{:ψ₂}, ℐ::Integer=1) = weyl_factor(c, 2)
conversion_factor(c::Conventions, ::Val{:ψ₃}, ℐ::Integer=1) = weyl_factor(c, 3)
conversion_factor(c::Conventions, ::Val{:ψ₄}, ℐ::Integer=1) = weyl_factor(c, 4)
conversion_factor(c::Conventions, ::Val{:φ₀}, ℐ::Integer=1) = faraday_factor(c, 0)
conversion_factor(c::Conventions, ::Val{:φ₁}, ℐ::Integer=1) = faraday_factor(c, 1)
conversion_factor(c::Conventions, ::Val{:φ₂}, ℐ::Integer=1) = faraday_factor(c, 2)
conversion_factor(c::Conventions, ::Val{:σ}, ℐ::Integer=1) = shear_factor(c, ℐ)
conversion_factor(c::Conventions, ::Val{:h}, ℐ::Integer=1) = strain_factor(c)
conversion_factor(c::Conventions, ::Val{:News}, ℐ::Integer=1) = news_factor(c)
