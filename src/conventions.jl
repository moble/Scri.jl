"""
    Conventions{S,R,L,M,Ψ,Σ,H,Φ,Ð,A}

Sign and scale factors setting conventions.  The parameters are

  - `c_s` for the signature of the metric
  - `c_R` for the Riemann tensor
  - `c_l` for the tetrad vector `l`
  - `c_m` for the tetrad vector `m`
  - `c_ψ` for the Weyl components
  - `c_σ` for the shear
  - `c_h` for the complex strain
  - `c_φ` for the Faraday components
  - `c_ð` for the ð operator
  - `c_α` for the supertranslation law

For precise definitions of each of these factors, see the [documentation](@ref
conventions-overview).  The default values are all 1, corresponding to the SXS conventions.

We constrain the possible values of these factors.  The signature `c_s`, Riemann definition
`c_R`, and supertranslation law `c_α` only ever differ by signs, so they must be ±1.  The
tetrad vector `l` is real, so `c_l` must be real.  The tetrad vector `m` is complex, but
must obey ``m ⋅ m̄ = c_s``, so `c_m` must have unit modulus (approximately).  And all of
these factors must be nonzero, but there are so many bizarre conventions that we permit them
to be complex if desired — though they will also usually just be signs.

Note that we *do* assume various structural conventions, such as ``l`` being the outgoing
null vector (so that ``ψ₄`` is the radiative field at future null infinity), ``m`` having a
certain handedness, and ``h`` having spin weight ``-2``.  Any users wishing to use different
structural conventions will have to handle the conversions themselves.

The types of all of these factors are retained as the type parameters of `Conventions`
listed above — `S,R,L,M,Ψ,Σ,H,Φ,Ð,A` — to allow for type stability and compile-time
optimization.  Most common conventions will compile away entirely.

Construct with keywords (validated), or by name from the table of published conventions:

    Conventions()                   # SXS defaults
    Conventions(; c_s=-1, c_l=-√2)  # a custom convention
    Conventions(:NP)                # a named preset; see the symbol-argument docstring

"""
struct Conventions{S,R,L,M,Ψ,Σ,H,Φ,Ð,A}
    c_s::S
    c_R::R
    c_l::L
    c_m::M
    c_ψ::Ψ
    c_σ::Σ
    c_h::H
    c_φ::Φ
    c_ð::Ð
    c_α::A
end

function Conventions(; c_s=1, c_R=1, c_l=1, c_m=1, c_ψ=1, c_σ=1, c_h=1, c_φ=1, c_ð=1, c_α=1)
    for (name, c) ∈ (("c_s", c_s), ("c_R", c_R), ("c_α", c_α))
        if c ≠ 1 && c ≠ -1
            throw(ArgumentError("$name is a sign and must be ±1; got $c"))
        end
    end
    for (name, c) ∈ (("c_l", c_l),)
        if !isreal(c)
            throw(ArgumentError("$name must be real; got $c"))
        end
    end
    if abs2(c_m) ≉ 1
        throw(
            ArgumentError(
                "c_m is the spin-phase factor e^(iΘ) of the m leg, " *
                "so it must have unit modulus; got $c_m",
            ),
        )
    end
    for (name, c) ∈ (
        ("c_s", c_s),
        ("c_R", c_R),
        ("c_l", c_l),
        ("c_m", c_m),
        ("c_ψ", c_ψ),
        ("c_σ", c_σ),
        ("c_h", c_h),
        ("c_φ", c_φ),
        ("c_ð", c_ð),
        ("c_α", c_α),
    )
        if iszero(c)
            throw(ArgumentError("$name must be nonzero; got $c"))
        end
    end
    # Exact ±1 → singleton (so common conventions elide); every other value is kept in its input
    # type verbatim, with no float conversion — extended-precision inputs stay exact, and it is
    # the caller's job to combine them only through precision-preserving operations (e.g. `x / y`
    # rather than `x * y^-1`, since an integer base to a negative power throws a `DomainError`).
    return Conventions(
        signify(c_s),
        signify(c_R),
        signify(c_l),
        signify(c_m),
        signify(c_ψ),
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
        Conventions(; c_s=-1, c_ψ=-1, c_l=(-s2), c_m=-1)
    elseif name === :ADLK
        Conventions(; c_σ=-1, c_l=(-s2))
    elseif name === :BR
        Conventions(; c_ψ=-1)
    elseif name === :C
        Conventions(; c_s=-1, c_ψ=-1, c_l=(-s2))
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
``ψₙ^{[c]} = c_s c_ψ c_R (c_l c_m)^{2-n} ψₙ^{[SXS]}``.
"""
function weyl_factor(c::Conventions, n::Integer)
    0 ≤ n ≤ 4 || throw(ArgumentError("Weyl component index n must be in 0:4; got $n"))
    q = dyad_factor(c)
    F = c.c_s * c.c_ψ * c.c_R
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
