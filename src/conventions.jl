"""
    Conventions{S,R,L,M,Ψ,Σ,Λ,H,Φ,Ð,A}

Sign and scale factors setting conventions.  The parameters are

  - `c_s` for the signature of the metric
  - `c_R` for the Riemann tensor
  - `c_l` for the tetrad vector `l`
  - `c_m` for the tetrad vector `m`
  - `c_ψ` for the Weyl components
  - `c_σ` for the shear
  - `c_λ` for the "incoming" shear
  - `c_h` for the complex strain
  - `c_φ` for the Faraday components
  - `c_ð` for the ð operator
  - `c_α` for the supertranslation law

For precise definitions of each of these factors, their allowed values, and the structural
conventions this package fixes once and for all (e.g. ``l`` outgoing, ``m``'s handedness,
``h`` of spin weight ``-2``), see the [documentation](@ref conventions-overview).  The
default values are all 1, corresponding to the SXS conventions.  The constructor validates
the values: `c_s`, `c_R`, and `c_α` must be ``±1``; `c_l` must be real; `c_m` must have
unit modulus; and every factor must be nonzero.

The types of all of these factors are retained as the type parameters of `Conventions`
listed above — `S,R,L,M,Ψ,Σ,Λ,H,Φ,Ð,A` — to allow for type stability and compile-time
optimization.  Most common conventions will compile away entirely.

Construct with keywords (validated), or by name from the table of published conventions:

    Conventions()                   # SXS defaults
    Conventions(; c_s=-1, c_l=-√2)  # a custom convention
    Conventions(:MB)                # a named preset; see the symbol-argument docstring

"""
struct Conventions{S,R,L,M,Ψ,Σ,Λ,H,Φ,Ð,A}
    c_s::S
    c_R::R
    c_l::L
    c_m::M
    c_ψ::Ψ
    c_σ::Σ
    c_λ::Λ
    c_h::H
    c_φ::Φ
    c_ð::Ð
    c_α::A
end

function Conventions(;
    c_s=1, c_R=1, c_l=1, c_m=1, c_ψ=1, c_σ=1, c_λ=1, c_h=1, c_φ=1, c_ð=1, c_α=1
)
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
        ("c_λ", c_λ),
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
        signify(c_λ),
        signify(c_h),
        signify(c_φ),
        signify(c_ð),
        signify(c_α),
    )
end

"""
    Conventions(name::Symbol; T=Float64)

Named presets.  Available: `:SXS` (the package default) and `:MB` (Moreschi/Boyle).  `T`
sets the type used for inexact values such as `√2`; exact values (`±1`, `2`) are stored
exactly regardless.  Values that are not clear from the references are left at the default.
"""
function Conventions(name::Symbol; T::Type{<:Real}=Float64)
    s2 = sqrt(T(2))
    return if name === :SXS
        Conventions()
    elseif name === :MB
        # See the footnote on the "Tetrad" conventions page: lowering the cotetrad of
        # Boyle (2015) identifies c_s = -1 and c_l = +√2 exactly, with the sign of the
        # m leg (relative to the polar spin frame) booked as c_m = -1.
        Conventions(; c_s=-1, c_l=s2, c_m=-1, c_h=2)
    else
        throw(
            ArgumentError("Unknown convention preset :$name; expected either :SXS or :MB")
        )
    end
end

"""
    dyad_factor(c::Conventions)

The dyad scaling ``c_l c_m`` — the only convention combination that enters the homogeneous
(peeling-tower) transformation laws.  Because each rung of a tower mixes neighboring
components, the mixing parameter converts between conventions by the ratio of neighboring
conversion factors — a purely multiplicative rule; no convention factor is ever conjugated.
On ``ℐ⁺`` the tower mixes *upward* (``ψₙ`` picks up ``ψₙ₊ₖ``), and ``Fₙ/Fₙ₊ₖ = (c_l c_m)^k``
gives ``b^{[c]} = (c_l c_m) b^{[SXS]}``, where ``b^{[SXS]} = ðu′/2κ``.  On ``ℐ⁻`` the tower
mixes *downward* (``ψₙ`` picks up ``ψₙ₋ₖ``), and ``Fₙ/Fₙ₋ₖ = (c_l c_m)^{-k}`` gives
``b̄^{[c]} = b̄^{[SXS]} / (c_l c_m)``, where ``b̄^{[SXS]} = ð̄v′/2κ``.  The bar is part of
the *name* of the ``ℐ⁻`` parameter: its SXS value is the conjugate of the SXS ``ℐ⁺`` value,
but ``b̄^{[c]}`` is **not** the conjugate of ``b^{[c]}`` — conjugating ``c_l c_m`` would
leave ``c_l`` upstairs.  See the "Convention dependence" section of the "BMS action on
fields" documentation page.
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
    shear_factor(c::Conventions)

Factor taking the asymptotic shear ``σ`` from the SXS convention to convention `c`.  This is
the ``ℐ⁺`` radiative shear, built from the outgoing ``l`` leg,

```math
σ = -c_σ mᵃmᵇ∇ₐl_b.
```

The tetrad legs are defined as *vectors* (see the "Tetrad" conventions page), so the two
``m`` legs contribute ``c_m²``, the covariant derivative is convention-independent (the
Christoffel symbols do not change), and the lowered ``l_b = g_{bc} l^c`` contributes ``c_s
c_l`` — one power of the metric survives.  Thus ``F_σ = c_s c_σ c_l c_m²``.

At ``ℐ⁻`` the relevant shear is a *different* NP spin coefficient, ``λ`` — the shear of the
``n`` congruence — with its own convention parameter; see [`lambda_factor`](@ref).
"""
shear_factor(c::Conventions) = c.c_s * c.c_σ * c.c_l * c.c_m^2

"""
    lambda_factor(c::Conventions)

Factor taking the ``ℐ⁻`` radiative shear ``λ`` from the SXS convention to convention `c`.
Just as the radiative Weyl component switches from ``ψ₄`` (on ``ℐ⁺``) to ``ψ₀`` (on
``ℐ⁻``), the radiative shear switches from ``σ`` to the NP coefficient ``λ`` — the shear of
the incoming ``n`` congruence, which is the one transverse to ``ℐ⁻``.  Following Moxon (and
thus SXS), it is defined *without* the minus sign that ``σ`` carries,

```math
λ = c_λ m̄ᵃm̄ᵇ∇ₐn_b        (spin weight -2).
```

Because the ``σ`` and ``λ`` definitions are structurally asymmetric (the minus in ``σ`` but
not in ``λ``), ``c_λ`` is an independent convention parameter, *not* tied to `c_σ`.  The two
``m̄`` legs are vectors contributing ``c̄_m² = 1/c_m²``, the covariant derivative is
convention-independent, and the lowered ``n_b = g_{bc} n^c`` contributes ``c_s / c_l`` (the
``n`` leg scales as ``1/c_l``, and lowering costs one ``c_s``).  Thus ``F_λ = c_s c_λ / (c_l
c_m²)``.
"""
lambda_factor(c::Conventions) = c.c_s * c.c_λ / (c.c_l * c.c_m^2)

"""
    strain_factor(c::Conventions)

Factor taking the strain ``h`` from the SXS convention to convention `c`: ``h^{[c]} = c_s
c_h h^{[SXS]}``.  The defining relation is ``h = c_h (h₊ - i h_×)``, and the orthonormal
polarization components inherit the signature flip of the metric perturbation
(``h_{ab}^{[c]} = c_s h_{ab}``, since both ``g_{ab}`` and ``η_{ab}`` flip), so ``F_h = c_s
c_h``.  No tetrad leg enters the definition, so there is no ``c_l`` or ``c_m`` dependence,
and the factor is the same on both null infinities.
"""
strain_factor(c::Conventions) = c.c_s * c.c_h

"""
    news_factor(c::Conventions)

Factor taking the News from the SXS convention to convention `c`.  The News is ``∂ᵤ`` of the
strain, and the coordinates (including ``u``) are shared by all conventions — only tetrad
and field *definitions* differ — so the News inherits the strain factor exactly.
"""
news_factor(c::Conventions) = strain_factor(c)

"""
    conversion_factor(c::Conventions, ::Val{S})

Factor taking data component `S` from the SXS convention to convention `c`, as in
``q^{[c]} = \\texttt{conversion\\_factor}(c, \\mathrm{Val}(S))\\, q^{[SXS]}``.  The
general convention-`X`-to-convention-`Y` factor is the ratio ``F_Y / F_X`` of two of these.
The radiative shear is ``σ`` on ``ℐ⁺`` and the distinct coefficient ``λ`` on ``ℐ⁻``, each
with its own factor; see [`shear_factor`](@ref) and [`lambda_factor`](@ref).
"""
conversion_factor(c::Conventions, ::Val{:ψ₀}) = weyl_factor(c, 0)
conversion_factor(c::Conventions, ::Val{:ψ₁}) = weyl_factor(c, 1)
conversion_factor(c::Conventions, ::Val{:ψ₂}) = weyl_factor(c, 2)
conversion_factor(c::Conventions, ::Val{:ψ₃}) = weyl_factor(c, 3)
conversion_factor(c::Conventions, ::Val{:ψ₄}) = weyl_factor(c, 4)
conversion_factor(c::Conventions, ::Val{:φ₀}) = faraday_factor(c, 0)
conversion_factor(c::Conventions, ::Val{:φ₁}) = faraday_factor(c, 1)
conversion_factor(c::Conventions, ::Val{:φ₂}) = faraday_factor(c, 2)
conversion_factor(c::Conventions, ::Val{:σ}) = shear_factor(c)
conversion_factor(c::Conventions, ::Val{:λ}) = lambda_factor(c)
conversion_factor(c::Conventions, ::Val{:h}) = strain_factor(c)
conversion_factor(c::Conventions, ::Val{:News}) = news_factor(c)
