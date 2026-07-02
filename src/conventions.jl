"""
    Conventions{T}

Sign and scale factors pinning down a choice of conventions for asymptotic field quantities,
following the conventions appendix (sign/scale parameters ``s₀, s₁, s₂, s₃, λ, Θ, ζ`` there,
renamed here to mnemonic ``c`` parameters).  Each field is a ``c``-parameter; the defaults
reproduce the **SpEC** convention together with the **Newman–Penrose** ``ð``, which is the
convention `Scri.jl` uses natively.

| field  | meaning                                             | default |
|:-------|:----------------------------------------------------|:--------|
| `c_s`  | metric signature (`+1` ⟹ ``−+++``, `−1` ⟹ ``+−−−``) | `+1`    |
| `c_l`  | normalization/scale of the ``l`` tetrad leg         | `1`     |
| `c_m`  | spin phase ``Θ`` of ``m`` (radians)                 | `0`     |
| `c_R`  | sign of the Riemann tensor                          | `+1`    |
| `c_Ψ`  | sign of the Weyl components                         | `+1`    |
| `c_σ`  | sign of the shear ``σ``                             | `+1`    |
| `c_h`  | scale of the strain ``h``                           | `1`     |
| `c_φ`  | sign of the Maxwell/Faraday components              | `+1`    |
| `c_ð`  | coefficient of ``ð`` (`1` ⟹ Newman–Penrose)         | `1`     |

The defining relations (with ``g, η`` the ``−+++`` metrics, so explicit factors of `c_s`
carry the signature; see the conventions appendix):

  - ``l_a = -(c_l/\\sqrt2)(dt - dr)_a``,
  - ``m_a = (e^{i\\,c_m}/\\sqrt2)(dθ + i\\,dφ)_a``,
  - ``c_R\\,R^a{}_{bcd} = ∂Γ - ∂Γ + ΓΓ - ΓΓ``,
  - ``Ψ_4 = c_Ψ\\,C_{abcd} n^a \\bar m^b n^c \\bar m^d``,
  - ``σ = -c_σ\\,m^a m^b ∇_a l_b``,
  - ``h = c_h^{-1}[\\tfrac12(h_{θθ} - h_{φφ}) - i h_{θφ}]``, with ``h_{ab} = c_s(g_{ab}-η_{ab})``,
  - ``ð\\,{}_sf = -c_ð\\,(\\sinθ)^s(∂_θ + i\\cscθ\\,∂_φ)[(\\sinθ)^{-s}\\,{}_sf]``.

Construct with keywords (validated), or by name from the table of published conventions:

    Conventions()                 # SpEC + Newman–Penrose ð (the package default)
    Conventions(; c_s=-1, c_l=-√2) # a custom convention
    Conventions(:NP)              # Newman–Penrose (1968)
    Conventions(:MB; T=Double64)  # Moreschi/Boyle, at extended precision

`c_φ` and `c_ð` are not fixed by the (gravitational-wave) appendix table, so the named
presets leave them at the Newman–Penrose default `1`; set them explicitly for, e.g., a
Geroch–Held–Penrose ``ð`` (`c_ð = 1/√2`, since ``ð_NP = √2 ð_GHP`` when restricted to a unit
round 2-sphere).
"""
struct Conventions{T<:Real}
    c_s::Int
    c_l::T
    c_m::T
    c_R::Int
    c_Ψ::Int
    c_σ::Int
    c_h::T
    c_φ::Int
    c_ð::T
end

function Conventions(;
    c_s::Integer=1,
    c_l=1,
    c_m=0,
    c_R::Integer=1,
    c_Ψ::Integer=1,
    c_σ::Integer=1,
    c_h=1,
    c_φ::Integer=1,
    c_ð=1,
)
    for (name, c) ∈ (("c_s", c_s), ("c_R", c_R), ("c_Ψ", c_Ψ), ("c_σ", c_σ), ("c_φ", c_φ))
        c == 1 || c == -1 || throw(ArgumentError("$name is a sign and must be ±1; got $c"))
    end
    (iszero(c_l) || iszero(c_h) || iszero(c_ð)) &&
        throw(ArgumentError("c_l, c_h, and c_ð must be nonzero"))
    T = promote_type(
        typeof(float(c_l)), typeof(float(c_m)), typeof(float(c_h)), typeof(float(c_ð))
    )
    return Conventions{T}(c_s, T(c_l), T(c_m), c_R, c_Ψ, c_σ, T(c_h), c_φ, T(c_ð))
end

"""
    Conventions(name::Symbol; T=Float64)

Named presets from the conventions-appendix table.  Available: `:SpEC`, `:MB`
(Boyle/Lehner), `:NP` (Newman–Penrose), `:ADLK` (Ashtekar et al.), `:BR` (Bishop et al.),
`:C` (Chandrasekhar).  Parameters listed as N/A in the source (e.g. the strain scale `c_h` in
conventions that do not define a strain) are left at the default `1`.
"""
function Conventions(name::Symbol; T::Type{<:Real}=Float64)
    s2 = sqrt(T(2))
    if name === :SpEC
        return Conventions(; c_l=one(T))
    elseif name === :MB
        return Conventions(; c_s=-1, c_l=(-s2), c_h=T(2))
    elseif name === :NP
        return Conventions(; c_s=-1, c_Ψ=-1, c_l=(-s2), c_m=T(π))
    elseif name === :ADLK
        return Conventions(; c_σ=-1, c_l=(-s2))
    elseif name === :BR
        return Conventions(; c_Ψ=-1, c_l=one(T))
    elseif name === :C
        return Conventions(; c_s=-1, c_Ψ=-1, c_l=(-s2))
    else
        throw(
            ArgumentError(
                "Unknown convention preset :$name; " *
                "expected one of :SpEC, :MB, :NP, :ADLK, :BR, :C",
            ),
        )
    end
end

"""
    weyl_factor(c::Conventions, n::Integer)

Factor relating the Weyl component ``ψ_n`` in convention `c` to its value in the SpEC
convention: ``ψ_n^{[c]} = \\texttt{weyl\\_factor}(c, n)\\, ψ_n^{[SpEC]}``.  From the appendix,
``c_s c_Ψ c_R (c_l e^{i c_m})^{2-n}``.
"""
function weyl_factor(c::Conventions, n::Integer)
    return c.c_s * c.c_Ψ * c.c_R * (c.c_l * cis(c.c_m))^(2 - n)
end

"""
    strain_factor(c::Conventions)

Factor relating the strain ``h`` in convention `c` to its SpEC value:
``h^{[c]} = \\texttt{strain\\_factor}(c)\\, h^{[SpEC]} = c_s c_h^{-1} e^{-2i c_m}\\, h^{[SpEC]}``.
"""
strain_factor(c::Conventions) = c.c_s * cis(-2 * c.c_m) / c.c_h

"""
    convert_weyl(ψₙ, n, from::Conventions, to::Conventions)
    convert_strain(h, from::Conventions, to::Conventions)

Convert a Weyl component ``ψ_n`` (or the strain ``h``) between two conventions, via the SpEC
reference.  (Conversions for the shear ``σ`` and the Faraday components are more involved and
not yet implemented here.)
"""
function convert_weyl(ψ, n, from::Conventions, to::Conventions)
    return ψ * (weyl_factor(to, n) / weyl_factor(from, n))
end
function convert_strain(h, from::Conventions, to::Conventions)
    return h * (strain_factor(to) / strain_factor(from))
end
