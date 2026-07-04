# BMS group elements: a supertranslation followed by a Lorentz transformation.
#
# Conventions: see docs/src/80-details/30-bms_group.md

@doc raw"""
    BMS{T<:Real,A,I}

An element of the BMS group: a supertranslation `α` followed by a (proper orthochronous)
Lorentz transformation `Λ`.  We write elements as pairs

```math
(Λ, α) = (Λ, 0) ∘ (1, α),
```

where `Λ ∈ Spin⁺(3,1)` is represented by a [`Lorentz`](@extref
:jl:type:`Quaternionic.Lorentz`) rotor and the supertranslation is a real-valued function on
the sphere, stored as its spin-weight-0 spherical-harmonic mode weights.

The second type parameter `A ∈ {+1, -1}` is the **supertranslation-sign convention**: the
element acts as `t′ = κ(t - A α)`.  This says how the element's *own* stored `α` is read as
a time shift, so it is fixed at construction (default `+1`, set with the `c_α` keyword).
Because `A` factors out of the group law, an element `(Λ, α)` with `A = -1` is the same
transformation as `(Λ, -α)` with `A = +1`.

The third type parameter `I ∈ {+1, -1}` is the **null-infinity representation**: it records
which labeling of the celestial sphere the stored `α` modes refer to — the outgoing labeling
of ``ℐ⁺`` (`I = +1`, the default) or the antipodal past-light-cone labeling of ``ℐ⁻`` (`I =
-1`; set with the `ℐ` keyword).  The groups ``BMS⁺`` and ``BMS⁻`` are isomorphic under the
antipodal matching, so this is purely a choice of *representation*: the element `(Λ, α)`
with `I = -1` is the same abstract transformation as `(Λ, α∘A)` with `I = +1`, where `A` is
the antipodal map — `(-1)^ℓ` on each mode.  See the conversion constructor `BMS(g; c_α, ℐ)`
below, and [Representations of the supertranslation](@ref bms_representations).

Both conventions are pure bookkeeping.  Operations combining two elements
([`compose`](@ref), `==`, `isapprox`) accept any mix: the differences are accounted for
exactly, and wherever the two inputs disagree the result uses the default (`+1`) convention.
Data also live on a fixed null infinity, so [`DataComponents`](@ref) carries `ℐ` as a type
parameter; the [`transform!`](@ref) bridge re-represents the element to match the data.  See
[the conventions page](@ref scri_pm_conventions).

See [the BMS group documentation](@ref bms_group) for details about the group structure,
action, and conventions.

## Action on points of ℐ

Treating points of `ℐ` as pairs `(t, 𝐤)` of a time coordinate and a null ray, the element
`(Λ, α)` acts (passively) as

```math
t'(𝐤) = κ(𝐤) \left(t(𝐤) - c_α α(𝐤)\right),
\qquad
𝐤' = Λ 𝐤,
```

where the conformal factor is `κ(𝐤) = k⁰ / (Λk)⁰`, and `c_α = ±1` is the element's
supertranslation-sign convention (its `A` type parameter) — it affects only the action,
never the group structure.  See the action functor `(g::BMS)(t, n̂)` and
[`conformal_factor`](@ref).

## Storage

- `Λ::Lorentz{T}`: the Lorentz rotor.  Note that `±Λ` represent the same transformation
  (the double cover); `==`, `isapprox`, and `hash` account for this.
- `α::Vector{Complex{T}}`: mode weights of the supertranslation, of length `(ℓₘₐₓ+1)²`,
  ordered by increasing `ℓ ∈ 0:ℓₘₐₓ`, then increasing `m ∈ -ℓ:ℓ` — that is, the mode
  `(ℓ, m)` is stored at index `ℓ² + ℓ + m + 1`.  The constructor enforces the reality
  condition `α_{ℓ,-m} = (-1)^m ᾱ_{ℓ,m}` (via [`impose_reality`](@ref)), so that the
  function on the sphere is real-valued.

## Constructors

```julia
BMS(Λ, α)             # Lorentz rotor (or plain Rotor for a pure rotation) and modes
BMS{T}(; kwargs...)   # assemble from parts; see below
BMS(g; c_α, ℐ)        # convert an element between representations; see its own docstring
```

The keyword constructor accepts any combination of

- `supertranslation`: mode weights as above (any perfect-square length);
- `frame_rotation`: a `Rotor`;
- `boost_velocity`: a `QuatVec` or 3-vector with `|v⃗| < 1`;
- `time_translation`: a real `δt`;
- `space_translation`: a `QuatVec` or 3-vector `δx⃗`;
- `c_α`, `ℐ`: the conventions described above (both default `+1`);
- `ℓₘₐₓ`: the resolution of the stored `α` (defaults to the smallest that fits).

The Lorentz part is assembled as the *passive* version of "rotate, then boost": `Λ =
inv(Boost(v⃗) * Lorentz(R))`, so that `frame_rotation` and `boost_velocity` have exactly the
same meaning as the `R` and `v⃗` arguments of [`transform!`](@ref).  In particular `κ(𝐤) =
1/(γ(1 - v⃗⋅n̂))` for a pure boost.  The translation parts contribute `α(n̂) = δt - δx⃗⋅n̂`
— i.e., the supertranslation of the spacetime translation by the 4-vector `(δt, δx⃗)`.  If
`supertranslation` is also supplied, any `time_translation` and/or `space_translation`
content is *added* to it to produce the total supertranslation.

## Operations

Group composition via `*` (or `∘`, or [`compose`](@ref) with explicit band limits),
inverse via `inv`, identity via `one`.  Apply to a null ray with `g(t, n̂)`.  Accessors:
[`ℓₘₐₓ`](@ref), [`lorentz`](@ref), [`supertranslation`](@ref), [`frame_rotation`](@ref),
[`boost_velocity`](@ref), [`translation`](@ref), [`time_translation`](@ref),
[`space_translation`](@ref), [`proper_supertranslation`](@ref).
"""
struct BMS{T<:Real,A,I}
    Λ::Lorentz{T}
    α::Vector{Complex{T}}
    function BMS{T,A,I}(Λ::Lorentz{T}, α::AbstractVector) where {T<:Real,A,I}
        # `A = ±1` is the supertranslation-sign convention: the action shifts time by `t′ =
        # κ(t - A α)`.  `I = ±1` is the null-infinity representation: which labeling of the
        # sphere the stored α modes refer to.  Both are part of how *this element's own* α
        # is interpreted, so they are fixed at construction (validated even when asserts are
        # disabled).
        (A === 1 || A === -1) ||
            throw(ArgumentError("c_α (supertranslation sign) must be +1 or -1; got $A"))
        (I === 1 || I === -1) || throw(
            ArgumentError("ℐ (null-infinity representation) must be +1 or -1; got $I")
        )
        N = length(α)
        L = isqrt(N)
        # Note that N ≥ 1 is required because accessors rely on the ℓ=0 mode existing;
        # an explicit throw (rather than @assert) because this invariant must hold even
        # when asserts are disabled.
        if L^2 != N || N < 1
            throw(
                ArgumentError(
                    "Supertranslation must have perfect-square length ≥ 1; got length $N"
                ),
            )
        end
        return new{T,A,I}(Λ, impose_reality(convert(Vector{Complex{T}}, α), L - 1, 1))
    end
end

function BMS(
    Λ::Lorentz{T1}, α::AbstractVector{T2}; c_α::Int=1, ℐ::Int=1
) where {T1<:Real,T2<:Number}
    T = promote_type(T1, real(T2))
    return BMS{T,Int(c_α),Int(ℐ)}(Lorentz{T}(Λ), α)
end
function BMS(R::Rotor{T1}, α::AbstractVector; c_α::Int=1, ℐ::Int=1) where {T1<:Real}
    return BMS(Lorentz(R), α; c_α, ℐ)
end
function BMS{T}(R::Rotor, α::AbstractVector; c_α::Int=1, ℐ::Int=1) where {T<:Real}
    return BMS{T,Int(c_α),Int(ℐ)}(Lorentz{T}(Lorentz(R)), α)
end
function BMS{T}(Λ::Lorentz, α::AbstractVector; c_α::Int=1, ℐ::Int=1) where {T<:Real}
    return BMS{T,Int(c_α),Int(ℐ)}(Lorentz{T}(Λ), α)
end
# Fully parameterized forms that accept a non-`Lorentz{T}` rotor and convert; the
# two-parameter form defaults to the ℐ⁺ representation.
BMS{T,A}(x, α::AbstractVector) where {T<:Real,A} = BMS{T,A,1}(x, α)
function BMS{T,A,I}(R::Rotor, α::AbstractVector) where {T<:Real,A,I}
    return BMS{T,A,I}(Lorentz{T}(Lorentz(R)), α)
end
function BMS{T,A,I}(Λ::Lorentz, α::AbstractVector) where {T<:Real,A,I}
    return BMS{T,A,I}(Lorentz{T}(Λ), α)
end

@doc raw"""
    BMS(g::BMS; c_α=+1, ℐ=+1)
    BMS{T,A,I}(g::BMS)

The same abstract BMS transformation as `g`, re-expressed with the requested conventions
(and, in the fully parameterized form, float type).

Changing `c_α` negates the stored modes (`t′ = κ(t - c_αα)` is unchanged under `(c_α, α) ↦
(-c_α, -α)`).  Changing `ℐ` relabels the sphere by the antipodal map `A` — the antipodal
matching that identifies ``BMS⁺`` and ``BMS⁻`` — so the modes transform as `α ↦ α∘A`, which
is `(-1)^ℓ` on each `ℓ` block, since `Yₗₘ(-n̂) = (-1)^ℓ Yₗₘ(n̂)`.  The Lorentz part is
unchanged in both cases.  Both operations are exact (pure sign flips) and involutive.  See
[Representations of the supertranslation](@ref bms_representations).

Note the distinction from constructing an element with non-default conventions from *raw
parts*: `BMS(Λ, α; ℐ=-1)` attaches the ℐ⁻ meaning to the raw modes `α`, whereas `BMS(BMS(Λ,
α); ℐ=-1)` converts the ℐ⁺ element `(Λ, α)` into its matched ℐ⁻ form — a different element
unless `α` has no odd-`ℓ` content.
"""
function BMS(g::BMS{T}; c_α::Int=1, ℐ::Int=1) where {T<:Real}
    return BMS{T,Int(c_α),Int(ℐ)}(g)
end

BMS{T,A,I}(g::BMS{T,A,I}) where {T<:Real,A,I} = g
function BMS{T,A,I}(g::BMS{S,A₀,I₀}) where {T<:Real,S<:Real,A,I,A₀,I₀}
    α = copy(g.α)
    if I != I₀
        for ℓ ∈ 1:2:ℓₘₐₓ(g)  # α ↦ α∘A: odd-ℓ blocks flip sign
            α[(ℓ ^ 2 + 1):((ℓ + 1) ^ 2)] .*= -1
        end
    end
    if A != A₀
        α .*= -1
    end
    return BMS{T,A,I}(Lorentz{T}(g.Λ), α)
end
# Converting with only the float type specified preserves both conventions.
BMS{T}(g::BMS{S,A,I}) where {T<:Real,S<:Real,A,I} = BMS{T,A,I}(g)
Base.convert(::Type{BMS{T}}, g::BMS) where {T<:Real} = BMS{T}(g)
# `convert` may change the representation: that is value-preserving, since `==` compares
# elements as abstract transformations.
Base.convert(::Type{BMS{T,A,I}}, g::BMS) where {T<:Real,A,I} = BMS{T,A,I}(g)
function Base.promote_rule(
    ::Type{BMS{T1,A,I}}, ::Type{BMS{T2,A,I}}
) where {T1<:Real,T2<:Real,A,I}
    # Only elements with the same conventions promote together; convert conventions
    # explicitly with `BMS(g; c_α, ℐ)`.
    return BMS{promote_type(T1, T2),A,I}
end
Quaternionic.basetype(::Type{BMS{T,A,I}}) where {T<:Real,A,I} = T
Quaternionic.basetype(::BMS{T}) where {T<:Real} = T

@doc raw"""
    translation_modes(T, δt, δx⃗, ℓₘₐₓ=1)

Mode weights (in the ordering used by [`BMS`](@ref)) of the supertranslation corresponding
to the spacetime translation by the 4-vector `(δt, δx⃗)`:

```math
α(n̂) = δt - δx⃗ ⋅ n̂.
```

In terms of (orthonormal, Condon-Shortley) spherical harmonics, the nonzero modes are

```math
α_{0,0} = 2\sqrt{π}\, δt,
\qquad
α_{1,0} = -2\sqrt{π/3}\, δz,
\qquad
α_{1,±1} = ±\sqrt{2π/3}\, (δx ∓ i\, δy).
```

The output is zero-padded up to `ℓₘₐₓ` (which must be at least 1 whenever `δx⃗ ≠ 0`).
"""
function translation_modes(::Type{T}, δt::Real, δx⃗::QuatVec, ℓₘₐₓ::Int=1) where {T<:Real}
    α = zeros(Complex{T}, (ℓₘₐₓ + 1)^2)
    α[1] = 2 * √T(π) * δt
    δx, δy, δz = vec(δx⃗)
    if ℓₘₐₓ ≥ 1
        α[2] = -√(2T(π) / 3) * (δx + im * δy)  # (ℓ,m) = (1,-1)
        α[3] = -2 * √(T(π) / 3) * δz           # (ℓ,m) = (1, 0)
        α[4] = +√(2T(π) / 3) * (δx - im * δy)  # (ℓ,m) = (1,+1)
    else
        @assert iszero(δx) && iszero(δy) && iszero(δz) "" *
            "ℓₘₐₓ must be at least 1 to represent a space translation"
    end
    return α
end

function BMS{T}(;
    supertranslation=nothing,
    frame_rotation=nothing,
    boost_velocity=nothing,
    time_translation=nothing,
    space_translation=nothing,
    c_α::Int=1,
    ℐ::Int=1,
    ℓₘₐₓ::Union{Nothing,Int}=nothing,
) where {T<:Real}
    R = isnothing(frame_rotation) ? one(Rotor{T}) : Rotor{T}(frame_rotation)
    v⃗ = as_quatvec(T, boost_velocity)
    if absvec(v⃗) ≥ 1
        throw(DomainError(absvec(v⃗), "boost speed must satisfy |v⃗| < 1 (c = 1)"))
    end
    # Passive version of "rotate by R, then boost by v⃗" — the same meaning these
    # parameters have in `transform!`.  (Rotor products renormalize, perturbing the last
    # ulps, so we skip the product when either factor is missing.)
    Λ = if isnothing(frame_rotation) && isnothing(boost_velocity)
        one(Lorentz{T})
    elseif isnothing(frame_rotation)
        inv(Boost(v⃗))
    elseif isnothing(boost_velocity)
        inv(Lorentz(R))
    else
        inv(Boost(v⃗) * Lorentz(R))
    end
    δt = isnothing(time_translation) ? zero(T) : T(time_translation)
    δx⃗ = as_quatvec(T, space_translation)
    ℓᵅ = isnothing(supertranslation) ? 0 : isqrt(length(supertranslation)) - 1
    if !isnothing(supertranslation)
        @assert (ℓᵅ + 1)^2 == length(supertranslation) "" *
            "Supertranslation must have perfect-square length; " *
            "got length $(length(supertranslation))"
    end
    translating = !isnothing(time_translation) || !isnothing(space_translation)
    L = something(ℓₘₐₓ, max(ℓᵅ, translating ? 1 : 0))
    if L < ℓᵅ
        throw(ArgumentError("ℓₘₐₓ=$L is too small for the given supertranslation"))
    end
    α = translation_modes(T, δt, δx⃗, L)
    if !isnothing(supertranslation)
        α[1:length(supertranslation)] .+= supertranslation
    end
    return BMS{T,Int(c_α),Int(ℐ)}(Λ, α)
end

function BMS(; kwargs...)
    # Infer the float type from the given parts; all-`nothing` (and bare-integer) kwargs
    # promote to Bool, which floats to Float64.
    T = float(reduce(promote_type, map(kwarg_basetype, values(values(kwargs))); init=Bool))
    return BMS{T}(; kwargs...)
end

# Helpers for the keyword constructors
kwarg_basetype(::Nothing) = Bool
kwarg_basetype(x::Real) = typeof(x)
kwarg_basetype(q::Rotor{T}) where {T} = real(T)
kwarg_basetype(q::QuatVec{T}) where {T} = T
kwarg_basetype(v::AbstractVector) = real(eltype(v))
kwarg_basetype(::Int) = Bool  # ℓₘₐₓ should not influence the float type
as_quatvec(::Type{T}, ::Nothing) where {T} = QuatVec{T}(0, 0, 0)
as_quatvec(::Type{T}, v⃗::QuatVec) where {T} = QuatVec{T}(v⃗)
function as_quatvec(::Type{T}, v⃗::AbstractVector) where {T}
    @assert length(v⃗) == 3 "Expected a 3-vector; got length $(length(v⃗))"
    return QuatVec{T}(v⃗[begin], v⃗[begin + 1], v⃗[begin + 2])
end

###
### Accessors
###

"""
    ℓₘₐₓ(g::BMS)

Largest `ℓ` of the spherical-harmonic modes stored in the supertranslation of `g`.
"""
ℓₘₐₓ(g::BMS) = isqrt(length(g.α)) - 1
ell_max(g::BMS) = ℓₘₐₓ(g)

"""
    c_α(g::BMS)

The supertranslation-sign convention of `g`: `+1` or `-1`, where the element acts as
`t′ = κ(t - c_α α)`.  See [`BMS`](@ref).
"""
c_α(::BMS{T,A}) where {T,A} = A

"""
    ℐ(g::BMS)

The null-infinity representation of `g`: `+1` if the stored supertranslation modes refer
to the outgoing labeling of ℐ⁺, `-1` for the antipodal past-light-cone labeling of ℐ⁻.
See [`BMS`](@ref), including its conversion constructor `BMS(g; c_α, ℐ)`.
"""
ℐ(::BMS{T,A,I}) where {T,A,I} = I

"""
    lorentz(g::BMS)

The Lorentz part `Λ` of `g = (Λ, α)`.
"""
lorentz(g::BMS) = g.Λ

"""
    supertranslation(g::BMS)

The supertranslation mode weights `α` of `g = (Λ, α)`.  This is the stored vector itself;
treat it as read-only.
"""
supertranslation(g::BMS) = g.α

"""
    frame_rotation(g::BMS)

The rotation `R` such that the Lorentz part of `g` represents the (passive) frame change
"rotate by `R`, then boost by [`boost_velocity`](@ref)`(g)`" — that is,
`lorentz(g) == inv(Boost(v⃗) * Lorentz(R))`, matching the meaning of the `R` argument of
[`transform!`](@ref).
"""
frame_rotation(g::BMS) = last(Quaternionic.vR(inv(g.Λ)))

"""
    boost_velocity(g::BMS)

The boost velocity `v⃗` of the frame change represented by the Lorentz part of `g`; see
[`frame_rotation`](@ref).
"""
boost_velocity(g::BMS) = first(Quaternionic.vR(inv(g.Λ)))

@doc raw"""
    translation(g::BMS) -> (δt, δx⃗)

Interpret the `ℓ ≤ 1` modes of the supertranslation of `g` as a spacetime translation,
returning the time translation `δt` and space translation `δx⃗`.  This inverts
[`translation_modes`](@ref); any `ℓ ≥ 2` ("proper supertranslation") content is ignored.
"""
function translation(g::BMS{T}) where {T<:Real}
    α = g.α
    δt = real(α[1]) / (2 * √T(π))
    if length(α) ≥ 4
        c = √(3 / (2T(π)))
        δx = +c * real(α[4])
        δy = -c * imag(α[4])
        δz = -real(α[3]) * √(3 / T(π)) / 2
        return (δt, QuatVec{T}(δx, δy, δz))
    else
        return (δt, QuatVec{T}(0, 0, 0))
    end
end

"""
    time_translation(g::BMS)

The time-translation part `δt` of the supertranslation of `g`; see [`translation`](@ref).
"""
time_translation(g::BMS) = first(translation(g))

"""
    space_translation(g::BMS)

The space-translation part `δx⃗` of the supertranslation of `g`; see [`translation`](@ref).
"""
space_translation(g::BMS) = last(translation(g))

"""
    proper_supertranslation(g::BMS)

A copy of the supertranslation modes of `g` with the `ℓ ≤ 1` (spacetime-translation) modes
zeroed out, leaving only the "proper" supertranslation content.
"""
function proper_supertranslation(g::BMS{T}) where {T<:Real}
    α = copy(g.α)
    α[1:min(4, length(α))] .= 0
    return α
end

###
### Equality, hashing, approximation, display
###

function padded_modes_equal(α₁, α₂)
    N = max(length(α₁), length(α₂))
    z = zero(eltype(α₁))
    for i ∈ 1:N
        a = i ≤ length(α₁) ? α₁[i] : z
        b = i ≤ length(α₂) ? α₂[i] : zero(eltype(α₂))
        a == b || return false
    end
    return true
end

"""
    g₁ == g₂

Mathematical equality of BMS elements — do they denote the same abstract transformation?
Convention differences (`A`, `I`) are accounted for exactly by converting both elements
to the default representation (the `BMS(g)` conversion constructor); then the Lorentz parts must agree
*up to overall sign* (since `±Λ` represent the same transformation), and the
supertranslation modes must agree after zero-padding both to a common `ℓₘₐₓ` (trailing
zero modes are physically vacuous).
"""
function Base.:(==)(g₁::BMS, g₂::BMS)
    c₁ = BMS(g₁)
    c₂ = BMS(g₂)
    (c₁.Λ == c₂.Λ || c₁.Λ == -c₂.Λ) || return false
    return padded_modes_equal(c₁.α, c₂.α)
end

function Base.hash(g::BMS, h::UInt)
    # Hash the default representation so that `==`-equal elements hash equally,
    # regardless of the conventions they are expressed in.
    c = BMS(g)
    # Strip trailing all-zero ℓ-blocks so padding does not affect the hash.
    L = isqrt(length(c.α))
    while L > 1 && all(iszero, @view c.α[((L - 1) ^ 2 + 1):(L ^ 2)])
        L -= 1
    end
    # Canonicalize the sign of Λ: flip so that the first nonzero entry of the flattened
    # components is positive.
    s = one(basetype(c))
    for x₀ ∈ components(c.Λ), x ∈ (real(x₀), imag(x₀))
        if !iszero(x)
            s = copysign(s, x)
            break
        end
    end
    h = hash(:BMS, h)
    h = hash(components(s < 0 ? -c.Λ : c.Λ), h)
    for i ∈ 1:(L ^ 2)
        h = hash(c.α[i], h)
    end
    return h
end

"""
    isapprox(g₁::BMS, g₂::BMS; kwargs...)

Approximate equality of BMS elements: both are converted to the default representation
(the `BMS(g)` conversion constructor; exact — pure sign flips), then the Lorentz parts are compared
up to overall sign and the supertranslation modes after zero-padding to a common length.
Keyword arguments are passed through to the underlying vector comparisons; as usual,
supply `atol` when comparing against (near-)zero quantities.
"""
function Base.isapprox(g₁′::BMS, g₂′::BMS; kwargs...)
    g₁ = BMS(g₁′)
    g₂ = BMS(g₂′)
    c₁ = collect(components(g₁.Λ))
    c₂ = collect(components(g₂.Λ))
    (isapprox(c₁, c₂; kwargs...) || isapprox(c₁, -c₂; kwargs...)) || return false
    N = max(length(g₁.α), length(g₂.α))
    T = promote_type(basetype(g₁), basetype(g₂))
    α₁ = zeros(Complex{T}, N)
    α₁[1:length(g₁.α)] = g₁.α
    α₂ = zeros(Complex{T}, N)
    α₂[1:length(g₂.α)] = g₂.α
    return isapprox(α₁, α₂; kwargs...)
end

function Base.show(io::IO, g::BMS{T,A,I}) where {T,A,I}
    return print(io, "BMS{", T, ",", A, ",", I, "}(", g.Λ, ", ", g.α, ")")
end
function Base.show(io::IO, ::MIME"text/plain", g::BMS{T,A,I}) where {T,A,I}
    signα = A == 1 ? "+" : "-"
    signᴵ = I == 1 ? "+" : "-"
    println(
        io, "BMS{", T, "} (c_α = ", signα, "1, ℐ = ", signᴵ, "1) with ℓₘₐₓ = ", ℓₘₐₓ(g), ":"
    )
    println(io, "  Λ = ", g.Λ)
    return print(io, "  α = ", g.α)
end

###
### Inline tests: construction, equality, translations
###

@testmodule BMSTestSetup begin
    import Random
    import LinearAlgebra: lu, dot, norm
    import Quaternionic
    import Quaternionic:
        QuatVec,
        Rotor,
        Lorentz,
        Boost,
        rotor,
        components,
        absvec,
        𝐤,
        from_spherical_coordinates
    import SphericalFunctions: ₛ𝐘, golden_ratio_spiral_rotors, D_matrices, WignerDindex
    using Scri: BMS, impose_reality
    using DoubleFloats: Double64

    const FloatTypes = (Float32, Float64, Double64, BigFloat)

    # Conditioning-scaled tolerance for band-limited synthesis/analysis round-trips on
    # the golden-ratio grid at working bandwidth ℓʷ.
    tol(::Type{T}, ℓʷ; c=100) where {T} = c * (ℓʷ + 1)^2 * eps(real(T))

    function random_direction(rng, ::Type{T}) where {T}
        x, y, z = randn(rng, T), randn(rng, T), randn(rng, T)
        n = hypot(x, y, z)
        return QuatVec(x / n, y / n, z / n)
    end

    random_rotation(rng, ::Type{T}) where {T} = randn(rng, Rotor{T})

    function random_boost(rng, ::Type{T}; βmax=1//2) where {T}
        return Boost(T(βmax) * rand(rng, T) * random_direction(rng, T))
    end

    function random_lorentz(rng, ::Type{T}; βmax=1//2) where {T}
        return Lorentz(random_rotation(rng, T)) * random_boost(rng, T; βmax)
    end

    function random_α(rng, ::Type{T}, ℓₘₐₓ) where {T}
        return impose_reality(randn(rng, Complex{T}, (ℓₘₐₓ + 1)^2), ℓₘₐₓ, 1)
    end

    function random_bms(rng, ::Type{T}; ℓₘₐₓ=2, βmax=1//2) where {T}
        return BMS(random_lorentz(rng, T; βmax), random_α(rng, T, ℓₘₐₓ))
    end

    ###
    ### Independent oracles — deliberately re-derived from raw Quaternionic geometry and
    ### SphericalFunctions primitives, NOT from the functions under test.
    ###

    # κ and the aberrated direction, straight from the 4-vector action on the section
    # 𝐤 = (1, ℐ n̂); ℐ = +1 (ℐ⁺) or -1 (ℐ⁻).  The trailing ℐ on the spatial part returns
    # the direction in the same n̂ labeling.
    function ray_map(Λ::Lorentz{T1}, n̂::QuatVec{T2}; ℐ=1) where {T1,T2}
        T = promote_type(T1, T2)
        nˣ, nʸ, nᶻ = vec(n̂)
        k′ = Λ(T[1, ℐ * nˣ, ℐ * nʸ, ℐ * nᶻ])
        return (1 / k′[1], QuatVec(ℐ * k′[2], ℐ * k′[3], ℐ * k′[4]) / k′[1])
    end

    # A rotor whose action takes 𝐤 to n̂ (any will do for spin weight 0).
    function rotor_pointing(n̂::QuatVec)
        _, x, y, z = components(n̂)
        return from_spherical_coordinates(atan(hypot(x, y), z), atan(y, x))
    end

    # Evaluate spin-0 mode weights at the given rotors / at a single direction.
    function α_eval(α::Vector{Complex{T}}, Rs::AbstractVector{<:Rotor}) where {T}
        return real.(ₛ𝐘(0, isqrt(length(α)) - 1, T, Rs) * α)
    end
    α_value(α, n̂::QuatVec) = α_eval(α, [rotor_pointing(n̂)])[1]

    # Reference action of (Λ, α) on a null ray, from raw geometry only.
    function act_ref(Λ, α, t, n̂; c_α=1, ℐ=1)
        κ, n̂′ = ray_map(Λ, n̂; ℐ)
        return (κ * (t - c_α * α_value(α, n̂)), n̂′)
    end

    # Pointwise supertranslation of the spacetime translation by (δt, δx⃗).
    α_translation(δt, δx⃗) = n̂ -> δt - dot(vec(δx⃗), vec(n̂))

    # Rotate spin-0 mode weights: if f′(n̂) = f(Q⁻¹ n̂), then (pinned empirically against
    # pointwise evaluation; see the "rotation conjugation" test items)
    #     f′_{ℓ,m′} = Σₘ conj(D^ℓ_{m′,m}(Q)) f_{ℓ,m}.
    function rotate_modes(α::Vector{Complex{T}}, Q::Rotor, ℓₘₐₓ) where {T}
        D = D_matrices(Q, ℓₘₐₓ)
        out = zeros(Complex{T}, (ℓₘₐₓ + 1)^2)
        for ℓ ∈ 0:ℓₘₐₓ, m′ ∈ (-ℓ):ℓ
            s = zero(Complex{T})
            for m ∈ (-ℓ):ℓ
                s += conj(D[WignerDindex(ℓ, m′, m)]) * α[ℓ ^ 2 + ℓ + m + 1]
            end
            out[ℓ ^ 2 + ℓ + m′ + 1] = s
        end
        return out
    end
end

@testitem "BMS: construction, reality, and padding" tags = [:unit, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_lorentz
    using Quaternionic: Lorentz, basetype
    using Scri: impose_reality, ℓₘₐₓ, lorentz, supertranslation

    rng = Random.Xoshiro(314)
    for T ∈ FloatTypes
        Λ = random_lorentz(rng, T)
        αraw = randn(rng, Complex{T}, 9)
        g = BMS(Λ, αraw)
        # The constructor imposes the reality condition, exactly as impose_reality does.
        @test g.α == impose_reality(αraw, 2, 1)
        @test ℓₘₐₓ(g) == 2
        @test lorentz(g) === g.Λ
        @test supertranslation(g) === g.α
        # Idempotence: re-wrapping the (already real) modes is exact.
        @test BMS(Λ, g.α).α == g.α
        # Real input modes are accepted.  Note that the reality condition averages the
        # input's (1,±1) modes to (1∓1)/2 = 0 here, leaving the m=0 modes alone.
        @test BMS(Λ, ones(T, 4)).α == Complex{T}[1, 0, 1, 0]
    end

    # Non-square and empty lengths are rejected.  The empty case matters separately:
    # isqrt(0)² == 0 would slip through the perfect-square check alone, but accessors
    # like `translation` rely on the ℓ=0 mode always being present.
    for N ∈ (0, 2, 3, 5, 8)
        @test_throws ArgumentError BMS(one(Lorentz{Float64}), zeros(ComplexF64, N))
    end

    # Type promotion across the two fields.
    g = BMS(one(Lorentz{Float32}), zeros(ComplexF64, 4))
    @test g isa BMS{Float64}
    @test BMS{Float32}(g) isa BMS{Float32}
    @test basetype(g) === Float64

    # The identity element.
    for T ∈ FloatTypes
        e = one(BMS{T})
        @test e isa BMS{T}
        @test isone(e)
        @test ℓₘₐₓ(e) == 0
        @test all(iszero, e.α)
    end
end

@testitem "BMS: equality, isapprox, hash" tags = [:unit, :fast] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: random_lorentz, random_α

    rng = Random.Xoshiro(159)
    Λ = random_lorentz(rng, Float64)
    α = random_α(rng, Float64, 2)

    g = BMS(Λ, α)
    # Double cover: ±Λ represent the same transformation.
    @test BMS(-Λ, α) == g
    @test hash(BMS(-Λ, α)) == hash(g)
    @test BMS(-Λ, α) ≈ g

    # Padding with zero modes is physically vacuous.
    gpad = BMS(Λ, [α; zeros(ComplexF64, 16 - 9)])
    @test gpad == g
    @test hash(gpad) == hash(g)
    @test gpad ≈ g

    # Genuinely different elements differ.
    @test BMS(Λ, 2α) != g
    @test !(BMS(Λ, 2α) ≈ g)
    Λ′ = random_lorentz(rng, Float64)
    @test BMS(Λ′, α) != g

    # Mixed-precision equality (values exactly representable in Float32; we set the
    # modes directly because `time_translation` would introduce a √π that rounds
    # differently in the two precisions).
    g32 = BMS{Float32}(; supertranslation=ComplexF32[1 // 2])
    g64 = BMS{Float64}(; supertranslation=ComplexF64[1 // 2])
    @test g32 == g64
    @test hash(g32) == hash(g64)

    # isapprox distinguishes near-but-not-equal elements through tolerances.
    gϵ = BMS(Λ, α .+ 1e-14)
    @test gϵ ≈ g
    @test !(isapprox(gϵ, g; atol=1e-20, rtol=0))
end

@testitem "BMS: kwarg constructor assembles parts" tags = [:unit, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_rotation, random_direction
    using Quaternionic: Quaternionic, QuatVec, Rotor, Lorentz, Boost, components
    using Scri: boost_velocity, frame_rotation, impose_reality, ℓₘₐₓ, translation_modes

    rng = Random.Xoshiro(265)
    for T ∈ (Float32, Float64, BigFloat)
        R = random_rotation(rng, T)
        v⃗ = T(3//10) * random_direction(rng, T)

        # Identity when no parts are given.
        @test isone(BMS{T}())

        # Pure rotation: Λ = inv(Lorentz(R)) = Lorentz(conj(R)).
        g = BMS{T}(; frame_rotation=R)
        @test g.Λ == Lorentz(conj(R))
        @test all(iszero, g.α)

        # Pure boost: Λ = inv(Boost(v⃗)) = Boost(-v⃗).
        g = BMS{T}(; boost_velocity=v⃗)
        @test g.Λ == Boost(-v⃗)

        # Combined: passive version of rotate-then-boost.
        g = BMS{T}(; frame_rotation=R, boost_velocity=v⃗)
        @test g.Λ == inv(Boost(v⃗) * Lorentz(R))

        # Round-trip through the accessors.
        v⃗′ = boost_velocity(g)
        R′ = frame_rotation(g)
        @test maximum(abs, components(v⃗′ - v⃗)) < 20eps(T)
        @test min(maximum(abs, components(R′ - R)), maximum(abs, components(R′ + R))) <
            20eps(T)

        # Supertranslation parts add exactly.
        αin = impose_reality(randn(rng, Complex{T}, 9), 2, 1)
        g = BMS{T}(; supertranslation=αin, time_translation=T(2), ℓₘₐₓ=3)
        @test ℓₘₐₓ(g) == 3
        expected = zeros(Complex{T}, 16)
        expected[1:9] .= αin
        expected .+= translation_modes(T, T(2), QuatVec{T}(0, 0, 0), 3)
        @test g.α == expected
    end

    # Speeds at or above c are rejected.
    @test_throws DomainError BMS{Float64}(; boost_velocity=[0.0, 0.0, 1.0])
    @test_throws DomainError BMS{Float64}(; boost_velocity=[0.8, 0.8, 0.0])
    # A too-small ℓₘₐₓ is rejected.
    @test_throws ArgumentError BMS{Float64}(; supertranslation=zeros(9), ℓₘₐₓ=1)
    # Zero velocity gives the identity Lorentz element exactly.
    @test isone(BMS{Float64}(; boost_velocity=QuatVec(0.0, 0.0, 0.0)))

    # The un-parameterized constructor infers the float type.
    @test BMS(; time_translation=1.5) isa BMS{Float64}
    @test BMS(; time_translation=1.5f0) isa BMS{Float32}
    @test BMS() isa BMS{Float64}
    @test BMS(; time_translation=2) isa BMS{Float64}
end

@testitem "BMS: translation modes round-trip and pointwise" tags = [
    :unit, :fast, :validation
] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_direction, α_value, α_translation
    using Quaternionic: QuatVec, components, absvec
    using Scri: translation, time_translation, proper_supertranslation

    rng = Random.Xoshiro(358)
    for T ∈ FloatTypes
        ϵ = 50eps(T)
        for _ ∈ 1:5
            δt = randn(rng, T)
            δx⃗ = QuatVec(randn(rng, T), randn(rng, T), randn(rng, T))
            g = BMS{T}(; time_translation=δt, space_translation=δx⃗)

            # Round-trip through the accessors is exact up to roundoff.
            δt′, δx⃗′ = translation(g)
            @test abs(δt′ - δt) < ϵ
            @test maximum(abs, components(δx⃗′ - δx⃗)) < ϵ
            @test abs(time_translation(g) - δt) < ϵ
            @test all(iszero, proper_supertranslation(g))

            # Pointwise: α(n̂) = δt - δx⃗⋅n̂.  This pins the spherical-harmonic
            # conventions of `translation_modes` once and for all.
            scale = max(abs(δt), absvec(δx⃗))
            for _ ∈ 1:5
                n̂ = random_direction(rng, T)
                @test abs(α_value(g.α, n̂) - α_translation(δt, δx⃗)(n̂)) < 200eps(T) * scale
            end
        end
    end
end

###
### Geometry primitives and the action on null rays
###

@doc raw"""
    transform_ray(Λ, n̂; ℐ=+1) -> (κ, n̂′)

Map the null ray `𝐤 = (1, ℐ n̂)` through the Lorentz transformation `Λ`, returning the
conformal factor and the new direction:

```math
κ = \frac{k⁰}{(Λk)⁰},
\qquad
n̂' = ℐ\,\frac{\overrightarrow{Λk}}{(Λk)⁰}.
```

The null-infinity sign `ℐ = ±1` selects the section `𝐤 = (1, ℐ n̂)` — outgoing
(future-cone) labeling at ``ℐ⁺`` (`ℐ = +1`) or the antipodal past-light-cone labeling at
``ℐ⁻`` (`ℐ = -1`); see [the conventions page](@ref scri_pm_conventions).  For a pure boost
this gives `κ = 1/(γ(1 - ℐ v⃗⋅n̂))`, and the trailing `ℐ` on `n̂′` undoes the section sign
so that the returned direction is in the same `n̂` labeling (rotations map `n̂ ↦ Rn̂`
regardless of `ℐ`).  The `ℐ = -1` map agrees with [`aberration`](@ref)`(…, -1)`.

This is the single place where the conformal-factor convention is defined; everything
else (the action, [`compose`](@ref), [`conformal_factor`](@ref)) uses it.  The direction
`n̂` must be a unit vector.
"""
function transform_ray(Λ::Lorentz{T1}, n̂::QuatVec{T2}; ℐ::Int=1) where {T1<:Real,T2<:Real}
    T = promote_type(T1, T2)
    nˣ, nʸ, nᶻ = vec(n̂)
    k′ = Λ(T[1, ℐ * nˣ, ℐ * nʸ, ℐ * nᶻ])
    κ = inv(k′[1])
    return (κ, QuatVec{T}(ℐ * k′[2] * κ, ℐ * k′[3] * κ, ℐ * k′[4] * κ))
end

"""
    conformal_factor(Λ, n̂; ℐ=+1)
    conformal_factor(g::BMS, n̂)

The conformal factor `κ(𝐤) = k⁰/(Λk)⁰` of the Lorentz transformation at the null ray
`𝐤 = (1, ℐ n̂)`; see [`transform_ray`](@ref).  For a raw `Lorentz` rotor — which is
representation-neutral — the null-infinity sign `ℐ` is supplied per call (default `+1`
for ``ℐ⁺``); for a [`BMS`](@ref) element it is the element's own `I`, so `n̂` is read in
the element's labeling.  Use the conversion constructor `BMS(g; ℐ=…)` to realize the
element on the other ``ℐ``.
"""
conformal_factor(Λ::Lorentz, n̂::QuatVec; ℐ::Int=1) = first(transform_ray(Λ, n̂; ℐ))
function conformal_factor(g::BMS{T,A,I}, n̂::QuatVec) where {T,A,I}
    return conformal_factor(g.Λ, n̂; ℐ=I)
end

"""
    rotor_from_direction(n̂)

A rotor `R` such that `R(𝐤) = n̂`, where `𝐤` is the unit vector in the `z` direction.
Any such rotor would do for evaluating spin-weight-0 functions; this one is built from
the spherical coordinates of `n̂`, using pole-safe two-argument `atan` forms.
"""
function rotor_from_direction(n̂::QuatVec)
    x, y, z = vec(n̂)
    return from_spherical_coordinates(atan(hypot(x, y), z), atan(y, x))
end

"""
    supertranslation_values(g::BMS, Rs)

Evaluate the supertranslation of `g` at the directions given by the rotors `Rs` (each
direction being `R(𝐤)`).  Returns a real vector.
"""
function supertranslation_values(g::BMS{T}, Rs::AbstractVector{<:Rotor}) where {T<:Real}
    Tp = promote_type(T, basetype(eltype(Rs)))
    return real.(ₛ𝐘(0, ℓₘₐₓ(g), Tp, Rs) * g.α)
end

@doc raw"""
    (g::BMS)(t, n̂) -> (t′, n̂′)

Act with `g = (Λ, α)` on the point of `ℐ` labeled by time `t` and (unit) direction `n̂`:

```math
t' = κ(𝐤)\,\left(t - c_α\,α(𝐤)\right),
\qquad
𝐤' = Λ𝐤,
```

with `κ` as in [`transform_ray`](@ref).  Both conventions are the element's own type
parameters, fixed at construction: the supertranslation sign `c_α` says how `α` enters the
time shift, and the null-infinity representation `ℐ` says which ``ℐ`` — and which
labeling of it — the coordinates `(t, n̂)` refer to.  To act on the other ``ℐ``, convert
the element first with the conversion constructor `BMS(g; ℐ=…)`.
"""
function (g::BMS{T,A,I})(t::Real, n̂::QuatVec) where {T<:Real,A,I}
    κ, n̂′ = transform_ray(g.Λ, n̂; ℐ=I)
    αₙ = supertranslation_values(g, [rotor_from_direction(n̂)])[1]
    return (κ * (t - A * αₙ), n̂′)
end

@testitem "BMS: conformal factor properties" tags = [:unit, :fast, :validation] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_direction, random_rotation, random_lorentz
    using Quaternionic: QuatVec, Lorentz, Boost, 𝐤
    using Scri: conformal_factor

    rng = Random.Xoshiro(979)
    for T ∈ FloatTypes
        n̂ = random_direction(rng, T)
        # Identity has κ ≡ 1, exactly.
        @test conformal_factor(one(Lorentz{T}), n̂) == 1
        # Rotations have κ ≡ 1 up to roundoff.
        for _ ∈ 1:5
            R = random_rotation(rng, T)
            @test abs(conformal_factor(Lorentz(R), random_direction(rng, T)) - 1) < 20eps(T)
        end
        # κ > 0 for any (orthochronous) Lorentz transformation.
        for _ ∈ 1:5
            Λ = random_lorentz(rng, T; βmax=9//10)
            @test conformal_factor(Λ, random_direction(rng, T)) > 0
        end
        # Boost along +z: (Λk)⁰ = e^{±η} at n̂ = ±ẑ, so κ = e^{∓η}.
        η = T(7//10)
        B = Boost(η, QuatVec{T}(0, 0, 1))
        @test abs(conformal_factor(B, QuatVec{T}(0, 0, 1)) - exp(-η)) < 20eps(T)
        @test abs(conformal_factor(B, QuatVec{T}(0, 0, -1)) - exp(η)) < 20eps(T)
        # And against the hand-rolled Doppler form 1/(γ(1 ∓ β)) for the active boost.
        β = tanh(η)
        γ = 1 / √(1 - β^2)
        @test abs(conformal_factor(B, QuatVec{T}(0, 0, 1)) - 1 / (γ * (1 + β))) < 20eps(T)
    end
end

@testitem "BMS: ℐ⁻ conventions via the ℐ representation" tags = [:unit, :fast, :validation] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup:
        FloatTypes, random_direction, random_rotation, random_bms, act_ref, tol
    using Quaternionic: QuatVec, Lorentz, Boost, components, 𝐤
    using Scri: conformal_factor, transform_ray, aberration, ℐ

    # ℐ is stored on the element as its third type parameter, defaulting to +1 (ℐ⁺); an
    # ℐ⁻-represented element is constructed with the `ℐ` keyword.
    rng = Random.Xoshiro(2718)
    @test BMS(one(Lorentz{Float64}), [1.0im]) isa BMS{Float64,1,1}
    @test ℐ(BMS(one(Lorentz{Float64}), [1.0im]; ℐ=-1)) == -1

    # The ℐ⁻ conformal factor is the antipodal (sign-flipped) Doppler factor: for a z-boost
    # the poles swap relative to ℐ⁺ (κ(+ẑ)=e^{+η}, κ(-ẑ)=e^{-η}), matching 1/(γ(1+v⃗·n̂)),
    # and the default ℐ=+1 reproduces the ℐ⁺ values.
    for T ∈ FloatTypes
        η = T(7//10)
        B = Boost(η, QuatVec{T}(0, 0, 1))
        g₊ = BMS(B, zeros(Complex{T}, 1))
        g₋ = BMS(B, zeros(Complex{T}, 1); ℐ=-1)
        @test abs(conformal_factor(g₊, QuatVec{T}(0, 0, 1)) - exp(-η)) < 20eps(T)
        @test abs(conformal_factor(g₋, QuatVec{T}(0, 0, 1)) - exp(η)) < 20eps(T)
        @test abs(conformal_factor(g₋, QuatVec{T}(0, 0, -1)) - exp(-η)) < 20eps(T)
        for _ ∈ 1:5
            n̂ = random_direction(rng, T)
            # The BMS convenience method must agree with the raw transform_ray at the same ℐ.
            @test conformal_factor(g₋, n̂) ≈ first(transform_ray(g₋.Λ, n̂; ℐ=-1))
        end
    end

    # The action functor realizes the element on its own ℐ: an ℐ⁻-represented element
    # uses the ℐ⁻ conformal factor and direction map, matching the independent oracle
    # act_ref(…; ℐ=-1).  (Attaching ℐ=-1 to the same raw modes reinterprets them in the
    # ℐ⁻ labeling — it does not convert them; `BMS(g; ℐ=-1)` does the conversion.)
    for T ∈ (Float32, Float64, BigFloat)
        for _ ∈ 1:5
            g = random_bms(rng, T; ℓₘₐₓ=2)
            t = 2 * randn(rng, T)
            n̂ = random_direction(rng, T)
            for ℐₖ ∈ (+1, -1)
                gₖ = BMS(g.Λ, g.α; ℐ=ℐₖ)
                t′, n̂′ = gₖ(t, n̂)
                t′ᵣ, n̂′ᵣ = act_ref(g.Λ, g.α, t, n̂; ℐ=ℐₖ)
                @test abs(t′ - t′ᵣ) < tol(T, 2)
                @test maximum(abs, components(n̂′ - n̂′ᵣ)) < 40eps(T)
            end
        end
    end

    # transform_ray at ℐ⁻ agrees with the independently-tested `aberration(…, -1)`
    # direction map (the same check test #9 makes at ℐ⁺), pinning the section convention.
    for _ ∈ 1:8
        R = random_rotation(rng, Float64)
        R′ₚ = random_rotation(rng, Float64)
        v⃗ = (0.7rand(rng)) * random_direction(rng, Float64)
        g = BMS{Float64}(; frame_rotation=R, boost_velocity=v⃗)
        n̂ᵣₑₛₜ = aberration(R * R′ₚ, v⃗, -1)(𝐤)
        _, n̂′ = transform_ray(g.Λ, n̂ᵣₑₛₜ; ℐ=-1)
        @test maximum(abs, components(n̂′ - R′ₚ(𝐤))) < 1e-12
    end
end

@testitem "BMS: representations of α (represent, mixed-convention ops)" tags = [
    :unit, :fast, :validation
] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_bms, random_direction, tol
    using Quaternionic: components
    using Scri: compose, c_α, ℐ

    rng = Random.Xoshiro(3535)

    # Conversion is exact, involutive, and the identity when the conventions already match.
    for T ∈ FloatTypes
        g = random_bms(rng, T; ℓₘₐₓ=3)
        @test BMS(g) === g
        for signs ∈ ((; c_α=-1), (; ℐ=-1), (; c_α=-1, ℐ=-1))
            g′ = BMS(g; signs...)
            @test BMS(g′).α == g.α  # exact round trip
            @test BMS(g′).Λ == g.Λ
            # Different conventions, same abstract transformation:
            @test g′ == g
            @test hash(g′) == hash(g)
            @test g′ ≈ g
        end
        # α ↦ α∘A on modes: odd-ℓ blocks flip, even-ℓ blocks don't.
        g₋ = BMS(g; ℐ=-1)
        @test g₋.α[1] == g.α[1]
        @test g₋.α[2:4] == -g.α[2:4]
        @test g₋.α[5:9] == g.α[5:9]
        @test g₋.α[10:16] == -g.α[10:16]
        # Reinterpreting the raw modes (constructing with ℐ=-1) is NOT the same element,
        # as long as there is odd-ℓ content.
        @test BMS(g.Λ, g.α; ℐ=-1) != g
    end

    # The matched element acts antipodally: if g takes (t, n̂) ↦ (t′, n̂′) on ℐ⁺, then
    # BMS(g; ℐ=-1) takes (t, -n̂) ↦ (t′, -n̂′) on ℐ⁻.  This is the pointwise
    # content of the isomorphism between BMS⁺ and BMS⁻.
    for T ∈ (Float32, Float64, BigFloat)
        for _ ∈ 1:5
            g = random_bms(rng, T; ℓₘₐₓ=2)
            g₋ = BMS(g; ℐ=-1)
            t = 2 * randn(rng, T)
            n̂ = random_direction(rng, T)
            t′, n̂′ = g(t, n̂)
            t′₋, n̂′₋ = g₋(t, -n̂)
            @test abs(t′₋ - t′) < tol(T, 2) * max(abs(t), one(T))
            @test maximum(abs, components(n̂′₋ + n̂′)) < 40eps(T)
        end
    end

    # Mixed-convention composition: conventions the inputs share are kept; where they
    # disagree, both are converted exactly and the result uses the default (+1).
    g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
    g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
    h = compose(g₂, g₁)
    # Disagreeing ℐ: bitwise exact, because the conversion round-trips exactly (pure
    # sign flips) before the identical composition pipeline runs.
    for (a, b) ∈ ((BMS(g₂; ℐ=-1), g₁), (g₂, BMS(g₁; ℐ=-1)))
        h′ = compose(a, b)
        @test ℐ(h′) == 1 && c_α(h′) == 1
        @test h′.α == h.α && h′.Λ == h.Λ
    end
    # Disagreeing c_α likewise.
    h′ = compose(BMS(g₂; c_α=-1), g₁)
    @test c_α(h′) == 1 && ℐ(h′) == 1
    @test h′.α == h.α && h′.Λ == h.Λ
    # Conventions shared by both inputs are kept, and the result is the SAME abstract
    # element as the default-representation composition — `represent` is a group
    # isomorphism.  The c_α flip is algebraically trivial, so it is bitwise exact; the ℐ
    # case runs the composition grid on the antipodal labeling, so the two truncated
    # results agree only to roundoff plus the (mild-boost) band-limit tail.
    h₋α = compose(BMS(g₂; c_α=-1), BMS(g₁; c_α=-1))
    @test c_α(h₋α) == -1
    @test h₋α.α == BMS(h; c_α=-1).α
    h₋ᴵ = compose(BMS(g₂; ℐ=-1), BMS(g₁; ℐ=-1); ℓₘₐₓ=8)
    @test ℐ(h₋ᴵ) == -1
    @test h₋ᴵ ≈ compose(g₂, g₁; ℓₘₐₓ=8) atol = 1e-5

    # inv preserves the element's conventions and inverts on its own ℐ.
    g₋ = BMS(g₁; ℐ=-1)
    g₋⁻¹ = inv(g₋; ℓₘₐₓ=10)
    @test ℐ(g₋⁻¹) == -1 && c_α(g₋⁻¹) == 1
    @test compose(g₋⁻¹, g₋; ℓₘₐₓ=10) ≈ one(g₋) atol = 1e-7
end

@testitem "BMS: κ cocycle via the 4-vector action" tags = [:unit, :fast, :validation] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_direction, random_lorentz, ray_map
    using Quaternionic: components

    # κ_{Λ₂Λ₁}(𝐤) = κ₂(Λ₁𝐤) κ₁(𝐤), and the direction maps compose accordingly.  This
    # uses only the oracle `ray_map` — no Scri code — so it validates the geometric
    # consistency that `compose` relies on, independently of the composition law.
    rng = Random.Xoshiro(112)
    for T ∈ FloatTypes
        for _ ∈ 1:5
            Λ₁ = random_lorentz(rng, T)
            Λ₂ = random_lorentz(rng, T)
            n̂ = random_direction(rng, T)
            κ₁, n̂₁ = ray_map(Λ₁, n̂)
            κ₂, n̂₂ = ray_map(Λ₂, n̂₁)
            κ₂₁, n̂₂₁ = ray_map(Λ₂ * Λ₁, n̂)
            @test abs(κ₂₁ - κ₂ * κ₁) < 40eps(T) * abs(κ₂₁)
            @test maximum(abs, components(n̂₂₁ - n̂₂)) < 40eps(T)
        end
    end
end

@testitem "BMS: direction map and rotor_from_direction" tags = [:unit, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_direction, random_rotation
    using Quaternionic: QuatVec, Lorentz, absvec, components, 𝐤
    using Scri: rotor_from_direction, transform_ray

    rng = Random.Xoshiro(213)
    for T ∈ FloatTypes
        # rotor_from_direction points where it should, including at the poles.
        for n̂ ∈ (
            random_direction(rng, T),
            random_direction(rng, T),
            QuatVec{T}(0, 0, 1),
            QuatVec{T}(0, 0, -1),
            QuatVec{T}(1, 0, 0),
        )
            R = rotor_from_direction(n̂)
            @test maximum(abs, components(R(𝐤) - n̂)) < 20eps(T)
        end
        # The direction map produces unit vectors.
        for _ ∈ 1:5
            Λ = BMSTestSetup.random_lorentz(rng, T; βmax=9//10)
            n̂ = random_direction(rng, T)
            _, n̂′ = transform_ray(Λ, n̂)
            @test abs(absvec(n̂′) - 1) < 40eps(T)
        end
        # For pure rotations, the direction map is just the rotation.
        for _ ∈ 1:5
            Q = random_rotation(rng, T)
            n̂ = random_direction(rng, T)
            _, n̂′ = transform_ray(Lorentz(Q), n̂)
            @test maximum(abs, components(n̂′ - Q(n̂))) < 40eps(T)
        end
    end
end

@testitem "BMS: boost_velocity convention matches compute_t′" tags = [
    :unit, :fast, :validation, :integration
] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_direction
    using LinearAlgebra: dot
    using Quaternionic: absvec
    using Scri: conformal_factor

    # The `boost_velocity` keyword must mean exactly what `v⃗` means in `transform!` and
    # `compute_t′`: there, the time scaling is t′ = κ(t - α) with 1/κ = γ(1 - v⃗⋅n̂)
    # evaluated at the rest-frame direction n̂.  This pins the sign with which the boost
    # enters the stored Lorentz rotor.
    rng = Random.Xoshiro(414)
    for T ∈ FloatTypes
        for _ ∈ 1:8
            v⃗ = (T(8//10) * rand(rng, T)) * random_direction(rng, T)
            β = absvec(v⃗)
            γ = 1 / √(1 - β^2)
            g = BMS{T}(; boost_velocity=v⃗)
            n̂ = random_direction(rng, T)
            κ = conformal_factor(g, n̂)
            κ_expected = 1 / (γ * (1 - dot(vec(v⃗), vec(n̂))))
            @test abs(κ - κ_expected) < 40eps(T) * abs(κ_expected)
        end
    end
end

@testitem "BMS: direction map vs aberration" tags = [
    :unit, :fast, :validation, :integration
] setup = [BMSTestSetup] begin
    import Random
    using Scri: aberration, transform_ray
    using .BMSTestSetup: random_direction, random_rotation
    using Quaternionic: QuatVec, Rotor, components, 𝐤

    # transform! evaluates rest-frame data at Rₚ = aberration(R * R′ₚ, v⃗) to fill the
    # transformed-frame pixel R′ₚ.  The BMS ray map must therefore take the rest-frame
    # direction Rₚ(𝐤) to the transformed-frame direction R′ₚ(𝐤).  This cross-validates
    # the stored Lorentz rotor (both the boost sign and the rotation/boost order) against
    # the independently-derived and independently-tested `aberration`.
    rng = Random.Xoshiro(515)
    for _ ∈ 1:10
        R = random_rotation(rng, Float64)
        R′ₚ = random_rotation(rng, Float64)
        v⃗ = (0.7rand(rng)) * random_direction(rng, Float64)
        g = BMS{Float64}(; frame_rotation=R, boost_velocity=v⃗)
        n̂ᵣₑₛₜ = aberration(R * R′ₚ, v⃗)(𝐤)
        _, n̂′ = transform_ray(g.Λ, n̂ᵣₑₛₜ)
        @test maximum(abs, components(n̂′ - R′ₚ(𝐤))) < 1e-12
    end
end

@testitem "BMS: action functor and c_α" tags = [:unit, :fast] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_direction, random_bms, act_ref, tol
    using Quaternionic: components
    using Scri: c_α, conformal_factor

    rng = Random.Xoshiro(616)
    for T ∈ (Float32, Float64, BigFloat)
        for _ ∈ 1:5
            g = random_bms(rng, T; ℓₘₐₓ=2)  # c_α = +1 by default
            t = 2 * randn(rng, T)
            n̂ = random_direction(rng, T)
            t′, n̂′ = g(t, n̂)
            t′ᵣ, n̂′ᵣ = act_ref(g.Λ, g.α, t, n̂; c_α=1)
            @test abs(t′ - t′ᵣ) < tol(T, 2)
            @test maximum(abs, components(n̂′ - n̂′ᵣ)) < 40eps(T)
            # c_α is the element's type parameter: the c_α = -1 element flips only the
            # supertranslation term, leaving the direction (and conformal factor) alone.
            g₋ = BMS(g.Λ, g.α; c_α=-1)
            @test c_α(g₋) == -1
            t′₋, n̂′₋ = g₋(t, n̂)
            @test n̂′₋ == n̂′
            κ = conformal_factor(g, n̂)
            @test abs((t′₋ + t′) - 2κ * t) < tol(T, 2) * max(abs(t), one(T))
        end
    end
end

###
### Group operations
###

# Modes are padded (or truncated) to the requested ℓₘₐₓ.
function resize_modes(α::Vector{Complex{T}}, ℓₘₐₓ::Int) where {T}
    N = (ℓₘₐₓ + 1)^2
    if N == length(α)
        return copy(α)
    elseif N > length(α)
        out = zeros(Complex{T}, N)
        out[1:length(α)] = α
        return out
    else
        return α[1:N]
    end
end

is_identity_rotor(Λ::Lorentz) = Λ == one(Λ) || Λ == -one(Λ)

@doc raw"""
    compose(g₂, g₁; ℓₘₐₓ, ℓʷ)
    g₂ * g₁
    g₂ ∘ g₁

Compose two BMS elements: the result acts as `g₁` *first*, then `g₂`:

```math
(Λ₂, α₂)\,(Λ₁, α₁) = \left(Λ₂ Λ₁,\; α₁(𝐤) + α₂(Λ₁𝐤)/κ₁(𝐤)\right).
```

Because `κ₁` and the mapped direction `Λ₁𝐤` are realized on the labeled celestial sphere,
the composed supertranslation depends on the labeling — which each element carries as its
`I` type parameter (see [Representations of the supertranslation](@ref
bms_representations)).  The inputs may use any mix of conventions
(`A`, `I`): each convention the two inputs *share* is kept for the result, and each
convention on which they *disagree* is resolved by converting both inputs — exactly — to
the default (`+1`).  `c_α` does not otherwise enter the group law.

The composed supertranslation is computed pointwise on a spherical grid of bandwidth `ℓʷ`
and then re-expanded in spherical harmonics, keeping modes up to `ℓₘₐₓ` (which defaults
to the larger of the two inputs' resolutions; `ℓʷ` defaults to `2ℓₘₐₓ + 1` or the largest
input resolution, whichever is greater).

!!! note "Band limits"
    The factor `1/κ₁` is *exactly* band-limited to `ℓ ≤ 1`, and for pure rotations
    `α₂∘Λ₁` has the same bandwidth as `α₂` — so when `Λ₁` involves no boost, the result
    is exact (up to roundoff) whenever `ℓʷ ≥ max(ℓ₁, ℓ₂)`, which the defaults guarantee.
    When `Λ₁` involves a boost, however, `α₂∘Λ₁` has unbounded bandwidth (with an
    exponentially decaying tail), so the result is *convergent in `ℓʷ`* rather than
    exact, and the modes above `ℓₘₐₓ` are simply discarded.  For large rapidities,
    increase `ℓʷ` (and possibly `ℓₘₐₓ`).

Two cases short-circuit the grid entirely and are exact: when `α₂ = 0` (so the result is
just `(Λ₂Λ₁, α₁)`), and when `Λ₁ = ±1` (so `κ₁ ≡ 1` and the supertranslations simply
add).
"""
function compose(
    g₂::BMS{T,A,I},
    g₁::BMS{T,A,I};
    ℓₘₐₓ::Int=max(ell_max(g₂), ell_max(g₁)),
    ℓʷ::Int=max(2ℓₘₐₓ + 1, ell_max(g₂), ell_max(g₁)),
) where {T<:Real,A,I}
    @assert ℓʷ ≥ ℓₘₐₓ "Working bandwidth ℓʷ=$ℓʷ must be at least ℓₘₐₓ=$ℓₘₐₓ"
    # Quaternionic's rotor products renormalize (by the spinor norm), which perturbs the
    # last ulps even when multiplying by the identity; short-circuiting ±1 factors keeps
    # the group identities exact.
    Λ = if is_identity_rotor(g₁.Λ)
        g₁.Λ == one(g₁.Λ) ? g₂.Λ : -g₂.Λ
    elseif is_identity_rotor(g₂.Λ)
        g₂.Λ == one(g₂.Λ) ? g₁.Λ : -g₁.Λ
    else
        g₂.Λ * g₁.Λ
    end
    α = if all(iszero, g₂.α)
        # (Λ₂, 0)(Λ₁, α₁) = (Λ₂Λ₁, α₁) — exact.
        resize_modes(g₁.α, ℓₘₐₓ)
    elseif is_identity_rotor(g₁.Λ)
        # κ₁ ≡ 1 and Λ₁𝐤 = 𝐤, so the supertranslations simply add — exact.
        resize_modes(g₁.α, ℓₘₐₓ) .+ resize_modes(g₂.α, ℓₘₐₓ)
    else
        Rₚ = golden_ratio_spiral_rotors(0, ℓʷ, T)
        α₁ₚ = real.(ₛ𝐘(0, ell_max(g₁), T, Rₚ) * g₁.α)
        R′ₚ = similar(Rₚ)
        κ₁ₚ = Vector{T}(undef, length(Rₚ))
        for p ∈ eachindex(Rₚ)
            κ₁, n̂′ = transform_ray(g₁.Λ, Rₚ[p](𝐤); ℐ=I)
            κ₁ₚ[p] = κ₁
            R′ₚ[p] = rotor_from_direction(n̂′)
        end
        α₂ₚ = real.(ₛ𝐘(0, ell_max(g₂), T, R′ₚ) * g₂.α)
        f = @. complex(α₁ₚ + α₂ₚ / κ₁ₚ)
        modes = lu(ₛ𝐘(0, ℓʷ, T, Rₚ)) \ f
        modes[1:((ℓₘₐₓ + 1) ^ 2)]
    end
    return BMS{T,A,I}(Λ, α)
end

function compose(g₂::BMS{T2,A₂,I₂}, g₁::BMS{T1,A₁,I₁}; kwargs...) where {T1,T2,A₁,A₂,I₁,I₂}
    # Mixed conventions: keep each convention the inputs share; where they disagree,
    # convert both — exactly — to the default (+1) and compose there.
    A = A₁ == A₂ ? A₁ : 1
    I = I₁ == I₂ ? I₁ : 1
    T = promote_type(T1, T2)
    return compose(BMS{T,A,I}(g₂), BMS{T,A,I}(g₁); kwargs...)
end

Base.:*(g₂::BMS, g₁::BMS) = compose(g₂, g₁)
Base.:∘(g₂::BMS, g₁::BMS) = compose(g₂, g₁)

function Base.one(::Type{BMS{T,A,I}}) where {T<:Real,A,I}
    return BMS{T,A,I}(one(Lorentz{T}), zeros(Complex{T}, 1))
end
Base.one(::Type{BMS{T,A}}) where {T<:Real,A} = one(BMS{T,A,1})
Base.one(::Type{BMS{T}}) where {T<:Real} = one(BMS{T,1,1})
Base.one(::Type{BMS}) = one(BMS{Float64,1,1})
Base.one(g::BMS{T,A,I}) where {T<:Real,A,I} = one(BMS{T,A,I})
Base.isone(g::BMS) = is_identity_rotor(g.Λ) && all(iszero, g.α)

@doc raw"""
    inv(g::BMS; ℓₘₐₓ, ℓʷ)

The inverse BMS element,

```math
(Λ, α)^{-1} = (1, -α) ∘ (Λ^{-1}, 0)
            = \left(Λ^{-1},\; -α(Λ^{-1}𝐤) / κ_{Λ^{-1}}(𝐤)\right),
```

computed through [`compose`](@ref); the result inherits `g`'s conventions (`A`, `I`),
and the null-infinity representation `I` governs the geometry internally.  `compose`'s
band-limit caveats (and `ℓₘₐₓ`/`ℓʷ` keywords) apply here too: the inverse of an element
whose Lorentz part involves a boost is exact only in the limit of large `ℓʷ` and `ℓₘₐₓ`.
The Lorentz part `inv(Λ)` is the exact GA reverse.
"""
function Base.inv(
    g::BMS{T,A,I}; ℓₘₐₓ::Int=ell_max(g), ℓʷ::Int=max(2ℓₘₐₓ + 1, ell_max(g))
) where {T<:Real,A,I}
    return compose(
        BMS{T,A,I}(one(Lorentz{T}), -g.α),
        BMS{T,A,I}(inv(g.Λ), zeros(Complex{T}, 1));
        ℓₘₐₓ,
        ℓʷ,
    )
end

###
### Inline tests: group operations
###

@testitem "BMS: exact group identities (fast paths)" tags = [:unit, :fast] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_bms, random_lorentz, random_α
    using Quaternionic: Lorentz
    using Scri: is_identity_rotor, time_translation

    rng = Random.Xoshiro(717)
    for T ∈ FloatTypes
        g = random_bms(rng, T; ℓₘₐₓ=2)
        e = one(BMS{T})

        # Identity laws hold exactly, thanks to the fast paths.
        @test g * e == g
        @test e * g == g

        # The defining decomposition (Λ, α) = (Λ, 0) ∘ (1, α), exactly.
        gΛ = BMS(g.Λ, zeros(Complex{T}, 1))
        gα = BMS(one(Lorentz{T}), g.α)
        @test gΛ * gα == g

        # Acting with a Lorentz transformation *after* a supertranslation never changes
        # the supertranslation, for any Λ.
        Λ₂ = random_lorentz(rng, T)
        h = BMS(Λ₂, zeros(Complex{T}, 1)) * g
        @test h.α == g.α
        @test h.Λ == Λ₂ * g.Λ

        # Supertranslations form an abelian subgroup; composition is exact addition.
        α₁ = random_α(rng, T, 2)
        α₂ = random_α(rng, T, 3)  # deliberately different ℓₘₐₓ
        s₁ = BMS(one(Lorentz{T}), α₁)
        s₂ = BMS(one(Lorentz{T}), α₂)
        s₂₁ = s₂ * s₁
        @test s₂₁ == s₁ * s₂
        @test s₂₁.α == [α₁; zeros(Complex{T}, 16 - 9)] .+ α₂
        @test is_identity_rotor(s₂₁.Λ)
        # Pure time translations compose additively, with no other modes appearing.
        # (Mode addition is exact; the δt accessor only rounds in the final division.)
        t₁ = BMS{T}(; time_translation=T(3//2))
        t₂ = BMS{T}(; time_translation=T(7//4))
        @test (t₂ * t₁) == BMS(one(Lorentz{T}), t₁.α .+ t₂.α)
        @test time_translation(t₂ * t₁) ≈ T(3//2) + T(7//4) atol = 8eps(T)
        @test all(iszero, (t₂ * t₁).α[2:end])
        # Inverses of the special subgroups are exact.
        @test inv(s₁) == BMS(one(Lorentz{T}), -α₁)
        @test inv(gΛ) == BMS(inv(g.Λ), zeros(Complex{T}, 1))
    end
end

@testitem "BMS: Lorentz sector is Quaternionic multiplication" tags = [
    :unit, :fast, :validation
] setup = [BMSTestSetup] begin
    import Random
    using .BMSTestSetup: FloatTypes, random_lorentz

    # The Lorentz sector of the group must reduce to plain rotor multiplication in
    # Quaternionic — which is tested independently over there.
    rng = Random.Xoshiro(818)
    for T ∈ FloatTypes
        Λ₁ = random_lorentz(rng, T)
        Λ₂ = random_lorentz(rng, T)
        g₁ = BMS(Λ₁, zeros(Complex{T}, 1))
        g₂ = BMS(Λ₂, zeros(Complex{T}, 1))
        h = g₂ * g₁
        @test h.Λ == Λ₂ * Λ₁
        @test all(iszero, h.α)
        @test inv(g₁).Λ == conj(Λ₁)
        @test inv(g₁).Λ == inv(Λ₁)
    end
end

@testitem "BMS: composed α matches the pointwise oracle" tags = [:unit, :validation] setup = [
    BMSTestSetup
] begin
    import Random
    using .BMSTestSetup: random_bms, random_direction, ray_map, α_value
    using Scri: compose

    # Evaluate the composed supertranslation at random (off-grid) directions and compare
    # with α₁(n̂) + α₂(n̂′)/κ₁(n̂) computed directly from raw geometry and the *input*
    # modes.  Mild boosts keep the truncation error far below the test tolerance.
    rng = Random.Xoshiro(919)
    for _ ∈ 1:5
        g₁ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
        g₂ = random_bms(rng, Float64; ℓₘₐₓ=2, βmax=1//5)
        h = compose(g₂, g₁; ℓₘₐₓ=10)
        for _ ∈ 1:5
            n̂ = random_direction(rng, Float64)
            κ₁, n̂′ = ray_map(g₁.Λ, n̂)
            expected = α_value(g₁.α, n̂) + α_value(g₂.α, n̂′) / κ₁
            @test abs(α_value(h.α, n̂) - expected) < 1e-6
        end
    end
end
