```@meta
CurrentModule = Scri
```

# Specifying a BMS Transformation

A general BMS element is the composition of three simpler
transformations, which this package applies in a fixed order:
supertranslation ``α`` first, then spatial rotation ``R``, then boost
``v⃗``.  Each of the three is specified in the coordinate system of the
initial, inertial observer ``A``, and each is described below as if it
were acting alone.  Throughout, we use the *passive* convention: the
spacetime events are fixed, and it is the coordinates used to describe
them that change (see [BMS Transformations](@ref) for the theory, and
[The BMS Group](@ref bms_group) for the group structure).

## The three parts

The **supertranslation** ``α`` shifts the (retarded or advanced) time
coordinate by a real, angle-dependent function of the direction in
``A``'s frame,

```math
t'(t, 𝐤) = t - c_α α(𝐤),
```

where ``c_α`` is the [sign convention](@ref conventions-type) of the
supertranslation law.  The null ray ``𝐤`` — and so the angular
coordinates — is untouched.

The **spatial rotation** ``R`` is a rotor ``R ∈ \mathrm{Spin}(3)`` that
rotates the angular coordinate system.  Observer ``B``'s spatial basis
is related to ``A``'s by ``𝐞' = R\, 𝐞\, R̃``, so a null ray at ``𝐤`` in
``A``'s frame appears at ``R\, 𝐤\, R̃`` in ``B``'s.  The time coordinate
is unchanged.

The **boost** is given by a velocity ``v⃗`` in ``A``'s frame.  Observer
``B``'s worldline is parallel to

```math
𝐭' = γ(𝐭 + v⃗),
```

and the boost rescales the time coordinate by the conformal factor,

```math
t' = κ(𝐤)\, t,
\qquad
κ(𝐤) = \frac{1}{γ(1 - ℐ\, v⃗ ⋅ k̂)},
```

with ``ℐ = ±1`` selecting [future or past null infinity](@ref
scri_pm_conventions).  Composing all three, the BMS element
``(Λ = \text{boost} ∘ R,\ α)`` acts on ``A``'s coordinates ``(t, 𝐤)``
as

```math
t' = κ(Λ, 𝐤)\,\bigl[t - c_α α(𝐤)\bigr],
\qquad
𝐤' = \frac{Λ𝐤}{(Λ𝐤)^0},
```

where ``κ(Λ, 𝐤) = n⁰/(Λ𝐤)⁰`` is the conformal factor of the combined
Lorentz transformation on the null ray ``𝐤 = (1, ℐ k̂)``.

## Calling `transform!`

There are two equivalent ways to apply this to waveform data.  The
direct form takes the three parts as positional arguments, with the
[`DataComponents`](@ref Scri.DataComponents) descriptor either supplied
explicitly or built from keywords:

```julia
# descriptor built from keywords
data′, t′ = transform!(data, t, v⃗, R, α; data_components=("h", "Psi4"), ℐ=+1)

# descriptor supplied directly
dc = DataComponents(:h, :ψ₄)
data′, t′ = transform!(data, t, v⃗, R, α, dc)
```

Here `v⃗` is a `QuatVec`, `R` is a `Rotor`, and `α` is the
supertranslation as a vector of spin-0 mode weights (its band limit may
be smaller than the data's, but not larger).  The supertranslation is
automatically symmetrized to enforce the reality condition ``α_{ℓ,-m} =
(-1)^m \bar α_{ℓ,m}``.

## Building a `BMS` element

For anything beyond a one-off call — composing transformations,
translating in space or time, or reusing an element — construct a
[`BMS`](@ref) value and hand it to the element form of `transform!`:

```julia
g = BMS(;
    boost_velocity = v⃗,
    frame_rotation = R,
    supertranslation = α,   # spin-0 mode weights, optional
    time_translation = δt,  # ℓ=0 supertranslation, optional
    space_translation = δx⃗, # ℓ=1 supertranslation, optional
)
data′, t′ = transform!(data, t, g, dc)
```

so that `transform!(data, t, g, dc)` reproduces the action of `g` on
the data.  Pure spacetime translations are just the ``ℓ = 0`` and
``ℓ = 1`` parts of a supertranslation, and the `time_translation` /
`space_translation` keywords set them by their physical parameters
``δt`` and ``δx⃗`` (see the [BMS decomposition](@ref
"Decomposition of BMS")).  Elements compose with
[`compose`](@ref Scri.compose), and their conformal factor on a
direction is available from [`conformal_factor`](@ref
Scri.conformal_factor).

The `BMS` element also carries the two convention signs `c_α` and `ℐ`
that fix how its supertranslation is [represented](@ref
bms_representations); operations that combine elements reconcile these
automatically.

```@docs; canonical=false
Scri.BMS
```
