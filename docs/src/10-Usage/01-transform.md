```@meta
CurrentModule = Scri
```

# Transforming Waveforms

The central function of this package is [`transform!`](@ref), which
applies a BMS transformation to asymptotic data.  This page describes
the data it expects, what it does to them, and what it returns; the
following pages describe how to [specify the transformation](@ref
"Specifying a BMS Transformation") and the [data components](@ref
"Data Components"), and how to [choose the angular resolution](@ref
"Choosing ``ℓ_\mathrm{max}``").

## The data array

The waveform data are passed as a single three-dimensional array of
complex numbers with dimensions ``(Nᵐ, Nᵗ, Nᵈ)``:

* The **first** dimension holds the spin-weighted spherical-harmonic
    *mode weights*, ordered by increasing ``ℓ`` starting from ``ℓ =
    0``, and by increasing ``m`` from ``-ℓ`` to ``+ℓ`` within each
    ``ℓ`` — so the mode ``(ℓ, m)`` sits at index ``ℓ² + ℓ + m + 1``.
    The band limit is implicit in the size: ``ℓ_\mathrm{max} =
    \sqrt{Nᵐ} - 1``.  Note that the ``ℓ < |s|`` modes are present
    *even for fields of nonzero spin weight* ``s``; they are ignored
    on input (they should be zero), and on output they hold the
    null-space diagnostic ``ξ`` of the [augmented SSHT](@ref
    "Augmented Direct SSHT") — but could also be safely ignored.
* The **second** dimension runs over the time samples, matching the
    `t` argument.
* The **third** dimension runs over the physical field components —
    strain, Weyl components, and so on — in the order declared by the
    [`DataComponents`](@ref) descriptor.

The array's element type must be complex, with a real type at least as
wide (referring to the number of bytes, or precision) as every other
input's; `transform!` throws an `ArgumentError` otherwise, since it
cannot widen the array it modifies in place.

!!! note

    Here, "first" means the most rapidly varying index in memory
    (the leftmost index in Julia's column-major order).  If you are
    transferring an array from Python, you probably won't be copying
    it; it will be effectively transposed for you, so the *last*
    dimension of the corresponding NumPy array will be the mode
    index, the middle dimension will be time, and the first dimension
    will be data component.

## In-place semantics

The exclamation mark is Julia's convention for functions that modify
their arguments: `transform!` overwrites `data` with the transformed
mode weights (and also returns it, for convenience).  This avoids
allocating a second copy of what is often a very large array.  If you
need to keep the original, copy it first:

```julia
data′, t′ = transform!(copy(data), t, v⃗, R, α, dc)
```

The band limit of the output is the band limit of the array you pass
in.  Because BMS transformations push power to higher ``ℓ``, you will
usually want to *pad first* — embed the data in a larger array before
transforming — as described under [Choosing
``ℓ_\mathrm{max}``](@ref).

## Calling forms

The fully explicit form takes the transformation as separate parts —
boost velocity `v⃗` (a `QuatVec`), frame rotation `R` (a `Rotor`), and
supertranslation mode weights `α` — along with the
[`DataComponents`](@ref) descriptor:

```julia
dc = DataComponents(:h, :ψ₄)
data′, t′ = transform!(data, t, v⃗, R, α, dc)
```

A keyword form builds the descriptor for you, accepting flexible
spellings of the component names and the `ℐ` and `conventions`
parameters:

```julia
data′, t′ = transform!(data, t, v⃗, R, α; data_components=("h", "Psi4"))
```

And a [`BMS`](@ref)-element form takes the whole transformation as a
single group element — the natural choice when composing
transformations or reusing them:

```julia
g = BMS(; boost_velocity=v⃗, frame_rotation=R, supertranslation=α)
data′, t′ = transform!(data, t, g, dc)
```

All three forms are described in detail on the next page,
[Specifying a BMS Transformation](@ref).

## The output time grid

Boosts and supertranslations mix time and direction, so a slice of
constant transformed time ``t'`` is "tilted" with respect to the
input slices, and slices near the ends of the input span are no
longer covered by data over the whole sphere.  `transform!` therefore
constructs a new time grid `t′` — with the same number of samples as
`t`, uniformly covering the largest span for which every direction on
the sphere has data — and returns it along with the data:

```julia
data′, t′ = transform!(data, t, v⃗, R, α, dc)
```

The transformed mode weights `data′[:, j, :]` are the data on the
slice of constant transformed time `t′[j]`.

If you need a specific output grid — most commonly in an optimization
loop over BMS parameters, where a fixed grid makes iterations
directly comparable — pass it with the `t′` keyword:

```julia
data′, _ = transform!(data, t, v⃗, R, α, dc; t′=my_grid)
```

The supplied grid must be strictly increasing and must stay within
the valid span (an `ArgumentError` explains the limits if not).

## Differentiability

`transform!` can be differentiated end-to-end with
[ForwardDiff](https://juliadiff.org/ForwardDiff.jl/) with respect to
the BMS parameters — boost velocity, rotation, and supertranslation
(hence also time and space translations).  Two things to know:

* **Fix the output grid.**  Pass the `t′` keyword; otherwise the
    default output grid itself moves with the parameters, and the
    derivative of the data at output index `j` mixes the field
    derivative with the motion of the grid.  With a fixed `t′`, you
    get clean partial derivatives of the field at fixed output times —
    which will usually be what you want in an optimization loop.  But
    this means choosing the range of `t′` values to be small enough
    that you are sure it will always be within the valid span for
    every boost or supertranslation parameter value the optimization
    loop will try.  See [`Scri.compute_t′_bounds`](@ref) for a utility
    that computes the valid span for a given input grid and BMS
    element.
* **Load ForwardDiff.**  Doing so activates a package extension that
    keeps dual numbers out of the parameter-independent pixel grid and
    its matrix factorizations; without it, derivatives with respect to
    the rotation come out as `NaN`.

The `data` array must have the dual type (e.g.,
`Complex{typeof(θ)}.(data)` for a differentiated parameter `θ` inside
the differentiated function), since `transform!` works in place.
Derivatives agree with finite differences at the `1e-8` level; see the
testitem `"transform!: ForwardDiff derivatives with respect to BMS
parameters"`.

## Checking the result

The transformation is spectrally accurate in angle but limited by the
band limit ``ℓ_\mathrm{max}`` and by interpolation error in time.
Use [`diagnostics`](@ref) to compute per-``ℓ`` power monitors of the
input and output, and see [Choosing ``ℓ_\mathrm{max}``](@ref) and
[Spline Errors](@ref) for how to interpret them.
