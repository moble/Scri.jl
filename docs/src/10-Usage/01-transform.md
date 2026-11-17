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
choosing_ell_max).

If you are optimizing the BMS parameters to minimize the difference
with another waveform, you may want to use the
[`transform_objective`](@ref) function instead, which is more
efficient in both time and memory usage.  See the [Automatic
Differentiation](@ref) page for details.  Still, that function is
based on `transform!` — at least in spirit — so it is helpful to
understand `transform!` before trying to use `transform_objective`.

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
"wide" (referring to the number of bytes, or precision) as every other
input's; `transform!` throws an `ArgumentError` otherwise, since it
cannot widen the array it modifies in place.

!!! note

    Here, the "first" dimension refers to the most rapidly varying
    index in memory (the leftmost index in Julia's column-major
    order).  If you are transferring an array from Python, you
    probably won't be copying it; it will be effectively transposed
    for you, so the *last* dimension of the corresponding NumPy array
    will be the mode index, the middle dimension will be time, and
    the first dimension will be data component.

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
transforming — as described under [Choosing ``ℓ_\mathrm{max}``](@ref
choosing_ell_max).

## Reusing the setup across calls

Part of every call to `transform!` depends on no BMS parameter: the
grid on which the transformed data are analyzed back into modes, and
the factorization of the analysis matrix there for each spin weight in
the data.  The cost of that setup grows steeply with
``ℓ_\mathrm{max}``, but not with the number of time samples, as the
rest of the call's cost does.  So the setup is a large share of the
cost at a high ``ℓ_\mathrm{max}`` for a short waveform, and a small
share for a long one.  When transforming many times at one band limit,
precision, and set of components, as in a scan over parameters, build
a [`TransformPlan`](@ref) once and pass it along:

```julia
plan = TransformPlan(data, dc)
for (v⃗, R, α) ∈ parameters
    data′, t′ = transform!(copy(data), t, v⃗, R, α, dc; plan)
    # ...
end
```

The result is identical to what the call would compute without the
plan.  How much time the plan saves depends on the length of the
waveform.  For the strain alone, on 8 threads, a call is this many
times faster with the plan than without it:

| ``ℓ_\mathrm{max}`` | 400 time samples | 4000 | 20,000 |
| :---: | :---: | :---: | :---: |
| 16 | 2.3 | 1.2 | 1.1 |
| 32 | 5.0 | 1.6 | 1.3 |

Building the plan takes about as long as the setup that it replaces:
6 ms at ``ℓ_\mathrm{max} = 16``, and 0.14 s at ``ℓ_\mathrm{max} = 32``.

Most of the time of a call at a high ``ℓ_\mathrm{max}`` starts with
modes, "synthesizes" to pixel values, (does some processing and) then
"analyzes" back to modes — which are dense matrix products and solves.
Those run on BLAS's own threads, however many threads Julia has been
started with, so starting Julia with more threads speeds up only the
work done pixel by pixel.  With 20,000 time samples and 8 BLAS
threads, that work is 2.3 to 4.5 times faster with 8 Julia threads
than with 1, but a whole call is only about 1.2 times faster, or 1.25
to 1.5 times with a plan.  With 8 Julia threads, the matrix products
and solves take 42% of the time at ``ℓ_\mathrm{max} = 8`` for the
strain alone, and 78% at ``ℓ_\mathrm{max} = 16`` for the strain and
all five Weyl components.  For those, the BLAS library and its number
of threads matter far more than Julia's.

## Signatures

The fully explicit form takes the transformation as separate parts —
boost velocity `v⃗` (a `QuatVec`), frame rotation `R` (a `Rotor`), and
supertranslation mode weights `α` — along with the
[`DataComponents`](@ref) descriptor.  Here `v⃗` is the velocity of the
new observer relative to the old one, and the new axes are the old
axes turned by `R` before the boost, both given in the old frame's
components: the new legs are ``𝐞'_μ = Λ\,𝐞_μ\,Λ̃`` with the frame
rotor ``Λ = B(v⃗)\,R``.  For example, with `R = Rotor(√2, 0, 0, √2)`,
a quarter turn about ``z``, the new ``x`` axis is the old ``y`` axis,
and the event with old coordinates ``(0, 1, 0, 0)`` has new
coordinates ``(0, 0, -1, 0)``.  The precise meaning of both is fixed
on the page [Frames and Components](@ref frames_and_components).

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

All three forms are described in detail on the next page, [Specifying
a BMS Transformation](@ref).

## The output time grid

Boosts and supertranslations mix time and space, so a slice of
constant transformed time ``t'`` is "tilted" with respect to the input
slices.  But each slice needs complete data over the whole sphere to
supply the spherical-harmonic synthesis; the first such slice of
``t'`` will be later than the first slice of ``t``, and the last slice
of ``t'`` will be earlier than the last slice of ``t``.  This is a
necessary feature of any transformation with a nonzero boost or a
supertranslation that varies over the sphere.  (A pure time
translation does not produce a tilt; it only relabels the slices.)
`transform!` therefore constructs a new time grid `t′` by *scaling*
the input `t` to cover the largest possible span.  In particular, `t′`
has the same number of samples as `t`.  This `t′` is returned along
with the data:

```julia
data′, t′ = transform!(data, t, v⃗, R, α, dc)
```

It is *possible* to pass a custom `t′` array with the corresponding keyword:

```julia
data′, _ = transform!(data, t, v⃗, R, α, dc; t′=my_grid)
```

However, the supplied grid must be strictly increasing, must stay
within the valid span of ranges, and *must have the same number of
samples as `t`*.  The latter is not a fundamental limitation, but the
general case is not yet implemented, because the motivation is not yet
clear.  The typical use case for this is probably covered by the
[`transform_objective`](@ref) function; see the [Automatic
Differentiation](@ref) page for details.

## Accuracy

The transformation is spectrally accurate in angle but limited by the
band limit ``ℓ_\mathrm{max}`` and by interpolation error in time.  Use
[`diagnostics`](@ref) to compute per-``ℓ`` power monitors of the input
and output, and see [Choosing ``ℓ_\mathrm{max}``](@ref
choosing_ell_max) and [Spline Errors](@ref) for how to interpret them.
