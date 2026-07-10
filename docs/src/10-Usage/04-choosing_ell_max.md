```@meta
CurrentModule = Scri
```

# Choosing ``ℓ_\mathrm{max}``

A BMS transformation mixes energy into higher ``ℓ`` modes.  The input
waveform will have a fairly low *initial* band limit
``ℓ_{\mathrm{max}, i}`` — often around 8 or less.  Any nontrivial BMS
transformation (other than a pure rotation or time translation) will
produce output with comparable amounts of energy at higher ``ℓ``.
Strictly speaking, the output is *not* band-limited at all; any
nonzero boost or supertranslation (other than a pure time translation)
will produce *some* energy at arbitrarily high ``ℓ``.  However,
spin-weighted spherical harmonics should converge exponentially, so we
just have to choose some finite ``ℓ_\mathrm{max}`` that is "good
enough" for our purposes.

Choosing a large enough ``ℓ_\mathrm{max}`` is probably the single most
important accuracy decision in using this package.  On the other hand,
the cost of the transform grows steeply with ``ℓ_\mathrm{max}`` —
scaling as ``(ℓ_\mathrm{max}+1)^4`` in memory, and possibly even worse
in time.  Therefore, we want to balance the accuracy against the cost
so that ``ℓ_\mathrm{max}`` is large enough but not too large.

## Specifying ``ℓ_\mathrm{max}``

You may have noticed that there is no explicit argument to most of the
functions in this package that specifies ``ℓ_\mathrm{max}``.  This is
because it is implicit in the size of the `data` array: if the first
dimension of the array has size ``Nᵐ``, then the maximum ``ℓ`` is read
off as ``ℓ_\mathrm{max} = \sqrt{Nᵐ} - 1``.

The [`transform!`](@ref) function works **in place**.  The band limit
of its output is exactly the band limit of the array you hand it — the
maximum ``ℓ`` is read off the first dimension as ``ℓ_\mathrm{max} =
\sqrt{Nᵐ} - 1`` — and the transform cannot grow that array to make
room for the mixing it produces.  Any energy that the transformation
pushes above ``ℓ_\mathrm{max}`` has nowhere to go: it is lost, and
(worse) it aliases back down onto the resolved modes, corrupting them.

So the workflow is to *pad first*.  Beginning with data at its native
band limit, embed it in a larger array whose ``ℓ_\mathrm{max}`` is big
enough to hold the transformed energy, and only then transform.  The
extra modes must be zero on input, and will be filled in during the
transform.

## How much is enough?

There is no single formula, but two effects set the scale.

* **Supertranslations.**  The effects of a supertranslation on
    asymptotic Bondi data can be understood in three ways.  First, the
    shear and strain are altered by an additive term like ``ð²α/2``
    (or its conjugate).  So, for a supertranslation band-limited at
    ``ℓ_α``, the effects on the shear and strain are bounded by
    ``ℓ_α`` — which will typically be smaller than the data's own band
    limit.  Second, powers of ``ðα`` multiply the various components
    of the Weyl tensor that mix into other components.  This *adds* up
    to ``4ℓ_α`` to the data's band limit — which may be quite
    significant — but only to the components that vanish most quickly.
    Finally, we can think of a supertranslation in terms of the Taylor
    expansion of the data in time.  That is, for a field ``f``, the
    supertranslation shifts it (assuming the Taylor series is a
    reasonable approximation) to

    ```math
    f(t - α) = f(t) - α \dot{f}(t) + \frac{α²}{2} \ddot{f}(t) - \ldots.
    ```

    This is a more complicated effect because progressively higher
    powers are multiplied by progressively higher time derivatives,
    which should have progressively lower amplitudes.  The size of
    this effect is intertwined with details about the waveform.

* **Boosts.**  A boost can similarly be thought of in distinct ways.
    First, the aberration changes the angular shape of the waves,
    which mixes modes.  However, this is typically relatively small
    (governed by the Doppler factor ``\sqrt{(1+β)/(1-β)}``), and
    spectral convergence of the spin-weighted spherical harmonics
    handles this with little overhead.  It turns out that this is
    typically much smaller than the second effect: the spatial
    translation induced by a boost, integrated over time.  For long
    waveforms, even a fairly small boost can accumulate a large
    translation.  This is equivalent to an ``ℓ=1`` supertranslation,
    which can be considered as above, but having a large effect size,
    so that the higher-order terms in the Taylor series are
    significant.

It should be clear that all of these effects are simply too
complicated to give any purely analytical formula.  The reliable
approach is empirical: pad, transform, and check the result with the
diagnostics below.  If your accuracy criteria are not satisfied,
increase ``ℓ_\mathrm{max}`` and repeat.

!!! danger "A separate limit: Interpolation in time"

    The ``ℓ_\mathrm{max}`` value controls the *angular* resolution.
    However, because boosts and supertranslations mix space and time,
    features in time become features in space — and specifically
    angular features.  Therefore, the time interpolation method
    imposes its own, independent accuracy limitations.  This is the
    subject of the "[Spline Errors](@ref)" page, and it is not
    affected by ``ℓ_\mathrm{max}``.

## Diagnostics

The accuracy criteria can usually be specified pretty generically in
terms of the per-``ℓ`` energy

```math
E(ℓ) = \sum_{m=-ℓ}^{ℓ} |f_{ℓm}|²
```

as a function of time.  `Scri` includes a utility function to compute
this quantity:

```julia
using Scri
d = Scri.diagnostics(data, dc)
d′ = Scri.diagnostics(data′, dc)
```

Here, we simply pass the original data and the transformed data, along
with the [`DataComponents`](@ref) descriptor, and the function returns
a dictionary mapping each component symbol (e.g., `:h` or `:ψ₄`) to an
``Nᵗ × (ℓ_\mathrm{max}+1)`` array of energies.

If you have the [`Plots.jl`](https://docs.juliaplots.org/) package
installed, loading it activates a package extension providing a
plotting method: pass the corresponding time array to get a quick
plot of the energy spectrum for each component:

```julia
using Plots  # activates the plotting extension
plots = Scri.diagnostics(t, data, dc)
plots′ = Scri.diagnostics(t′, data′, dc)
```

A fully resolved transform has ``E(ℓ)`` falling to the noise (or
acceptable error) floor at ``ℓ ≤ ℓ_\mathrm{max}``.  If
``E(ℓ_\mathrm{max})`` is comparable to or greater than the energy in
the preceding modes, energy has piled up against the ceiling and the
result is [aliased](https://en.wikipedia.org/wiki/Aliasing) — for
which the cure is simply to increase ``ℓ_\mathrm{max}``.

The [augmented SSHT](@ref "Augmented Direct SSHT") supplies a second,
cheaper check.  Because the pixel grid is larger than the mode basis,
the analysis step also returns the energy ``ξ`` that lives in the null
space — the part of the pixel data that *no* band-limited field can
represent.  Analytically ``ξ ≈ 0``; a large ``ξ`` signals exactly the
same thing as a large ``E(ℓ_\mathrm{max})`` — a transformation extreme
enough to push energy past the band limit.  For ``σ``, ``h``, or
``ψ₄`` data, these are the ``ℓ = 0`` and ``ℓ = 1`` modes — shown, for
example, in red on the [spline-error plots](@ref "Spline Errors").
