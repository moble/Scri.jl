```@meta
CurrentModule = Scri
```

# Scri

This is a Julia package for working with gravitational and
electromagnetic waveforms at future or past null infinity, including
transformations under the BMS group.  The package is designed to be
fast, accurate, and easy to use, and is intended to be a useful tool
for researchers in numerical relativity and gravitational-wave
astronomy.

## Installation

From the Julia REPL, press `]` to enter the package manager, then run:

```julia
pkg> add Scri
```

Or equivalently, from any Julia session:

```julia
using Pkg
Pkg.add("Scri")
```

## Example Usage

The main functionality of the package is provided by the
[`transform!`](@ref) function, which takes as input a set of waveforms
represented as a three-dimensional array of complex numbers, the
corresponding times, and a set of parameters specifying the BMS
transformation to be applied, and returns the transformed waveform and
the corresponding retarded times in the new frame.

```@example quickstart
using Random  # hide
Random.seed!(1234)  # hide
using Quaternionic
using Scri

# You have to tell the transform what components are in the data array
data_components = ("h", "Psi4")

# Construct random test data
ℓₘₐₓ = 8
Nᵗ = 10_000
Nᵐ = (ℓₘₐₓ + 1)^2
Nᵈ = length(data_components)
data = randn(ComplexF64, Nᵐ, Nᵗ, Nᵈ)
t = collect(LinRange(-50.0, 5000.0, Nᵗ))

# Construct a random-ish BMS transformation
v⃗ = 1e-3 * normalize(randn(QuatVecF64))
R = randn(RotorF64)
α = 1e-3 * randn(ComplexF64, Nᵐ)

# Perform the transformation
data′, t′ = transform!(data, t, v⃗, R, α; data_components)
size(data′), extrema(t′)
```

A few very important points:

  1. `transform!` modifies its input `data` **in place** (and also
      returns it); [copy](@extref Julia Base.copy) first if you need
      the original.
  2. The `data` array is three-dimensional — modes × times ×
      components — with the mode weights ordered by increasing ``ℓ``
      from 0 (even for nonzero spin weight), then increasing ``m``,
      and the third axis matching `data_components` in length and
      order.
  3. The supertranslation `α` above was constructed randomly, so it
      doesn't represent a real-valued function; `transform!`
      automatically imposes the reality condition by averaging each
      mode with its complex-conjugate partner.

See [Transforming Waveforms](@ref) for the full description of the
data layout, the calling forms, and the returned time grid.

## Performance

The transformation is heavily multithreaded.  Start Julia with threads
enabled — e.g., `julia -t auto` — or `transform!` will run on a single
core.  The dominant costs scale steeply with the angular band limit,
so read [Choosing ``ℓ_\mathrm{max}``](@ref) before padding your data
to very high ``ℓ``.

## Name and Pronunciation

"Scri" (rhymes with "sky") is the pronunciation of the character ``ℐ``
(script I; given by the latex `\mathscr{I}`), which is the notation
for [null infinity](https://en.wikipedia.org/wiki/Null_infinity) in
the theory of asymptotically flat spacetimes, introduced by [Penrose
in 1963](@cite Penrose_1963) and now standard throughout the
literature.  We usually distinguish between future null infinity,
denoted ``ℐ⁺``, and past null infinity, denoted ``ℐ⁻``.  This package
defaults to working with future null infinity, since the author's
primary interest is in the emission of gravitational waves from
isolated systems, but there are options to model past null
infinity, to investigate the detection of incoming radiation.
