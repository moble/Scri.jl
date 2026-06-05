```@meta
CurrentModule = Scri
```

# Scri

!!! warning "Documentation in progress"

    This documentation is still being written.  Some sections may be
    incomplete.  Please bear with me while I fill things in.

This is a Julia package for working with gravitational and
electromagnetic waveforms at future or past null infinity, including
transformations under the BMS group.  The package is designed to be
fast, accurate, and easy to use, and is intended to be a useful tool
for researchers in numerical relativity and gravitational-wave
astronomy.

The main functionality of the package is provided by the
[`transform!`](@ref) function, which takes as input a set of waveforms
represented as a three-dimensional array of complex numbers, the
corresponding times, and a set of parameters specifying the BMS
transformation to be applied, and returns the transformed waveform and
the corresponding retarded times in the new frame.

```julia
using Quaternionic
using Scri

# Construct random test data
ℓₘₐₓ = 8
data_components = ("h", "Psi4")
Nᵗ = 10_000
Nᵐ = (ℓₘₐₓ + 1)^2
Nᵈ = length(data_components)
data = randn(ComplexF64, Nᵐ, Nᵗ, Nᵈ)
t = collect(LinRange(-50, 5000, Nᵗ))

# Construct a random-ish BMS transformation
v⃗ = 1e-3 * normalize(randn(QuatVecF64))
R = randn(RotorF64)
α = 1e-3 * randn(ComplexF64, Nᵐ)

# Perform the transformation
data′, t′ = Scri.transform!(data, t, v⃗, R, α; data_components)
```

Note a few very important points:

  1. The `transform!` function modifies its input `data` in place,
      and also returns the modified version.  (The exclamation mark is
      the Julia convention for indicating that the function modifies
      its arguments in place.)  This is done for performance reasons,
      to avoid unnecessary allocations.  If you want to keep the
      original data, you can make a copy before calling `transform!`
      with the [`Base.copy`](@extref Julia) function.
  2. The input `data` must be shaped as a three-dimensional array,
      with the first dimension corresponding to the mode weights, the
      second dimension corresponding to time, and the third dimension
      corresponding to the different field components.
  3. The first index of the `data` array represents mode weights
      starting with ``ℓ=0``, even for fields with nonzero spin weight.
      The mode weights are ordered in the standard way, with the
      ``m=-ℓ`` mode first and the ``m=+ℓ`` mode last, then
      incrementing ``ℓ``.
  4. The third index of the `data` array must correspond to the
      `data_components` argument.  That is, it must have the same
      length, and describe the components in the same order.
  5. In this example, the supertranslation `α` was constructed
      randomly, and doesn't represent a real-valued function.
      Internally, `transform!` automatically imposes the reality
      condition on `α` by averaging each mode with its
      complex-conjugate partner.
