# [Overview](@id conventions-overview)

## Motivation

As with all of GR, the literature on spacetime asymptotics is,
unfortunately, full of subtly varying conventions for objects that are
fundamentally the same.  The metric signature, the tetrad
normalization, curvature quantities, and even the ``ð`` operator all
differ between references.  These choices are not just internal; they
actually affect the transformation laws of the fields.  For example,
the Newman–Penrose Weyl components ``ψ_n`` are defined as contractions
of the Weyl tensor with various elements of the tetrad, ``(l, m,
\bar{m}, n)``.  But these are null vectors, so the normalizations are
not fixed and naturally (even accounting for different letters used to
represent the elements) different authors choose different
normalizations for the tetrad elements.  For example, two reasonable
choices found in the literature require

```math
\begin{aligned}
l \qquad &↔ \qquad l\sqrt{2}, \\
\hphantom{\sqrt{2}}n \qquad &↔ \qquad n/\sqrt{2},
\end{aligned}
```

which changes the transformation law for ``ψ₃`` in a nontrivial way:

```math
ψ₃' = \frac{e^{-iλ}}{κ³} \left[ψ₃ + \frac{ðα}{2κ} ψ₄\right]
\qquad ↔ \qquad
ψ₃' = \frac{e^{-iλ}}{κ³} \left[ψ₃ + \frac{ðα}{2κ} \frac{ψ₄}{\sqrt{2}}\right].
```

Assuming conventions that differ from those used to construct the
input data can lead to incorrect results, so it is important to know
which conventions are being used.  This package provides a
`Conventions` type to specify the conventions of the input data, so
that the correct transformation laws are applied.

## [The `Conventions` Type](@id conventions-type)

To keep track of the conventions in effect, we provide a simple
`Conventions` type, which contains fields that specify the conventions
a user assumes, compared to the default conventions used natively by
this package.  This type frequently allows the conventions to be
compiled away, so that the overhead of using non-default conventions
is very small.  (We assume units where ``G = c = 1``, if relevant.)

The conventions are fixed by parameters found in each of the key
definitions in the following pages, as summarized in the table below.
The definitions are meant to be interpreted lexically, so that any
other author would write down any one of those same formulas entirely
within their own conventions, but might differ over the value of the
``c_X`` parameter in that formula.

| Parameter                              | Quantity                  | Default | Space         |
|:---------------------------------------|:--------------------------|:-------:|:-------------:|
| [``c_s``](@ref "Metric and Curvature") | Metric signature          | ``1``   | ``\{-1,+1\}`` |
| [``c_R``](@ref "Metric and Curvature") | Riemann definition        | ``1``   | ``\{-1,+1\}`` |
| [``c_l``](@ref "Tetrad")               | Tetrad ``l`` scale        | ``1``   | ``ℝ^×``       |
| [``c_m``](@ref "Tetrad")               | Tetrad ``m`` phase        | ``1``   | ``e^{iℝ}``    |
| [``c_ψ``](@ref "Tensor Components")    | Weyl component definition | ``1``   | ``ℂ^×``       |
| [``c_σ``](@ref "Tensor Components")    | Shear definition          | ``1``   | ``ℂ^×``       |
| [``c_h``](@ref "Tensor Components")    | Complex strain definition | ``1``   | ``ℂ^×``       |
| [``c_φ``](@ref "Tensor Components")    | Faraday definition        | ``1``   | ``ℂ^×``       |
| [``c_ð``](@ref the-operator-eth)       | ``ð`` definition          | ``1``   | ``ℂ^×``       |
| [``c_α``](@ref "BMS Transformations")  | Supertranslation sign     | ``1``   | ``\{-1,+1\}`` |

The "Space" column indicates the set of values that each parameter can
take.  The superscript ``×`` indicates the multiplicative group of
nonzero elements of the field — meaning that the parameters cannot be
zero.

```@docs; canonical=false
Scri.Conventions
```
