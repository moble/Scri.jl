# Tensor Components

For simplicity, we will now assume that the tetrad is regular at null
infinity.  That is, we will interpret the tetrad legs ``(l, m, m̄,
n)`` on this page as the rescaled versions ``(l̃, m̃, m̃̄, ñ)``
defined on [the previous page](@ref "Tetrad").  For the Faraday
components below, this is enough to make the components finite — and
generally nonzero — at null infinity.  The Weyl and strain components
require one further power of ``ω``: the curvature must be read as
``K_{abcd} = ω^{-1} C̃_{abcd}`` (where ``C̃`` is the Weyl tensor of
the conformal metric), and the strain as the leading coefficient of
``r\,h``.  We leave those rescalings implicit, because they are shared
by all conventions, and so cannot affect the conversion factors that
are our only concern on this page.

The Newman-Penrose Weyl components are defined as

```math
\begin{aligned}
ψ_0 &= c_ψ C_{abcd} l^a m^b l^c m^d, \\
ψ_1 &= c_ψ C_{abcd} l^a n^b l^c m^d, \\
ψ_2 &= c_ψ C_{abcd} l^a m^b m̄^c n^d, \\
ψ_3 &= c_ψ C_{abcd} l^a n^b m̄^c n^d, \\
ψ_4 &= c_ψ C_{abcd} n^a m̄^b n^c m̄^d.
\end{aligned}
```

Combining the Weyl-tensor conversion ``{}^{[A]}C_{abcd} = c_s c_R
C_{abcd}`` from the [previous pages](@ref "Metric and Curvature") with
the leg scalings ``c_l``, ``1/c_l``, ``c_m``, and ``c̄_m = 1/c_m``, we
can convert between conventions according to the formula

```math
{}^{[A]}ψ_n = c_ψ c_s c_R (c_l c_m)^{2-n}\, ψ_n.
```

The (Maxwell/Faraday) field-strength tensor ``F_{ab}`` is similarly
decomposed into components as

```math
\begin{aligned}
φ_0 &= c_φ F_{ab}\, l^a m^b, \\
φ_1 &= c_φ \tfrac{1}{2} F_{ab}\, (l^a n^b + m̄^a m^b), \\
φ_2 &= c_φ F_{ab}\, m̄^a n^b.
\end{aligned}
```

These are consistent with the Weyl components above: each ``φ_n``
carries spin weight ``1-n`` (just as ``ψ_n`` carries ``2-n``), and
each step down either tower trades a leg scaling as ``c_l`` or ``c_m``
for one scaling as its inverse.  Paralleling the Weyl conversion ``∝
(c_l c_m)^{2-n}``, the inter-convention factor is

```math
{}^{[A]}φ_n = c_φ\, (c_l c_m)^{1-n}\, φ_n,
```

with the dyad scaling ``c_l c_m`` raised to the spin weight ``1-n``,
and neither a signature nor a Riemann-sign factor: no metric or
curvature enters the definition of ``F_{ab}``.

The metric perturbation is defined as

```math
h_{ab} = g_{ab} - η_{ab},
```

where ``η_{ab}`` is the Minkowski metric.  The strain components are
defined as

```math
\begin{aligned}
h_+ &= \frac{1}{2} (h_{θ̂θ̂} - h_{ϕ̂ϕ̂}), \\
h_× &= h_{θ̂ϕ̂},
\end{aligned}
```

where the hats indicate orthonormal components in the spherical basis.
These are combined into a single complex strain component[^1]

```math
h = c_h(h_+ - i h_×).
```

Alternatively, we could define the strain in terms of the tetrad legs,
as

```math
h = \frac{c_h}{c̄_m^2} h_{ab} m̄^a m̄^b;
```

the ``c̄_m^2`` in the denominator exactly cancels the scaling of the
two ``m̄`` legs, so the two definitions agree identically.  In either
case, the only convention factors that survive are ``c_h`` itself and
the signature flip of the metric perturbation, ``{}^{[A]}h_{ab} = c_s
h_{ab}`` (both ``g_{ab}`` and ``η_{ab}`` flip with the signature), so
we have

```math
{}^{[A]}h = c_s c_h\, h.
```

Near *future* null infinity, the radiative parts of these fields are
related.  In the default convention, the relation takes the familiar
form

```math
ψ_4 ∼ -ḧ,
```

where the dots denote derivatives with respect to the retarded time
``u``, and both sides are understood as the finite radiative data
described at the top of this page.  The relation is itself
convention-dependent: dividing the conversion factors above gives

```math
{}^{[A]}ψ_4 ∼ -\frac{c_ψ c_R}{c_h (c_l c_m)^2}\, {}^{[A]}ḧ
```

(the ``c_s`` factors cancel).

[^1]: Note that defining ``h`` this way forces it to have spin weight
    ``-2``, which is essentially universal in the modern GW astronomy
    community.  There are sources that define the strain to have spin
    weight ``+2``, which is equivalent to taking the complex conjugate
    of our definition.  That choice cannot be accommodated
    automatically by this package.

Finally, we also define the shear on ``ℐ⁺`` as

```math
σ = -c_σ m^a m^b ∇_a l_b,
```

which has spin weight ``+2``.  The conversion factor follows the same
reasoning as above: the two ``m`` legs are vectors, contributing
``c_m^2``; the covariant derivative is convention-independent (the
Christoffel symbols do not change, as shown under "[Metric and
Curvature](@ref)"); and ``l_b = g_{bc} l^c`` is the *lowered* tetrad
vector, contributing ``c_s c_l`` — one power of the metric survives.
So

```math
{}^{[A]}σ = c_s c_σ c_l c_m^2\, σ.
```

On ``ℐ⁺``, the radiative degrees of freedom lie along the transverse
``l`` congruence.  On ``ℐ⁻`` the transverse congruence is instead
``n``, and the relevant shear is a *different* Newman–Penrose spin
coefficient, ``λ``:

```math
λ = c_λ m̄^a m̄^b ∇_a n_b,
```

which has spin weight ``-2``.  This is exactly analogous to the Weyl
tower: just as the radiative Weyl component switches from the outgoing
``ψ₄`` on ``ℐ⁺`` to the incoming ``ψ₀`` on ``ℐ⁻``, the radiative shear
switches from ``σ`` to ``λ``.  The conversion factor for ``λ`` follows
the same logic:

```math
{}^{[A]}λ = \frac{c_s c_λ}{c_l c_m^2}\, λ.
```
