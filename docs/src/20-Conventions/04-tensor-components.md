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

The two shears defined at the end of this page need no further power
of ``ω``, but they do need a connection, because they are coefficients
of one and not components of a tensor.  Written with the rescaled legs
and with the covariant derivative of the conformal metric ``ĝ``, both
are finite at null infinity, and both agree with the shears computed
from the physical metric and its own legs: the terms by which the two
connections differ are each proportional to ``m̃ ⋅ m̃`` or to the
inner product of ``m̃`` with the leg being differentiated, and both of
those vanish, so nothing survives the contraction.  The quantity being
transformed is therefore the familiar peeling coefficient,
``\lim_{r → ∞} r^2 σ`` at ``ℐ⁺`` and ``\lim_{r → ∞} r^2 λ`` at
``ℐ⁻``.

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

These are consistent with the Weyl components above: each ``φ_n`` has
spin weight ``1-n`` (just as ``ψ_n`` has ``2-n``), and each step down
either tower trades a leg scaling as ``c_l`` or ``c_m`` for one
scaling as its inverse.  Paralleling the Weyl conversion ``∝ (c_l
c_m)^{2-n}``, the inter-convention factor is

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

We assume Bondi gauge near null infinity, and the angular part of the
metric is

```math
r^2 h_{AB} = r^2 \left(γ_{AB} + \frac{C_{AB}}{r} + …\right).
```

Here, ``r`` is the areal radius, ``γ_{AB}`` is the unit round sphere,
and ``C_{AB}`` the trace-free "shear tensor" (its time-derivative
being the "News tensor", though each is a tensor in a very limited
sense).  The radiative strain is then the limit implied at the top of
this page,

```math
h = c_h \lim_{r → ∞} r\,(h_+ - i h_×)
  = \frac{c_h}{2}\, C_{AB}\, q̄^A q̄^B,
```

where ``q^A = \sqrt{2}\, m̃^A / c_m`` is the dyad of the previous page
with its convention factor divided out.  The dyad appears twice here,
so its intrinsic sign is irrelevant.

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
    ``-2``, which is *nearly* universal in the modern GW astronomy
    community.  There are sources that define the strain to have spin
    weight ``+2``, which is equivalent to taking the complex conjugate
    of our definition.  That choice cannot be accommodated
    automatically by this package.

Finally, we also define the shear on ``ℐ⁺`` as

```math
σ = -c_σ m^a m^b ∇_a l_b,
```

which has spin weight ``+2``.  With the default ``c_σ = 1`` this is
Eq. (85c) of [Moxon_2020](@citet), the definition used by the SpECTRE
CCE code that produces SXS shear data.  The conversion factor follows
the same reasoning as above: the two ``m`` legs are vectors,
contributing ``c_m^2``; the covariant derivative is
convention-independent (the Christoffel symbols do not change, as
shown under "[Metric and Curvature](@ref)"); and ``l_b = g_{bc} l^c``
is the *lowered* tetrad vector, contributing ``c_s c_l`` — one power
of the metric survives.  So

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

In the default convention, both radiative shears are tied to the
strain.  The relation follows from the definitions above.  At ``ℐ⁺``
the conformal angular metric is ``h_{AB} = γ_{AB} + ω\,C_{AB} + …``
and the cross term of the conformal metric is ``ĝ_{uω} = 1``, so the
only Christoffel symbol needed is

```math
Γ̂^u{}_{AB} = -\frac{1}{2} ∂_ω h_{AB} = -\frac{C_{AB}}{2}.
```

Contracting it with the lowered leg ``l̃_a = -c_s c_l (du)_a/\sqrt{2}``
of the previous page, and with the dyad twice, gives

```math
σ = -c_σ\, m̃^a m̃^b ∇̂_a l̃_b
  = -c_σ\, \frac{c_s c_l}{\sqrt{2}}\, m̃^A m̃^B Γ̂^u{}_{AB}
  = c_s c_σ c_l c_m^2\, \frac{J^{(1)}}{2\sqrt{2}},
\qquad
J^{(1)} = \frac{1}{2} C_{AB}\, q^A q^B,
```

the prefactor being exactly the conversion factor derived above, as it
must be.  Here ``J^{(1)}`` is the leading coefficient of the Bondi
metric function ``J`` of [Moxon_2020](@citet), and it equals ``h̄`` at
the default ``c_h = 1`` by the strain expression given earlier.  The
same computation at ``ℐ⁻``, with ``ñ`` in place of ``l̃``, gives the
conjugate relation:

```math
σ = \frac{h̄}{2\sqrt2} \quad (ℐ⁺),
\qquad
λ = \frac{h}{2\sqrt2} \quad (ℐ⁻),
```

with all three understood as the finite radiative data described at
the top of this page.
This agrees with Eq. (86c) of [Moxon_2020](@citet) in the limit.  That
equation gives the shear at finite radius, in terms of the same ``J``
and of ``K = \sqrt{1 + J J̄}``, as

```math
σ = -\frac{(1 + K)\, ∂_r J}{4 \sqrt{2}\, K}
    + \frac{J^2\, ∂_r J̄}{4 \sqrt{2}\, K (1 + K)}.
```

As ``r → ∞`` we have ``J = J^{(1)}/r + …``, so ``K → 1`` and the
second term falls off like ``r^{-4}``, leaving only the first.  There
``∂_r J = -J^{(1)}/r^2 + …``, whose sign cancels the one in front, and
``(1 + K)/K → 2`` turns the ``1/4\sqrt{2}`` into ``1/2\sqrt{2}``:

```math
σ = \frac{J^{(1)}}{2 \sqrt{2}\, r^2} + O(r^{-3}),
```

which is the rescaled value above.

!!! warning "The dyad is null only with respect to its own metric"
    The round dyad used above is null at ``ℐ``, where the conformal
    angular metric is the unit sphere, but it is *not* null on the
    perturbed metric at finite radius, where ``m ⋅ m = c_m^2
    J^{(1)}/r``.  Repeating the computation there — with the physical
    metric, the physical connection, and ``r^2 h_{AB}`` in place of
    ``h_{AB}`` — therefore requires the correction that restores
    nullity, which is written out in Eq. (81a) of
    [Moxon_2020](@citet) and reads ``m^A = c_m(q^A - J q̄^A/2)/
    \sqrt{2} r`` to first order in ``J``.  That correction is not a
    small one: it contributes ``-4J^{(1)}`` against the ``+2J^{(1)}``
    of the term that survives without it, so omitting it reverses the
    sign of the relation between the shear and the strain.
