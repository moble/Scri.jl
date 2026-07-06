# Tensor Components

For simplicity, we will now assume that the tetrad is regular at null
infinity.  That is, we will interpret the tetrad legs ``(l, m, m̄,
n)`` on this page as the rescaled versions ``(l̃, m̃, m̃̄, ñ)``
defined on [the previous page](@ref "Tetrad").  This means that the
tensor components defined here may generally be nonzero and finite at
null infinity.

The Newman-Penrose Weyl components are defined as

```math
\begin{aligned}
ψ_0 &= c_ψ C_{abcd} l^a m^b l^c m^d, \\
ψ_1 &= c_ψ C_{abcd} l^a n^b l^c m^d, \\
ψ_2 &= c_ψ C_{abcd} l^a m^b \bar{m}^c n^d, \\
ψ_3 &= c_ψ C_{abcd} l^a n^b \bar{m}^c n^d, \\
ψ_4 &= c_ψ C_{abcd} n^a \bar{m}^b n^c \bar{m}^d.
\end{aligned}
```

We can convert between conventions according to the formula

```math
{}^{[A]}ψ_n = c_ψ c_s c_R (c_l c_m)^{2-n}\, ψ_n,
```

The (Maxwell/Faraday) field-strength tensor ``F_{ab}`` is similarly
decomposed into components as

```math
\begin{aligned}
φ_0 &= c_φ F_{ab}\, ℓ^a m^b, \\
φ_1 &= c_φ \tfrac{1}{2} F_{ab}\, (ℓ^a n^b + \bar{m}^a m^b), \\
φ_2 &= c_φ F_{ab}\, \bar{m}^a n^b.
\end{aligned}
```

These are consistent with the Weyl components above: each ``φ_n``
carries spin weight ``1-n`` (just as ``ψ_n`` carries ``2-n``), and
stepping down the tower swaps an ``ℓ``-type dyad slot (``o``) for an
``n``-type one (``ι``).  Paralleling the Weyl conversion ``∝ (c_l
c_m)^{2-n}``, the inter-convention factor is

```math
{}^{[A]}φ_n = c_φ\, (c_l c_m)^{1-n}\, φ_n,
```

with the dyad scaling ``c_l c_m`` raised to the spin weight ``1-n``
and no Riemann-sign factor (the field strength does not see the
curvature convention).

The metric perturbation is defined as

```math
h_{ab} = g_{ab} - \eta_{ab},
```

where ``\eta_{ab}`` is the Minkowski metric.  The strain components
are defined as

```math
\begin{aligned}
h_+ &= \frac{1}{2} (h_{\hat{\theta}\hat{\theta}} - h_{\hat{\phi}\hat{\phi}}), \\
h_\times &= h_{\hat{\theta}\hat{\phi}},
\end{aligned}
```

where the hats indicate orthonormal components in the spherical basis.
These are combined into a single complex strain component[^1]

```math
h = c_h(h_+ - i h_\times).
```

Alternatively, we could define the strain in terms of the tetrad legs, as

```math
h = \frac{c_h}{{c̄}_m^2} h_{ab} m̄^a m̄^b.
```

In either case, we have

```math
{}^{[A]}h = c_h\, h.
```

Near *future* null infinity, we have the asymptotic relation

```math
\Psi_4 \sim -\ddot{h},
```

where the dots indicate time derivatives.

[^1]: Note that defining ``h`` this way forces it to have spin weight
    ``-2``, which is essentially universal in the modern GW astronomy
    community.  There are sources that define the strain to have spin
    weight ``+2``, which is equivalent to taking the complex conjugate
    of our definition.  That choice cannot be accomodated
    automatically by this package.
