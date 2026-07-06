# Tensor Components

The Newman-Penrose Weyl components are defined as

```math
\begin{aligned}
\Psi_0 &= C_{abcd} \ell^a m^b \ell^c m^d, \\
\Psi_1 &= C_{abcd} \ell^a n^b \ell^c m^d, \\
\Psi_2 &= C_{abcd} \ell^a m^b \bar{m}^c n^d, \\
\Psi_3 &= C_{abcd} \ell^a n^b \bar{m}^c n^d, \\
\Psi_4 &= C_{abcd} n^a \bar{m}^b n^c \bar{m}^d.
\end{aligned}
```

The (Maxwell/Faraday) field-strength tensor ``F_{ab}`` is similarly
decomposed into components as

```math
\begin{aligned}
φ_0 &= F_{ab}\, ℓ^a m^b, \\
φ_1 &= \tfrac{1}{2} F_{ab}\, (ℓ^a n^b + \bar{m}^a m^b), \\
φ_2 &= F_{ab}\, \bar{m}^a n^b.
\end{aligned}
```

These are consistent with the Weyl components above: each ``φ_n``
carries spin weight ``1-n`` (just as ``Ψ_n`` carries ``2-n``), and
stepping down the tower swaps an ``ℓ``-type dyad slot (``o``) for an
``n``-type one (``ι``).  The overall sign is the ``c_φ`` of the
[convention table](@ref "Convention parameters") below, and —
paralleling the Weyl conversion ``∝ (c_l c_m)^{2-n}`` — the
inter-convention factor is

```math
φ_n^{[X]} = c_φ\, (c_l c_m)^{1-n}\, φ_n^{[\mathrm{SpEC}]},
```

with the dyad scaling ``c_l c_m`` raised to the spin weight
``1-n`` and no Riemann-sign factor (the field strength does not see the
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
h_\times &= h_{\hat{\theta}\hat{\phi}}, \\
h &= h_+ - i h_\times.
\end{aligned}
```

where the hats indicate orthonormal components in the spherical basis.
Near *future* null infinity, we have the asymptotic relation

```math
\Psi_4 \sim -\ddot{h},
```

where the dots indicate time derivatives.
