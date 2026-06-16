# The standard tetrad

We define the standard tetrad in Minkowski spacetime for simplicity.
The tetrad in a general spacetime will approach this form only
asymptotically, but it is impossible to write down the full tetrad in
closed form.  See [Sachs_1962a](@cite) for the complete construction.
The transformation laws are simpler to write in Minkowski, and the
asymptotic limit is all we need for the transformations at null
infinity.  Begin with the usual Minkowski coordinates ``(t, r, θ, φ)``
and the associated coordinate basis vectors ``(∂ₜ, ∂ᵣ, ∂_θ, ∂_φ)``.
The standard null tetrad is then

```math
\begin{aligned}
l^a &= \frac{1}{\sqrt{2}} \left(∂ₜ + ∂ᵣ\right)ᵃ, \\
m^a &= \frac{1}{r\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
\bar{m}^a &= \frac{1}{r\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
n^a &= \frac{1}{\sqrt{2}} \left(∂ₜ - ∂ᵣ\right)ᵃ.
\end{aligned}
```

Of course, ``t`` and ``r`` are poor choices for coordinates near null
infinity, so we transform to systems using the inverse radial
coordinate ``ω = 1/r`` and retarded or advanced time coordinates ``u =
t - r`` or ``v = t + r``.

In the ``u`` coordinate system relevant to ``ℐ⁺``, we have

```math
∂ₜ = ∂ᵤ
\qquad
∂ᵣ = -∂ᵤ - ω² \, ∂_ω,
```

so the tetrad is

```math
\begin{aligned}
l^a &= \frac{1}{\sqrt{2}} \left(-ω² \, ∂_ω\right)ᵃ, \\
m^a &= \frac{ω}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
\bar{m}^a &= \frac{ω}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
n^a &= \frac{1}{\sqrt{2}} \left(2∂ᵤ + ω² \, ∂_ω\right)ᵃ.
\end{aligned}
```

Note that the tetrad is singular at ``r = ∞`` (``ω = 0``), so we
rescale to obtain a tetrad that remains regular at null infinity:

```math
\begin{aligned}
\frac{1}{ω²} l^a &\to l̃^a = -\frac{1}{\sqrt{2}} \left(∂_ω\right)ᵃ, \\
\frac{1}{ω} m^a &\to m̃^a = \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
\frac{1}{ω} \bar{m}^a &\to \bar{m̃}^a = \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
n^a &\to ñ^a = \sqrt{2} \left(∂ᵤ\right)ᵃ.
\end{aligned}
```

In the ``v`` coordinate system relevant to ``ℐ⁻``, we have

```math
∂ₜ = ∂ᵥ
\qquad
∂ᵣ = ∂ᵥ - ω² \, ∂_ω.
```

The tetrad is

```math
\begin{aligned}
l^a &= \frac{1}{\sqrt{2}} \left(2∂ᵥ - ω² \, ∂_ω\right)ᵃ, \\
m^a &= \frac{ω}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
\bar{m}^a &= \frac{ω}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
n^a &= \frac{1}{\sqrt{2}} \left(ω² \, ∂_ω\right)ᵃ.
\end{aligned}
```

Again, we construct a regular tetrad at null infinity by rescaling —
though this time the scaling is opposite for the ``l`` and ``n`` legs:

```math
\begin{aligned}
l^a &\to l̃^a = \sqrt{2} \left(∂ᵥ\right)ᵃ, \\
\frac{1}{ω} m^a &\to m̃^a = \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
\frac{1}{ω} \bar{m}^a &\to \bar{m̃}^a = \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_φ\right)ᵃ, \\
\frac{1}{ω²} n^a &\to ñ^a = \frac{1}{\sqrt{2}} \left(∂_ω\right)ᵃ.
\end{aligned}
```
