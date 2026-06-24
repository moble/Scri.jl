# The standard tetrad

We define the standard tetrad in Minkowski spacetime for simplicity.
The tetrad in a general spacetime will approach this form only
asymptotically, but it is impossible to write down the full tetrad in
closed form.  See [Sachs_1962a](@cite) for the complete construction,
which requires solution of differential equations to find the tetrad.
The transformation laws are simpler to write in Minkowski, and the
asymptotic limit is all we need for the transformations at null
infinity.  Begin with the usual Minkowski coordinates ``(t, r, θ, ϕ)``
and the associated coordinate basis vectors ``(∂ₜ, ∂ᵣ, ∂_θ, ∂_ϕ)``.
The standard null tetrad is then

```math
\begin{aligned}
lᵃ &= \frac{1}{\sqrt{2}} \left(∂ₜ + ∂ᵣ\right)ᵃ, \\
mᵃ &= \frac{1}{r\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
m̄ᵃ &= \frac{1}{r\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &= \frac{1}{\sqrt{2}} \left(∂ₜ - ∂ᵣ\right)ᵃ.
\end{aligned}
```

The inner products are ``l \cdot n = -1`` and ``m \cdot m̄ = 1``, with
all other combinations zero.  In particular, each of these vectors is
null.  This simplicity is what make this tetrad so useful.

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
lᵃ &= \frac{ω²}{\sqrt{2}} \left(-∂_ω\right)ᵃ, \\
mᵃ &= \frac{ω}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
m̄ᵃ &= \frac{ω}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &= \frac{1}{\sqrt{2}} \left(2∂ᵤ + ω² \, ∂_ω\right)ᵃ.
\end{aligned}
```

The inner products are the same as before, but the tetrad is still
singular at ``r = ∞`` (``ω = 0``).  We can rescale to obtain a tetrad
that remains regular at null infinity:

```math
\begin{aligned}
\frac{1}{ω²} lᵃ &\to l̃ᵃ = -\frac{1}{\sqrt{2}} \left(∂_ω\right)ᵃ, \\
\frac{1}{ω} mᵃ &\to m̃ᵃ = \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
\frac{1}{ω} m̄ᵃ &\to m̃̄ᵃ = \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &\to ñᵃ = \sqrt{2} \left(∂ᵤ\right)ᵃ.
\end{aligned}
```

This form is the one we use in a neighborhood of ``ℐ⁺`` in the
conformally compactified spacetime.  Importantly, note that ``m̃ᵃ``,
``m̃̄ᵃ``, and ``ñᵃ`` are all in the tangent space of ``ℐ⁺``, while
``l̃ᵃ`` is transverse to it.  Still, with the conformal metric, the
inner products of the tetrads retain their usual values.

In the ``v`` coordinate system relevant to ``ℐ⁻``, we have

```math
∂ₜ = ∂ᵥ
\qquad
∂ᵣ = ∂ᵥ - ω² \, ∂_ω.
```

The tetrad is

```math
\begin{aligned}
lᵃ &= \frac{1}{\sqrt{2}} \left(2∂ᵥ - ω² \, ∂_ω\right)ᵃ, \\
mᵃ &= \frac{ω}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
m̄ᵃ &= \frac{ω}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &= \frac{ω²}{\sqrt{2}} \left(∂_ω\right)ᵃ.
\end{aligned}
```

Again, we construct a regular tetrad at null infinity by rescaling —
though this time the scaling is opposite for the ``l`` and ``n`` legs:

```math
\begin{aligned}
lᵃ &\to l̃ᵃ = \sqrt{2} \left(∂ᵥ\right)ᵃ, \\
\frac{1}{ω} mᵃ &\to m̃ᵃ = \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
\frac{1}{ω} m̄ᵃ &\to m̃̄ᵃ = \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
\frac{1}{ω²} nᵃ &\to ñᵃ = \frac{1}{\sqrt{2}} \left(∂_ω\right)ᵃ.
\end{aligned}
```

This form is the one we use in a neighborhood of ``ℐ⁻`` in the
conformally compactified spacetime.  Importantly, note that ``l̃ᵃ``,
``m̃ᵃ``, and ``m̃̄ᵃ`` are all in the tangent space of ``ℐ⁻``, while
``ñᵃ`` is transverse to it.
