# Tetrad

We define the standard tetrad in Minkowski spacetime for simplicity.
The tetrad in a general spacetime will approach this form only
asymptotically, but it is impossible to write down the full tetrad in
closed form.  See [Sachs_1962a](@cite) for the complete construction,
which requires solution of differential equations to find the tetrad.
The transformation laws are simpler to write in Minkowski, and the
asymptotic limit is all we need for the transformations at null
infinity.

Begin with the usual Minkowski coordinates ``(t, r, θ, ϕ)`` and the
associated coordinate basis vectors ``(∂ₜ, ∂ᵣ, ∂_θ, ∂_ϕ)``.  The
standard null tetrad is then[^1]

```math
\begin{aligned}
lᵃ &= \frac{c_l}{\sqrt{2}} \left(∂ₜ + ∂ᵣ\right)ᵃ, \\
mᵃ &= \frac{c_m}{r\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
m̄ᵃ &= \frac{\bar{c}_m}{r\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &= \frac{c_l^{-1}}{\sqrt{2}} \left(∂ₜ - ∂ᵣ\right)ᵃ.
\end{aligned}
```

(Note that we require ``|c_m| = 1``.)  The inner products are ``l
\cdot n = -c_s`` and ``m \cdot m̄ = c_s``, with all other combinations
zero.  In particular, each of these vectors is null.  This simplicity
is what makes this tetrad so useful.

[^1]: The asymptotic-transformation paper [Boyle_2015](@citet) works
    instead with the *co*-tetrad defined as ``l_a = (du)_a``, ``m_a =
    -\frac{\sqrt{2}}{1+zz̄}(dz̄)_a``, ``n_a = -(dω)_a``, where ``z =
    e^{iϕ} \cot(θ/2)`` is the stereographic coordinate.  That paper is
    unfortunately not explicit about its signature choice, but is
    consistent with the opposite signature ``{+}{-}{-}{-}``.  Lowering
    the indices of the regular ``ℐ⁺`` tetrad below with the conformal
    metric of that signature shows that this cotetrad is the member of
    our family with ``c_s = -1`` and ``c_l = \sqrt{2}``: the ``l`` and
    ``n`` normalizations then match exactly, with no residual sign.
    (The extensions of ``n`` off ``ℐ⁺`` differ at ``O(ω²)``, which is
    immaterial at ``ℐ⁺`` itself.)  The ``m`` legs agree up to the
    standard phase relating the polar spin frame used here to the
    stereographic one, together with an overall sign, which we book as
    ``c_m = -1``.  The resulting dyad scaling ``c_l c_m = -\sqrt{2}``,
    together with the Newman–Penrose-vs-GHP eth (see "[The Operator
    ð](@ref the-operator-eth)"), accounts for the factor of
    ``\sqrt{2}`` and the sign between that paper's component-mixing
    parameter and the one derived in "[BMS action on fields](@ref)".

The conventions found in the literature are more widely varied than we
have accounted for here.  Beyond simple name differences[^2], some
authors use different normalizations, so that ``l \cdot n = -2``, for
example.  Here, we impose our simpler (and far more common)
normalization to reduce the number of conventions; with this choice, a
hypothetical ``c_n`` is restricted to equal ``1/c_l``, and
``c_{\bar{m}}`` is restricted to equal ``\bar{c}_m``.

[^2]:  [Sachs_1962a](@Citet), for example, uses entirely different
    names ``(k, t, \bar{t}, m)`` for what we would call ``(l, m,
    \bar{m}, n)``.  Numerous authors also swap the names of ``l`` and
    ``n`` or of ``m`` and ``\bar{m}``.

## Asymptotic coordinates

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
lᵃ &= c_l \frac{ω²}{\sqrt{2}} \left(-∂_ω\right)ᵃ, \\
mᵃ &= c_m \frac{ω}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
m̄ᵃ &= \bar{c}_m \frac{ω}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &= \frac{1}{c_l} \frac{1}{\sqrt{2}} \left(2∂ᵤ + ω² \, ∂_ω\right)ᵃ.
\end{aligned}
```

The inner products are the same as before, but the tetrad is still
singular at ``r = ∞`` (``ω = 0``).  We can rescale to obtain a tetrad
that remains regular at null infinity; evaluating at ``ℐ⁺`` (``ω =
0``), the result is

```math
\begin{aligned}
\frac{1}{ω²} lᵃ &\to l̃ᵃ = -c_l \frac{1}{\sqrt{2}} \left(∂_ω\right)ᵃ, \\
\frac{1}{ω} mᵃ &\to m̃ᵃ = c_m \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
\frac{1}{ω} m̄ᵃ &\to m̃̄ᵃ = \bar{c}_m \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &\to ñᵃ = \frac{1}{c_l} \sqrt{2} \left(∂ᵤ\right)ᵃ.
\end{aligned}
```

This form is the one we use at ``ℐ⁺`` in the conformally compactified
spacetime.  (The ``n`` leg needed no rescaling; we have simply dropped
its ``ω² ∂_ω`` term, which vanishes at ``ℐ⁺``.) Importantly, note that
``m̃ᵃ``, ``m̃̄ᵃ``, and ``ñᵃ`` are all in the tangent space of ``ℐ⁺``,
while ``l̃ᵃ`` is transverse to it.  Still, with the *conformal*
metric, the inner products of the tetrad at ``ℐ⁺`` retain their usual,
convenient values.

In the ``v`` coordinate system relevant to ``ℐ⁻``, we have

```math
∂ₜ = ∂ᵥ
\qquad
∂ᵣ = ∂ᵥ - ω² \, ∂_ω.
```

The tetrad is

```math
\begin{aligned}
lᵃ &= c_l \frac{1}{\sqrt{2}} \left(2∂ᵥ - ω² \, ∂_ω\right)ᵃ, \\
mᵃ &= c_m \frac{ω}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
m̄ᵃ &= \bar{c}_m \frac{ω}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
nᵃ &= \frac{1}{c_l} \frac{ω²}{\sqrt{2}} \left(∂_ω\right)ᵃ.
\end{aligned}
```

Again, we construct a regular tetrad at null infinity by rescaling —
though this time the scaling is opposite for the ``l`` and ``n`` legs
— and evaluating at ``ℐ⁻`` (``ω = 0``):

```math
\begin{aligned}
lᵃ &\to l̃ᵃ = c_l \sqrt{2} \left(∂ᵥ\right)ᵃ, \\
\frac{1}{ω} mᵃ &\to m̃ᵃ = c_m \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
\frac{1}{ω} m̄ᵃ &\to m̃̄ᵃ = \bar{c}_m \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right)ᵃ, \\
\frac{1}{ω²} nᵃ &\to ñᵃ = \frac{1}{c_l} \frac{1}{\sqrt{2}} \left(∂_ω\right)ᵃ.
\end{aligned}
```

This form is the one we use at ``ℐ⁻`` in the conformally compactified
spacetime.  (Here it is the ``l`` leg that needed no rescaling, and
whose ``ω² ∂_ω`` term vanishes at ``ℐ⁻``.) Importantly, note that
``l̃ᵃ``, ``m̃ᵃ``, and ``m̃̄ᵃ`` are all in the tangent space of ``ℐ⁻``,
while ``ñᵃ`` is transverse to it.
