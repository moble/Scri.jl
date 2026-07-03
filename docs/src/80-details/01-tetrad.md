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
null.  This simplicity is what make this tetrad so useful.[^boyle2015]

[^boyle2015]:
    The asymptotic-transformation paper [Boyle_2015](@citet) works
    instead with the *co*-tetrad defined as ``l_a = (du)_a``, ``m_a =
    -\frac{\sqrt{2}}{1+z\bar{z}}(d\bar{z})_a``, ``n_a = -(dω)_a``,
    where ``z`` is the stereographic coordinate.  These are one-forms
    rather than vectors, and differ from the tetrad here by index
    placement, an overall sign on ``n``, and a sign and ``\sqrt{2}``
    in ``m`` (the magnitude of the normalization ``m ⋅ m̄`` still
    matches ours, up to the signature-dependent sign).  That paper is
    unfortunately not explicit about its signature choice, but is
    consistent with the opposite signature ``{+}{-}{-}{-}``.  Together
    with the Newman–Penrose-vs-GHP eth (see "[The eth
    operator](@ref)"), these convention differences account for the
    factor of ``\sqrt{2}`` and the sign between that paper's
    component-mixing parameter and the one derived in "[BMS action on
    fields](@ref)".

## Convention parameters on the tetrad

The tetrad above is the package's *default*, but various sources in
the literature use different scaling and phases for the legs.  Three
of the [convention parameters](@ref "Convention parameters") act
directly on the tetrad — ``c_l`` (the scale of the ``l`` leg), ``c_m``
(the scale of the ``m`` leg), and ``c_s`` (the signature) — giving
the general null tetrad

```math
\begin{aligned}
ℓ^a &= \frac{c_l}{\sqrt{2}}\left(∂ₜ + ∂ᵣ\right)^a, &\qquad
m^a &= \frac{c_m}{r\sqrt{2}}\left(∂_θ + \frac{i}{\sin θ}∂_ϕ\right)^a, \\
n^a &= \frac{1}{c_l\sqrt{2}}\left(∂ₜ - ∂ᵣ\right)^a, &
\bar m^a &= \frac{\bar{c_m}}{r\sqrt{2}}\left(∂_θ - \frac{i}{\sin θ}∂_ϕ\right)^a,
\end{aligned}
```

with inner products ``ℓ ⋅ n = -c_s`` and ``m ⋅ \bar m = c_s`` (all
others zero).  The ``n`` leg include ``1/c_l`` precisely so that ``ℓ ⋅
n`` is preserved under the rescaling, and ``c_s = +1`` is the ``-+++``
signature (the inner products take their familiar values ``-1`` and
``+1``).  The defaults ``c_l = 1``, ``c_m = 0``, ``c_s = +1`` recover
the standard tetrad above, which is what `Scri.jl` uses by default; we
write the explicit forms with the defaults below, reinstating ``c_l,
c_m, c_s`` only where they survive into a transformation law (see
["BMS action on the tetrad"](@ref "BMS action on the tetrad") and
["BMS action on fields"](@ref "BMS action on fields")).  These three,
together with the eth coefficient ``c_ð``, are the only convention
parameters that touch the tetrad and its derivatives; the
curvature-sign parameters ``c_R, c_Ψ, c_σ, c_h, c_φ`` enter only when
the field *components* are assembled.

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
