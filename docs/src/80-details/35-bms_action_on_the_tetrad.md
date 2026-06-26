# BMS action on the tetrad

The next question is how [the standard tetrad](@ref "The standard
tetrad") transforms under coordinate transformations — and
specifically BMS transformations.

## Future null infinity

One *very* important point is that the tetrad is *defined* in terms of
the coordinates.  In particular, if we have an unprimed coordinate
system and a primed coordinate system, we actually have *two distinct*
tetrads.  For example, we have these tetrads that are regular near
``ℐ⁺``:

```math
\begin{aligned}
l̃ &= -\frac{1}{\sqrt{2}} ∂_ω &\qquad\qquad
l̃' &= -\frac{1}{\sqrt{2}} ∂_{ω'} \\
m̃ &= \frac{1}{\sqrt{2}} \left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right) &
m̃' &= \frac{1}{\sqrt{2}} \left(∂_{θ'} + \frac{i}{\sin θ'} ∂_{ϕ'}\right) \\
m̄̃ &= \frac{1}{\sqrt{2}} \left(∂_θ - \frac{i}{\sin θ} ∂_ϕ\right) &
m̄̃' &= \frac{1}{\sqrt{2}} \left(∂_{θ'} - \frac{i}{\sin θ'} ∂_{ϕ'}\right) \\
ñ &= \sqrt{2} ∂ᵤ &
ñ' &= \sqrt{2} ∂_{u'}
\end{aligned}
```

The crucial point is that even though the primes on the coordinates
indicate a different coordinate system, the primes on the tetrads also
indicate *fundamentally different vectors*, not just the coordinate
transformation of a single tetrad.

Each tetrad is orthonormal in its *own* conformal frame: ``l̃ ⋅ ñ =
-1`` and ``m̃ ⋅ m̄̃ = 1``, with all other combinations zero, in the
unprimed frame — and the same values for the primed tetrad in the
*primed* conformal frame.  This distinction is the source of the
conformal factors below: the vanishing products are conformally
invariant, but the two normalizations are not — each holds only in the
frame where its tetrad is the regular one, so they are exactly the
frame-dependent data that a BMS transformation rescales.

It is cleanest to track this through the conformal factor itself.  The
unprimed and primed frames complete the *same* physical metric ``g``
to conformal metrics ``ĝ = ω²g`` and ``ĝ' = ω'²g``, and the BMS
conformal factor relates the two radial coordinates by ``ω' = κω``.
The conformal metrics therefore scale as

```math
ĝ' = κ² \, ĝ.
```

A null tetrad orthonormal in ``ĝ'`` is thus just ``1/κ`` times one
orthonormal in ``ĝ``, so *every* leg — the generator ``ñ``, the dyad
``m̃``, ``m̄̃``, and the transverse ``l̃`` alike — picks up the same
weight ``1/κ``.  The results ``ñ' = ñ/κ``, etc., derived below all
follow from this one weight.  (The induced metric on ``ℐ⁺`` itself is
degenerate, with ``ñ`` its degenerate direction and the *cut metric*
on the screen ``m̃``, ``m̄̃`` its non-degenerate part; ``ω' = κω``
scales that cut metric by ``κ²`` and the affine parameter ``u`` along
the generators by ``κ`` — the same ``κ``, which largely defines the
BMS group as opposed to an arbitrary conformal map of the sphere.)

This uniform ``1/κ`` is *not* related to the differing powers ``ω⁻²``,
``ω⁻¹``, ``ω⁰`` by which [the standard tetrad](@ref "The standard
tetrad") was made regular.  Those are ``ω⁻⁽¹⁺ᵇ⁾`` with ``b`` the boost
weight (``+1``, ``0``, ``-1`` for ``l̃``, ``m̃``, ``ñ``).  Under BMS
each leg is regularized with ``ω' = κω`` in place of ``ω`` *and*
transformed by the Lorentz boost, and the two contributions combine to
the same ``1/κ`` for every leg:

| leg | regularize | boost | net |
|-----|:----------:|:-----:|:---:|
| ``l̃`` | ``κ⁻²`` | ``κ⁺¹`` | ``1/κ`` |
| ``m̃`` | ``κ⁻¹`` | ``κ⁰`` | ``1/κ`` |
| ``ñ`` | ``κ⁰`` | ``κ⁻¹`` | ``1/κ`` |

The differing regularization powers ``κ⁻⁽¹⁺ᵇ⁾`` are exactly cancelled
by the boost weights ``κᵇ``, leaving the uniform Weyl factor ``1/κ``.

Another important set of relations is the action of the vectors on
their coordinates, such as

```math
\begin{aligned}
ñ(u) &= \sqrt{2} &\qquad\qquad
ñ'(u') &= \sqrt{2} \\
m̃(u) &= 0 &\qquad\qquad
m̃'(u') &= 0,
\end{aligned}
```

and so on.  These relations allow us to easily determine how the
tetrads transform under BMS transformations.

Recall that the tangent space of ``ℐ⁺`` is spanned by ``m̃``, ``m̄̃``,
and ``ñ``, while ``l̃`` is transverse to it.  Importantly, *the same
statement applies with the primed tetrad*.  In particular, ``ñ'``
must be expressible as a linear combination of ``m̃``, ``m̄̃``, and
``ñ``.  We can now begin to use the specific form of the BMS
transformations to determine the coefficients in this expansion:

```math
u' = κ(θ, ϕ) \left[u - εᵅα(θ, ϕ)\right], \qquad
θ' = θ'(θ, ϕ), \qquad
ϕ' = ϕ'(θ, ϕ).
```

The point of these last two expressions is to emphasize that the
angular coordinates are independent of the time.  Now, we can expand
in the unprimed coordinates:

```math
ñ' = a ñ + b m̃ + c m̄̃,
```

but we know that since ``ñ'`` is real, we must have ``c = b̄``.  We
can then determine the coefficients by applying both sides to the
coordinates.  For example, applying both sides to ``u``, we have

```math
\begin{aligned}
ñ'(u) &= a ñ(u) + b m̃(u) + b̄ m̄̃(u) \\
&= a \sqrt{2}.
\end{aligned}
```

But since ``u = u' / κ + εᵅα``, we know that ``ñ'(u) = \sqrt{2} /
κ``, so we have ``a = 1/κ``.  We can similarly apply both sides to
``θ`` and ``ϕ`` to determine that ``b = b̄ = 0``, so that

```math
ñ' = \frac{1}{κ} ñ.
```

Next, we expand

```math
m̃' = d m̃ + e m̄̃ + f ñ.
```

The orthogonality relation ``m̃' ⋅ m̃' = 0`` implies that ``de = 0`` — either ``d`` or ``e`` must be zero.  Continuity and the condition that ``m̃' = m̃`` for the trivial transformation implies that we must have ``e = 0``, so

```math
m̃' = d m̃ + f ñ.
```

Applying both sides to ``u'``, we get

```math
m̃'(u') = 0 = d m̃(u') + f \sqrt{2} κ,
```

which we solve to find

```math
\frac{f}{d} = -\frac{m̃(u')}{\sqrt{2} κ} = \frac{ðu'}{2 κ}.
```

We can write

```math
m̃' = d \left(m̃ + \frac{ðu'}{2 κ} ñ\right).
```

It remains to fix ``d``.  The Lorentz part of the BMS transformation
is a single rotor ``L ∈ \mathrm{Spin}⁺(3,1)``.  Carried to the pole
and factored [as for the Lorentz group](@ref "Iwasawa and Hopf"), its
three pieces act on ``𝐦`` one at a time.  The null rotation ``N``
fixes ``ñ`` and reproduces the shift ``m̃ ↦ m̃ + (ðu'/2κ)\,ñ`` found
above, and is just the [null rotation worked out earlier](@ref "Null
rotations") with the roles of ``\boldsymbol{ℓ}`` and ``𝐧`` exchanged.
Its coefficient is linear in ``u`` because ``m̃`` and ``m̃'`` are
tangent to *different* cuts — of constant ``u`` and of constant ``u'``
— and the angle-dependent rescaling ``u' = κu`` tilts those cuts apart
at a rate set by the lever arm ``u``.  The shift therefore vanishes on
the boost's fixed cut ``u = 0`` and grows linearly up the generator,
as the null rotation's shear should.  That leaves the boost ``A`` and
the spin ``K`` to supply ``d``, and each is a one-line computation.

The ``K`` factor is the Hopf-fiber rotor ``R_γ =
\exp\left[\tfrac{γ}{2} 𝐱𝐲\right]``, a rotation about the null
direction.  On the pole dyad ``𝐦 = (𝐱 + 𝐈₃𝐲)/\sqrt{2}`` it acts as
a pure phase,

```math
R_γ\, 𝐦\, R̃_γ = e^{𝐱𝐲γ}\, 𝐦,
```

because ``𝐈₃𝐦 = 𝐱𝐲𝐦``: on the dyad the spatial pseudoscalar
``𝐈₃`` and the screen pseudoscalar ``𝐱𝐲`` coincide, so ``𝐱𝐲`` is
the dyad's own ``i``.  A quantity of spin weight ``s`` thus acquires a
phase ``e^{s𝐱𝐲γ}``.

The ``A`` factor is the boost ``R_{φₐ} = \exp\left[\tfrac{φₐ}{2}
𝐭𝐳\right]`` along the radial null direction.  It rescales the generators,

```math
R_{φₐ}\, \boldsymbol{ℓ}\, R̃_{φₐ} = e^{φₐ}\, \boldsymbol{ℓ},
\qquad
R_{φₐ}\, 𝐧\, R̃_{φₐ} = e^{-φₐ}\, 𝐧.
```

This is the very rescaling behind ``ñ' = ñ/κ``, which identifies the
conformal factor as

```math
κ = e^{φₐ}.
```

For the generator this boost is the whole story — ``ñ`` involves no
power of ``ω`` (it was regularized with ``ω⁰``), so ``ñ' = ñ/κ`` is
the boost acting alone.  The dyad is different: the boost leaves the
spacetime multivector ``𝐦`` fixed, up to the reinterpretation of
``𝐈₃``, so its ``1/κ`` comes instead from regularizing with ``ω' =
κω``.  Both routes are just the one uniform Weyl factor above, so
``|d| = 1/κ``.  That conformal origin is why the factor is invisible
to the normalization ``m̃' ⋅ m̄̃' = 1``, which holds in the primed
frame and would wrongly give ``|d| = 1`` in the unprimed metric.

Collecting the boost and the spin,

```math
d = \frac{e^{𝐱𝐲γ}}{κ} = e^{-φₐ + 𝐱𝐲γ},
```

whose exponent ``-φₐ + 𝐱𝐲γ`` is precisely the spin-boost weight of the dyad.
Three checks: a pure supertranslation has ``L = 1``, so ``φₐ = γ = 0`` and
``d = 1``, leaving only the null-rotation shift; a spatial rotation has ``φₐ =
0``, so ``κ = 1`` and ``d = e^{𝐱𝐲γ}`` is a pure phase — the dyad simply
spins; and a boost along the line of sight has ``γ = 0``, giving the real ``d
= 1/κ`` of aberration.  The dyad therefore transforms as

```math
m̃' = \frac{e^{iγ}}{κ}\left(m̃ + \frac{ðu'}{2 κ}\, ñ\right).
```

With ``ñ'`` and ``m̃'`` in hand, the transverse leg ``l̃'`` follows from
orthogonality alone — no further input from the transformation is needed.
Being transverse, it requires the full basis,

```math
l̃' = A\, l̃ + B\, ñ + C\, m̃ + D\, m̄̃,
```

and the four coefficients are fixed by demanding that ``l̃'`` satisfy
the standard inner products with the primed tetrad.  The orthogonality
conditions are conformally invariant, but the normalizations hold in
``ĝ'`` — so in the unprimed metric they read ``-1/κ²`` and ``1/κ²``,
since ``ĝ' = κ²ĝ``.  Using ``l̃ ⋅ ñ = -1``, ``m̃ ⋅ m̄̃ = 1``, and the
vanishing of every other unprimed product:

- ``l̃' ⋅ ñ' = -1/κ²`` in ``ĝ``: with ``ñ' = ñ/κ``, only the ``l̃`` term
  survives, so ``-A/κ = -1/κ²`` and ``A = 1/κ``;
- ``l̃' ⋅ m̃' = 0`` (frame-independent): with ``m̃' = d\,(m̃ +
  \tfrac{ðu'}{2κ}ñ)``, this gives ``D = A\,\tfrac{ðu'}{2κ} =
  \tfrac{ðu'}{2κ²}``;
- ``l̃' ⋅ m̄̃' = 0``: the conjugate condition gives ``C =
  \tfrac{\overline{ðu'}}{2κ²}``;
- ``l̃' ⋅ l̃' = 0``: the surviving terms ``-2AB + 2CD`` give ``B = CD/A =
  \tfrac{|ðu'|²}{4κ³}``.

Collecting these,

```math
l̃' = \frac{1}{κ}\, l̃ + \frac{\overline{ðu'}}{2κ²}\, m̃ + \frac{ðu'}{2κ²}\, m̄̃
      + \frac{|ðu'|²}{4κ³}\, ñ
    = \frac{1}{κ}\left( l̃ + \frac{\overline{ðu'}}{2κ}\, m̃ + \frac{ðu'}{2κ}\, m̄̃
      + \frac{|ðu'|²}{4κ²}\, ñ \right).
```

This is exactly the ``\boldsymbol{ℓ} ↔ 𝐧`` image of the null rotation
``𝐑\,𝐧\,𝐑̄ = 𝐧 + \boldsymbol{ξ} + \tfrac12 ξ²\boldsymbol{ℓ}``: writing the
screen vector ``\boldsymbol{ξ} = \tfrac{\overline{ðu'}}{2κ}m̃ +
\tfrac{ðu'}{2κ}m̄̃`` (so that ``\tfrac12 ξ² = |ðu'|²/4κ²``), the transverse
leg picks up a linear screen term and a quadratic ``ñ`` term — and, like
every other leg, the uniform Weyl factor ``1/κ`` out front.

## Past null infinity

The story at ``ℐ⁻`` is the mirror image, and almost everything carries
over verbatim under the exchange ``l̃ ↔ ñ`` together with ``u → v``
(advanced time).  The [regular tetrad](@ref "The standard tetrad")
here is

```math
l̃ = \sqrt{2}\, ∂_v, \qquad
m̃ = \frac{1}{\sqrt{2}}\left(∂_θ + \frac{i}{\sin θ} ∂_ϕ\right), \qquad
ñ = \frac{1}{\sqrt{2}}\, ∂_ω,
```

so now ``l̃`` is the generator tangent to ``ℐ⁻`` and ``ñ`` is
transverse — the opposite of ``ℐ⁺``.  The BMS action keeps the same
shape,

```math
v' = κ(θ, ϕ)\left[v - εᵅ α(θ, ϕ)\right],
\qquad
κ = \frac{1}{γ(1 + v⃗ ⋅ n̂)},
```

with two points of convention worth flagging:

- **Antipodal labeling.**  Following [PenroseRindler_1984](@citet),
  ``ℐ⁻`` is labeled by the observer's past light cone — the directions
  radiation *arrives from* — so its section is ``𝐧 = (1, -n̂)``,
  antipodal to ``ℐ⁺``.  This is the single sign ``εᴵ = -1``, and it is
  why ``κ`` includes ``1 + v⃗⋅n̂`` rather than ``1 - v⃗⋅n̂`` (see
  [Future and past null infinity](@ref scri_pm_conventions)).
- **Same uniform weight.**  The conformal frames still satisfy ``ĝ' =
  κ²ĝ``, so every leg again picks up ``1/κ``; nothing in the weight
  argument changes.

Repeating the ``ℐ⁺`` derivation under ``l̃ ↔ ñ`` then gives

```math
\begin{aligned}
l̃' &= \frac{1}{κ}\, l̃, \\
m̃' &= \frac{e^{iγ}}{κ}\left(m̃ + \frac{ðv'}{2κ}\, l̃\right), \\
ñ' &= \frac{1}{κ}\left( ñ + \frac{\overline{ðv'}}{2κ}\, m̃
        + \frac{ðv'}{2κ}\, m̄̃ + \frac{|ðv'|²}{4κ²}\, l̃ \right).
\end{aligned}
```

The only structural change is that the null rotation now fixes ``l̃``
(the ``ℐ⁻`` generator) and shears ``m̃`` and ``ñ`` *along* it.  That
same exchange reverses the peeling tower of the [field
components](@ref "BMS action on fields"): the component left unmixed by
the null rotation is the one assembled from the most factors of the
generator — ``ψ₄`` at ``ℐ⁺`` (fixed by the rotation about ``ñ``), but
``ψ₀`` at ``ℐ⁻`` (fixed by the rotation about ``l̃``).
