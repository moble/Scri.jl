# BMS and the tetrad

The next question is how [the standard tetrad](@ref "The standard
tetrad") transforms under coordinate transformations — and
specifically BMS transformations.  One *very* important point is that
the tetrad is *defined* in terms of the coordinates.  In particular,
if we have an unprimed coordinate system and a primed coordinate
system, we actually have *two distinct* tetrads.  For example, we have
these tetrads that are regular near ``ℐ⁺``:

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

With the conformal metric, the inner products of the tetrads retain
their usual values: ``l̃ ⋅ ñ = -1`` and ``m̃ ⋅ m̄̃ = 1``, with all
other combinations zero in the unprimed frame, and equivalent
statements in the primed frame.  Another important set of relations is
the action of the vectors on their coordinates, such as

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

The dyad is a unit tangent vector of the celestial sphere, so under the same
dilation it carries the opposite weight and scales by ``1/κ`` — the same
factor, for the same reason, as the generator.  (The boost leaves the
spacetime multivector ``𝐦`` itself fixed, up to the reinterpretation of
``𝐈₃``; the ``1/κ`` is the conformal weight of the regularized leg, which is
exactly why it stays hidden from the naive normalization ``m̃' ⋅ m̄̃' = 1`` —
that product is taken in the *primed* conformal metric.)

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
