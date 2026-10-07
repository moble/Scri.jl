# [BMS Action on the Tetrad — a Geometric Route](@id tetrad_geometric)

The derivation on [the previous page](@ref "BMS Action on the Tetrad")
obtains the transformation of the tetrad by applying each leg to the
coordinates and solving for the coefficients.  That route is short and
easy to verify, but it has two shortcomings.  The conformal factor and
the spin phase are obtained by separate arguments, one from the
conformal weight of the legs and one from the action of a rotation on
the dyad, and their identification with the ``A`` and ``K`` factors of
[the ``KAN`` decomposition](@ref "Iwasawa's ``KAN`` decomposition") is
stated rather than derived.  The null-rotation parameter, meanwhile,
appears as a coordinate derivative, with no evident relation to the
group structure at all.

This page derives the same law along a different route, in which each
factor of the answer is produced by a structural feature of the group.
The ``KAN`` decomposition supplies the spin and the boost, and the
semidirect-product structure of BMS supplies the null rotation.  The
two derivations agree, and the coordinate computation of the previous
page remains the most direct verification; what is gained here is an
account of *why* the answer takes the form it does.

## A Lorentz transformation distinguishes one cut

It is worth recalling why the ``KAN`` decomposition, which is such a
natural description of the Lorentz group, does not immediately describe
the BMS action on the tetrad.

A Lorentz transformation is a generalized rotation about the spacetime
origin.  A spatial rotation is taken about a spatial origin, and a boost
is likewise taken about the time origin, so the transformation as a
whole is a rotation about the full spacetime origin.  Consequently the
objects it describes most naturally are those attached to that origin,
and in particular the null cone through it.  The ``KAN`` decomposition
inherits this: it is stated relative to a null direction *at a point*.

The same fact appears at ``ℐ`` in the time law.  A Lorentz
transformation acts as ``t' = κ(k̂)\,t``, so the cuts it preserves
satisfy ``t\left(κ - 1\right) = 0``.  Unless ``κ ≡ 1`` — that is, unless
the transformation is a rotation — the only such cut is ``t = 0``.  A
boost therefore distinguishes exactly one cut of ``ℐ``, and it is the
cut on which the two observers' origins coincide.

This suggests splitting the problem in two.  First, describe the action
of the transformation *on its own distinguished cut*, where we may hope
that ``KAN`` applies directly.  Second, account for the displacement
from that cut to whichever cut the data actually occupies.  The
remainder of the page follows that plan.

## Cuts as the homogeneous space

The second half of the plan requires knowing how BMS relates the
distinguished cuts of different elements, and the semidirect-product
structure answers this.

Recall from [the BMS page](@ref bms_group) that ``\text{BMS} = ℒ ⋉_ρ 𝒮``
with ``ρ(Λ)α = α/κ``, that ``𝒮`` is a normal subgroup, and that ``ℒ``
is not.  The failure of ``ℒ`` to be normal is precisely what we need:
conjugating a Lorentz transformation by a supertranslation produces an
element that is not a Lorentz transformation.  Explicitly, using the
group law of that page,

```math
\left(1, β\right)\left(Λ, 0\right)\left(1, β\right)^{-1}
= \left(Λ,\; \frac{β ∘ Λ}{κ} - β\right).
```

The supertranslation content on the right vanishes only when ``β ∘ Λ =
κβ``, which for a boost fails for every nonzero ``β``.  So there is a
whole family of Lorentz subgroups of BMS, conjugate to one another by
supertranslations, and each distinguishes its own cut.  The situation is
identical in structure to the Poincaré group, where the Lorentz
subgroups are conjugate by translations and each distinguishes its own
origin event.

The analogy extends to the homogeneous spaces.  Minkowski spacetime is
the coset space of the Poincaré group by a Lorentz subgroup, its points
labeling the possible choices of origin.  In the same way,

```math
\left\{\text{cuts of } ℐ \right\} \;=\; \text{BMS} / ℒ ,
```

the cuts labeling the possible choices of *origin cut*.  A choice of
origin in special relativity is replaced, in the BMS setting, by a
choice of cut.

!!! note "Why there is no ``KAN`` decomposition of BMS"

    Iwasawa's theorem applies to semisimple groups, and BMS is not
    semisimple: the supertranslations are an abelian normal subgroup.
    There is therefore no useful ``KAN``-like factorization of BMS
    itself.  But this is again the Poincaré situation, and the same
    remedy applies.  One does not seek an Iwasawa decomposition of the
    Poincaré group; one uses the semidirect structure to reduce to the
    stabilizer of a point, and that stabilizer — the Lorentz group — is
    semisimple.  Here the stabilizer of a cut is likewise a Lorentz
    group, so the correct generalization is not "``KAN`` for BMS" but
    "``KAN`` relative to a chosen cut", together with the conjugation
    above to move between cuts.

## The celestial sphere and its conformal frames

Before using ``KAN`` we need to know what it measures.  The answer is
that it measures the action on a particular four-dimensional space
attached to the celestial sphere, and that both the conformal factor and
the spin phase are components of a single object there.

Recall from [the conformal-factor
section](@ref "Lorentz transformations ``ℒ`` and the conformal factor
``κ``") that a Lorentz transformation acts on the celestial sphere by a
conformal map, with ``{dΩ'}² = κ²\, dΩ²``.  Being conformal, the
differential of that map at any point is a *similarity* of the tangent
plane: a positive rescaling composed with a rotation.  Identifying the
tangent plane with ``ℂ`` by way of the dyad, a similarity is
multiplication by a single complex number.  Its modulus is ``κ``, and its
argument is an angle that we will call ``-γ``, the sign chosen for later
convenience.  A conformal map of the sphere therefore supplies exactly
two real numbers at each point, and there is no room for more.

The group-theoretic counterpart is the following.  Consider the set of
pairs

```math
\left(𝐤,\; \left[𝐦\right]\right),
```

where ``𝐤`` is a future-pointing null *vector* — not merely a ray — and
``𝐦`` is a complex null vector with ``𝐦 ⋅ 𝐦̄ = 1`` and ``𝐦 ⋅ 𝐤 = 0``,
taken modulo the shift ``𝐦 ↦ 𝐦 + c\,𝐤``.  The null vectors form a
three-dimensional set, and for fixed ``𝐤`` the admissible ``[𝐦]``
differ only by a phase, so the pairs form a four-dimensional space.  We
will call it the **bundle of conformal frames**, since a pair records a
point of the sphere, a scale, and an orientation in the tangent plane.

``\mathrm{Spin}^+(3,1)`` acts transitively on this bundle, and the
stabilizer of a pair is easy to identify.  An element fixing ``𝐤``
cannot contain a boost along it, since a boost rescales ``𝐤``; and an
element fixing ``[𝐦]`` cannot contain a rotation about ``𝐤``, since
such a rotation changes the phase.  What remains is exactly the
nilpotent subgroup, and therefore

```math
\left\{\text{conformal frames}\right\}
\;=\; \mathrm{Spin}^+(3,1) / N_{\boldsymbol{ℓ}} ,
```

whose dimension is ``6 - 2 = 4`` as required.  Reading the ``KAN``
factorization of a rotor against this identification gives the
statement we were after:

- the ``𝕊²`` part of ``K`` determines the point of the sphere;
- the ``M`` part of ``K``, the Hopf fiber, determines the orientation in
  the tangent plane;
- ``A`` determines the scale, since it is the factor that rescales
  ``𝐤``; and
- ``N`` acts trivially, since it is the stabilizer.

So the modulus and argument of the conformal map's differential are the
``A`` and ``M`` factors of ``KAN``, and they are the only quantities
that a Lorentz transformation can contribute to any object built from a
point of the sphere together with a dyad.  In particular, writing the
boost as ``𝐑_{φₐ} = \exp\left[\tfrac{φₐ}{2}𝐭𝐳\right]`` and the
Hopf-fiber rotation as ``𝐑_γ = \exp\left[\tfrac{γ}{2}𝐱𝐲\right]``, we
have

```math
κ = e^{φₐ},
\qquad
𝐑_γ\, 𝐦\, 𝐑̃_γ = e^{iγ}\, 𝐦 .
```

The second of these deserves a word.  On the dyad the spatial
pseudoscalar and the screen pseudoscalar agree, ``𝐈₃𝐦 = 𝐱𝐲\,𝐦``, so
``𝐱𝐲`` is the dyad's own unit imaginary and the rotation acts as a pure
phase; a quantity of spin weight ``s`` accordingly acquires ``e^{isγ}``.
Note that this generator is the negative of the one written for ``M`` on
[the Lorentz page](@ref IwasawaHopfLevi), a sign chosen here so that the
dyad phase comes out positive.

!!! details "Why the null rotation cannot contribute a weight"

    That ``N`` is the stabilizer of a conformal frame gives an
    immediate reason for a fact established by other means on [the
    Lorentz page](@ref which_null_direction): no quantity constructed
    from a point of the sphere and a dyad can transform under ``N`` by
    any factor at all, because ``N`` moves neither.  A null rotation can
    act only on data that involve more than the conformal frame — that
    is, on data that involve the transverse leg, and hence the choice of
    cut.  This is why its effect must appear as a mixing of components
    rather than as a weight, and it is also why the mixing is bound up
    with the slicing.

## The law on the distinguished cut

We can now dispose of the first half of the plan.

On the cut distinguished by the transformation, the coordinate
derivative that produced the mixing on the previous page vanishes: with
``t' = κ\,t`` and ``κ`` independent of time, ``ð t' = t\, ðκ``, which is
zero at ``t = 0``.  The transition between the two tetrads on that cut
therefore contains no null rotation, and by the previous section it is
precisely the ``M`` and ``A`` factors of ``KAN``.

It remains to convert this into a statement about the legs that are
regular at ``ℐ``.  Each unregularized leg has a boost weight ``b``,
equal to ``+1``, ``0`` and ``-1`` for ``\boldsymbol{ℓ}``, ``𝐦`` and
``𝐧``, and the ``A`` factor multiplies it by ``κ^{b}``; the ``M``
factor contributes the phase ``e^{iγ}`` on the dyad alone.  [The regular
tetrad](@ref "Tetrad") is obtained by multiplying the leg of weight
``b`` by ``ω^{-(1+b)}``, and the primed observer regularizes with ``ω' =
κω``.  The ratio of the primed regular leg to the unprimed one is
therefore

```math
\left(\frac{ω'}{ω}\right)^{-(1+b)} κ^{b}
= κ^{-(1+b)}\, κ^{b}
= \frac{1}{κ}
\qquad\text{for every } b .
```

The graded regularization exactly cancels the boost grading, and every
leg acquires the same factor ``1/κ``.  This is the uniform Weyl factor
of the previous page, obtained there from the relation ``ĝ' = κ²ĝ``
between the two conformal metrics; the two accounts are equivalent, but
the present one exhibits the cancellation rather than the result.  The
spin phase involves no power of ``ω`` and so survives untouched.  Hence,
on the distinguished cut,

```math
ñ' = \frac{1}{κ}\, ñ,
\qquad
m̃' = \frac{e^{iγ}}{κ}\, m̃,
\qquad
l̃' = \frac{1}{κ}\, l̃ .
```

There is a second route to the dyad relation which is worth recording,
because it identifies ``γ`` independently of the group theory.  The
legs ``m̃`` and ``m̃'`` are *coordinate* dyads, built from the
respective angular coordinates, so they are related by the inverse
Jacobian of the angular map.  That map is conformal, so its Jacobian is
the single complex number of the previous section, of modulus ``κ`` and
argument ``-γ``; its inverse is ``e^{iγ}/κ``.  The two derivations
agree, which fixes the sign convention for ``γ`` and confirms that the
angle appearing in the tetrad law is the Hopf-fiber angle of ``K``.

## What a difference of cuts does to the tetrad

Away from the distinguished cut the two observers no longer agree on
which two-surface through a given point is the cut, and by the previous
section that disagreement is the only difference between their tetrads
that remains to be accounted for.  Two questions follow.  What kind of
transformation of the tetrad does a difference of cuts induce, and by
what quantity is the difference measured?  The first is settled by the
geometry of the surface, the second by an elementary slope.

Both tetrads are adapted to ``ℐ``, meaning that each has its generator
leg tangent to the generators of the surface, and the screen of each is
the tangent plane to that observer's own cut.  Two cross-sections of a
single null hypersurface can differ only by a displacement along its
generators, since that is the one direction the surface itself
distinguishes.  Their tangent planes therefore differ by a **shear
along the generator**, by which we mean the following linear map of the
three-dimensional tangent space of ``ℐ``.  It fixes the generator, and
it sends each screen vector to itself plus a multiple of the generator,
that multiple being a linear functional of the screen vector:

```math
ñ \;\longmapsto\; ñ,
\qquad
𝐯 \;\longmapsto\; 𝐯 + λ(𝐯)\, ñ
\quad\text{for } 𝐯 \text{ in the screen} .
```

The map is thus the identity modulo the generator, and the nilpotent
piece added to it takes the screen into the generator line and
annihilates the generator.  This is a shear in the ordinary sense —
each screen vector slides along the generator by an amount proportional
to itself, as the cards of a deck slide by an amount proportional to
their height — rather than a uniform displacement, which would add the
same multiple of the generator to every vector and would not be linear.
The functional ``λ`` has two real components, or one complex component
on the dyad, so the shear has two parameters.

A shear of the screen along a null direction is, by definition, a null
rotation about that direction.  At ``ℐ⁺`` the generator is ``ñ``, so the
transition lies in the four-parameter stabilizer ``M A N_{𝐧}``, and
applying the map to the dyad gives

```math
m̃ \;\longmapsto\; m̃ + b\, ñ,
\qquad
b = λ\!\left(m̃\right),
```

for some complex ``b``, which it remains to determine.

!!! warning "Two meanings of "shear" "

    The word is used here in the linear-algebra sense just defined, and
    not in the Newman–Penrose sense, where the shear ``σ`` of a null
    congruence is the trace-free part of its transverse deformation.
    The two are related — a reslicing does change the Newman–Penrose
    shear of the cut, and that change is [the inhomogeneous term in the
    transformation of ``σ``](@ref "Strain, shear, and news") — but they
    are different objects, and only the linear-algebra meaning is
    intended in this section.

The parameter is fixed by the one property that singles the primed cut
out, namely that it is a level set of ``t'``.  A vector lies in the
primed screen if and only if it annihilates ``t'``, so

```math
0 = \left(m̃ + b\, ñ\right)\!\left(t'\right)
  = m̃\!\left(t'\right) + b\, ñ\!\left(t'\right),
\qquad\text{hence}\qquad
b = -\frac{m̃\!\left(t'\right)}{ñ\!\left(t'\right)} .
```

This is nothing more than the formula for the slope of one graph over
another.  Both cuts are transverse to the generator, so within the
three-dimensional tangent space of ``ℐ`` each of their tangent planes is
a graph over the other, with the generator as the fiber direction.  The
numerator is then the rate at which ``t'`` changes along the unprimed
screen, the denominator the rate at which it changes along the
generator, and the ratio is the distance one must slide along the
generator, per unit angular displacement, in order to remain on the
primed cut.  **The tilt of one cut relative to the other is therefore
exactly what the dyad records**, and this is the link between the
coordinate transformation and the tetrad.

Both derivatives are immediate.  The generator leg is normalized so that
``ñ(t) = \sqrt2``, and it annihilates the angles, while ``t'`` differs
from ``κ t`` only by a function of the angles; hence

```math
ñ\!\left(t'\right) = \sqrt2\, κ .
```

The dyad, being tangent to the unprimed cut, annihilates ``t``, so only
the angular dependence of ``t'`` contributes.  For any function of the
angles alone the regular dyad gives ``m̃(f) = -ð f/\sqrt2``, so

```math
m̃\!\left(t'\right) = -\frac{ð t'}{\sqrt2} .
```

Together these give

```math
b = \frac{ð t'}{2 κ} ,
```

the factor of two arising from the two leg normalizations, one from the
dyad and one from the generator.

Two remarks before continuing.  First, no term along the conjugate dyad
appears, and the reason is worth isolating.  The induced metric on ``ℐ``
is degenerate along ``ñ``, so a *cut metric* is well defined on the
quotient of the tangent space by the generator; the shear is the identity
on that quotient and so distorts nothing there, while the residual map
of the quotient onto itself is an orientation-preserving similarity,
since a proper orthochronous transformation acts on the sphere
conformally and preserves orientation.  An orientation-preserving
similarity of a two-plane is complex-linear, so it multiplies the dyad by
a single complex number and cannot produce a conjugate term.  Second,
the tetrad is now seen to depend on the two coordinate systems through
the single function ``ð t'`` and nothing else.  What remains is to
determine that function, and the answer is structural.

The transverse leg requires no further input, but it does require
stepping off the surface.  The shear above was a map of the tangent
space of ``ℐ``, which does not contain ``l̃``; the null rotation that
extends it to the full tangent space of the spacetime is the one
[worked out earlier](@ref "Null rotations"), with the two null vectors
exchanged, and it acts on the remaining leg as

```math
\boldsymbol{ℓ} \;\longmapsto\;
\boldsymbol{ℓ} + \boldsymbol{ξ} + \tfrac12 ξ²\, 𝐧 .
```

Here the displacement is no longer purely linear: the transverse leg
acquires a screen term at first order in the parameter and a generator
term at second order.  Expanding the screen vector in the dyad as
``\boldsymbol{ξ} = \bar b\, 𝐦 + b\, 𝐦̄``, so that ``𝐦 ⋅
\boldsymbol{ξ} = b`` and ``\tfrac12 ξ² = |b|²``, those two terms are the
linear and quadratic coefficients of the law collected below.  Imposing
orthogonality on the primed tetrad, as on the previous page, gives the
same result.

## The tilt of the cuts is a commutator

It remains to evaluate ``ð t'``, and the point of doing so within the
group is that the result then requires no computation in coordinates.

Suppose the data occupy a cut other than the one the transformation
distinguishes.  To describe the action there we conjugate by the
translation that moves the origin cut, and, as recorded above, that
conjugation does not return a Lorentz transformation.  Taking ``β = c``
constant — an ordinary time translation — the conjugate is

```math
\left(1, c\right)\left(Λ, 0\right)\left(1, c\right)^{-1}
= \left(Λ,\; c\left(\frac{1}{κ} - 1\right)\right),
```

whose supertranslation content is nonzero, and angle-dependent through
``κ``, precisely because translations and boosts do not commute.
Substituting into the time law ``t' = κ\left(t - c_α α\right)`` gives

```math
t' = κ\left(t + c_α c\right) - c_α c ,
```

which confirms that the distinguished cut has moved to ``t = -c_α c``,
the image of ``t = 0`` under the translation, as it must.  The angular
derivative is then

```math
ð t' = \left(t + c_α c\right) ðκ .
```

The null-rotation parameter is therefore proportional to the *offset
between the cut carrying the data and the cut the transformation
distinguishes*.  It vanishes on the distinguished cut and grows linearly
away from it, and the growth rate is fixed by ``ðκ`` alone.  Nothing
here was computed from the coordinates; the linear growth is the twist
``ρ(Λ)α = α/κ`` of the semidirect product, evaluated on the offset.

Repeating the conjugation with a general ``β`` produces the second term
of the general formula, and it is more compact simply to quote the
result.  For an arbitrary element ``(Λ, α)``, with ``t' = κ(t - c_α
α)``,

```math
ð t' = \left(ðκ\right)\left[t - c_α α\right] - κ\, c_α\, ðα ,
```

in which the first term is the offset just discussed and the second is
the tilt a supertranslation imprints on the cuts directly.  The second
term survives when ``Λ = 1``, where there is no Lorentz transformation
at all and hence no ``KAN`` factorization to appeal to — which is the
structural reason the mixing parameter cannot be a factor of a rotor.

!!! note "An analogy worth keeping"

    The dependence of ``ð t'`` on the choice of origin cut is the same
    in form as the dependence of angular momentum on the choice of
    origin point: the value shifts by the offset multiplied by a
    fixed quantity, ``ðκ/2κ`` here in place of the momentum.  This is
    not merely an analogy of form.  The supertranslation ambiguity of
    BMS angular momentum is the same phenomenon, arising from the same
    failure of ``ℒ`` to be normal.

## The assembled law

Collecting the uniform Weyl factor, the spin phase, and the shear:

```math
\begin{aligned}
ñ' &= \frac{1}{κ}\, ñ, \\
m̃' &= \frac{e^{iγ}}{κ}\left(m̃ + b\, ñ\right), \\
l̃' &= \frac{1}{κ}\left(l̃ + \bar b\, m̃ + b\, m̄̃ + |b|²\, ñ\right),
\end{aligned}
\qquad
b = \frac{ð t'}{2κ},
```

which is the law of [the previous page](@ref "BMS Action on the
Tetrad").  Each ingredient now has a stated origin: the factor ``1/κ``
from the cancellation of the boost grading against the graded
regularization, with ``κ = e^{φₐ}`` the ``A`` factor; the phase from the
Hopf fiber of ``K``, equivalently from the rotation part of the
conformal map's differential; and the shear from the offset between the
data cut and the cut the transformation distinguishes, by way of the
semidirect twist.

Three limiting cases confirm the accounting.  A pure supertranslation
has ``Λ = 1``, hence ``κ = 1`` and ``γ = 0``, and only the shear
survives, with ``b = -c_α\, ðα/2``.  A spatial rotation has ``κ ≡ 1``
and ``ðκ = 0``, so ``b = 0`` and the dyad merely acquires its phase;
consistently, a rotation distinguishes every cut rather than one.  A
boost along the line of sight has ``γ = 0`` and, on its distinguished
cut, ``b = 0``, leaving the real factor ``1/κ`` of aberration.

## Past null infinity

Everything above carries over to ``ℐ⁻`` under the exchange of the two
null legs together with the replacement of retarded by advanced time.
The generator there is ``l̃`` rather than ``ñ``, so the shear is a null
rotation about ``l̃`` and the transition lies in ``M A
N_{\boldsymbol{ℓ}}``; the conformal factor involves ``1 + v⃗ ⋅ k̂``
rather than ``1 - v⃗ ⋅ k̂``, following the [antipodal
labeling](@ref scri_pm_conventions) of that surface.  The argument for
the shear is unchanged, since it used only that both tetrads are adapted
to one null hypersurface, and the identification of the conformal factor
and spin phase with the ``A`` and ``M`` factors is unchanged, since it
used only the action on the conformal frames of the sphere.  The
resulting law is the one collected at the end of [the previous
page](@ref "BMS Action on the Tetrad").

## Relation to the coordinate derivation

The two derivations are complementary rather than competing.  The
coordinate derivation is the more direct and the easier to check: it
requires only the action of the legs on the coordinates, and it produces
the null-rotation parameter immediately, in the form in which the code
computes it.  The route taken here is longer, but it answers questions
the coordinate route leaves open — why exactly two real numbers appear
in the spin–boost part, why they are the ``A`` and ``M`` factors of
``KAN``, why the mixing parameter is linear in the retarded time, and
why the mixing cannot be read off from any factorization of the Lorentz
rotor.  For the [convention dependence](@ref
convention_dependence_tetrad_future) of the result, and for the
consequences for the [field components](@ref "BMS Action on Fields"),
the previous page and its successor remain the references.
