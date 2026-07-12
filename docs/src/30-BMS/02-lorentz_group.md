# The Lorentz Group

Here, our focus is on using Geometric Algebra (GA) to understand and
implement Lorentz transformations — and to decompose those
transformations into conventional pieces.  In the context of spacetime
and Lorentz transformations, GA allows us to work with null tetrads
and the Lorentz group in a more intuitive and geometrically meaningful
way.  This treatment ties in surprisingly neatly with the approach via
stereographic coordinates and conformal transformations
[PenroseRindler_1984, PenroseRindler_1986, ODonnell_2003](@cite),
while retaining clearer connections to the underlying geometry.

## The standard basis and the null tetrad

We start with a standard basis of spacetime vectors, which we denote
as ``(𝐭, 𝐱, 𝐲, 𝐳)``.  We use signature ``{-}{+}{+}{+}``, meaning
that

```math
𝐭² = -1, \qquad 𝐱² = 𝐲² = 𝐳² = +1,
```

and all other products are zero.  We then define the *null* tetrad

```math
\begin{aligned}
\boldsymbol{ℓ} &= \frac{𝐭+𝐳}{\sqrt{2}}, \\
𝐦 &= \frac{𝐱+𝐈₃𝐲}{\sqrt{2}}, \\
𝐦̄ &= \frac{𝐱-𝐈₃𝐲}{\sqrt{2}}, \\
𝐧 &= \frac{𝐭-𝐳}{\sqrt{2}},
\end{aligned}
```

where ``𝐈₃ = 𝐈𝐭 = 𝐱𝐲𝐳`` is the spatial pseudoscalar, and
replaces the ``i`` used in the Newman-Penrose formalism.  See
[below](@ref "Reinterpreting ``i``") for details about this
replacement, but the important features here are that ``𝐈₃²=-1``, and
it commutes with spatial vectors but anticommutes with ``𝐭``.  This
means that the null tetrad is indeed null,

```math
\boldsymbol{ℓ}² = 𝐦² = 𝐦̄² = 𝐧² = 0,
```

and the only nonzero inner products are

```math
\boldsymbol{ℓ} \cdot 𝐧 = -1, \qquad 𝐦 \cdot 𝐦̄ = 1.
```

This simple tetrad is aligned with the axes of the standard basis, but
we can apply any Lorentz transformation to it to get a more general
null tetrad.  In particular, we can apply an ordinary rotation to
"point" ``\boldsymbol{ℓ}`` in any direction we like.  This suggests
the standard factorization of the Lorentz group, beginning with
rotations tied to the spherical coordinates, followed by rotations
about the radial direction, followed by boosts in the radial
direction, and finally followed by null rotations about the rotated
``\boldsymbol{ℓ}``.  We will discuss these transformations in more
detail below, but the important point is that we decompose a general
Lorentz transformation *specifically with respect to* either
``\boldsymbol{ℓ}`` or ``𝐧``.

## The Lorentz group

Before we decompose the (proper orthochronous) Lorentz group, we need
to understand how the Lorentz group shows up in Geometric Algebra.  The
answer is already implicit in [the primer](@ref "A Primer on Geometric
Algebra"): every proper orthochronous Lorentz transformation is a rotor
``R ∈ \mathrm{Spin}⁺(3,1)`` — an even, unit-norm product of vectors —
acting on a vector by conjugation,

```math
𝐯' = R\, 𝐯\, R̃,
```

with composition of transformations given by multiplication of rotors.
Each rotor is the exponential of a bivector, and the *type* of
transformation is fixed by the sign of that bivector's square — exactly
the three cases met in the primer.

A **rotation** comes from a *spatial* bivector, which squares to
``-1``.  Writing ``𝐢 = 𝐳𝐲``, ``𝐣 = 𝐱𝐳``, ``𝐤 = 𝐲𝐱`` for the three
spatial bivectors, a rotation by angle ``θ`` in the plane of a unit
bivector ``𝐁`` is

```math
R = \exp\left[\frac{θ}{2} 𝐁\right] = \cos\frac{θ}{2} + 𝐁 \sin\frac{θ}{2}.
```

These involve no factor of ``𝐭`` and generate the maximal compact
subgroup ``\mathrm{Spin}(3)``.

A **boost** comes from a *timelike* bivector, which squares to ``+1``.
A boost along ``𝐳`` with rapidity ``φ`` is

```math
R = \exp\left[\frac{φ}{2} 𝐭𝐳\right] = \cosh\frac{φ}{2} + 𝐭𝐳 \sinh\frac{φ}{2},
```

the trigonometric functions turned hyperbolic by ``(𝐭𝐳)² = +1``.  The
speed is ``β = \tanh φ`` and the Lorentz factor is ``γ = \cosh φ``.

Rotations and boosts together already generate the whole proper
orthochronous group, but there is a third, degenerate case that will
organize everything below: the **null rotation**, whose generating
bivector squares to *zero*.  We take that up next, then use all three
to decompose a general transformation in two complementary ways — one
adapted to a time axis, one adapted to a null direction.

## Null rotations

A null rotation is a very particular type of Lorentz transformation
(sometimes called a [parabolic
transformation](https://en.wikipedia.org/wiki/Lorentz_group#Parabolic))
that leaves a chosen null vector invariant.  The basic idea is that
the null rotation simultaneously boosts in a direction orthogonal to
the null vector, and rotates in the plane defined by that direction
and the null vector — with the rotation being exactly what is needed
to leave the null direction unchanged.

For simplicity, let us choose the null vector

```math
\boldsymbol{ℓ} = \frac{𝐭+𝐳}{\sqrt{2}}
```

as the invariant null vector.  Its complementary null vector is

```math
𝐧 = \frac{𝐭-𝐳}{\sqrt{2}},
```

which will be important.  Now, for any (not necessarily unit) vector
``\boldsymbol{ξ}`` in the ``𝐱``-``𝐲`` plane, the bivector
``\boldsymbol{ℓ ξ} = -\boldsymbol{ξ ℓ}`` generates a null rotation.
More specifically, this bivector generates a boost in the
``\boldsymbol{ξ}`` direction, and *simultaneously* a rotation in the
``\boldsymbol{ξ}``-``𝐳`` plane.  Specifically, we expand the product

```math
\boldsymbol{ℓ ξ} = \frac{𝐭\boldsymbol{ξ} + 𝐳\boldsymbol{ξ}}{\sqrt{2}}.
```

The ``𝐭\boldsymbol{ξ}`` term generates a boost in the
``\boldsymbol{ξ}`` direction, and the ``𝐳\boldsymbol{ξ}`` term
generates a rotation in the ``\boldsymbol{ξ}``-``𝐳`` plane.

Define the spinor

```math
𝐑 = \exp\left[ \frac{1}{2} \boldsymbol{ℓ ξ} \right].
```

Because ``\boldsymbol{ℓ}² = 0``, the exponential series terminates
after the second term, and we have

```math
𝐑 = 1 + \frac{1}{2} \boldsymbol{ℓ ξ},
```

which makes calculations particularly simple.  The results on the
basis ``(𝐭, 𝐱, 𝐲, 𝐳)`` are not enlightening, but the results on
the ``(\boldsymbol{ℓ}, 𝐧, 𝐦, 𝐦̄)`` basis are very interesting:

```math
\begin{aligned}
𝐑 \boldsymbol{ℓ} \bar{𝐑} &= \boldsymbol{ℓ}, \\
𝐑 𝐦 \bar{𝐑} &= 𝐦 + 𝐱\boldsymbol{ξℓ}, \\
𝐑 𝐦̄ \bar{𝐑} &= 𝐦̄ - 𝐱\boldsymbol{ξℓ}, \\
𝐑 𝐧 \bar{𝐑} &= 𝐧 + \boldsymbol{ξ} + \frac{1}{2} ξ² \boldsymbol{ℓ}.
\end{aligned}
```

This is just the usual conformal transformation of the tetrad, though
exhibited in a simpler and more geometric form.

## Cartan (polar) decomposition

There are two natural ways to factor a general rotor ``Λ ∈
\mathrm{Spin}⁺(3,1)``, and we will use both.  The first is the one
closest to physical intuition: every proper orthochronous Lorentz
transformation is a boost followed by a rotation,

```math
Λ = B(v⃗)\, R,
```

with ``B(v⃗) = \exp\left[\tfrac{φ}{2} 𝐭v̂\right]`` a pure boost of
velocity ``v⃗ = (\tanh φ)\, v̂`` and ``R ∈ \mathrm{Spin}(3)`` a pure
rotation.  This is the *polar decomposition* of the rotor — the
group-theoretic *Cartan decomposition* ``G = \exp(𝔭)\, K``, with the
boosts the noncompact part ``𝔭`` and the rotations the compact part
``K`` — and it is unique.  It is the sense in which two boosts in
different directions compose to "a boost plus a Wigner rotation":
multiplying two pure boosts yields a rotor whose polar decomposition has
a nontrivial ``R`` (see the [aberration page](@ref "Wigner rotation:
composition of non-collinear boosts")).

The factorization is cheap, with no transcendental functions and no
square roots of the transformation.  `Quaternionic.jl` represents a
Lorentz rotor as a quaternion whose four components are themselves
complex; in that representation a pure boost has real scalar part and
imaginary vector part, while a pure rotation is entirely real.  The
rotation is therefore fixed by the real parts alone,

```math
R = \frac{\mathrm{Re}\, Λ}{\lVert \mathrm{Re}\, Λ \rVert},
\qquad
B = Λ\, R̃,
```

with ``\mathrm{Re}`` taken component-by-component and the norm the
ordinary Euclidean one (the normalization works because
``\mathrm{Re}\, Λ = \cosh\tfrac{φ}{2}\, R``).  This is implemented as
`Quaternionic.vR`, which returns the pair ``(v⃗, R)`` in exactly the
``Λ = B(v⃗)\, R`` convention used here and by the [aberration
tests](@ref "Aberration of Gravitational Waves").

The polar decomposition is adapted to a choice of *time* axis: it splits
a transformation into the boost an observer feels and the rotation they
see.  For radiation we instead want a decomposition adapted to a *null*
direction, so that the factors line up with the spin and boost weights
of the fields.  That is the Iwasawa ``KAN`` decomposition, to which we
now turn.

## Iwasawa's ``KAN`` decomposition

[Knapp_1996](@Citet) describes [Iwasawa's
decomposition](https://en.wikipedia.org/wiki/Iwasawa_decomposition) —
a factorization of the Lorentz group or its double cover
``\mathrm{Spin}^+(3,1)`` (or any semisimple Lie group) into three
subgroups:

```math
G = KAN.
```

This decomposition is most useful when we have a preferred time axis
and preferred spatial direction.  We use those to pick out boosts
along that spatial direction, and null rotations that fix the
corresponding null vector.  The decomposition results in the following
subgroups, though we assume the particular forms given:

- **``K``** is the maximal compact subgroup ("K" from the German
  "kompakt"), which in the case of the Lorentz group is the rotation
  group ``\mathrm{SO}(3)``, and in the case of ``\mathrm{Spin}(3,1)``
  is ``\mathrm{Spin}(3)`` — which we choose to be given by the usual
  construction via the spatial bivectors, exponentiating ``𝐢 =
  𝐳𝐲``, ``𝐣 = 𝐱𝐳``, and ``𝐤 = 𝐲𝐱``.  Note that this subgroup
  does not involve any factors of ``𝐭``; it is purely spatial.
- **``A``** is the abelian subgroup of boosts in a fixed direction.
  We take this direction to be the z-axis, so that we have ``A =
  \left\{ \exp\left[\tfrac{φₐ}{2} \, 𝐭𝐳 \right] \mid φₐ ∈ ℝ
  \right\}``.
- **``N``** is the (nilpotent, hence the "N") subgroup of null
  rotations, which consists of all Lorentz transformations that can be
  represented as null rotations about a fixed null vector.  We take
  this vector to be the null vector ``\boldsymbol{ℓ}``, so that we
  have ``N = \left\{ \exp\left[\tfrac{1}{2} \boldsymbol{ℓ ξ}\right]
  \mid \boldsymbol{ξ} = ξˣ𝐱 + ξʸ𝐲 \right\}``.  Nilpotency means that
  the generator ``\boldsymbol{ℓ ξ}`` raised to an integer power is
  zero — in this case, ``(\boldsymbol{ℓ ξ})² = 0`` because
  ``\boldsymbol{ℓ}`` and ``\boldsymbol{ξ}`` anticommute and
  ``\boldsymbol{ℓ}² = 0``.

Note that ``φₐ`` here is *not* the rapidity of the overall boost if we
factor a transformation as a boost and a rotation.  Rather, it is the
rapidity of the particular boost picked out as ``A`` by this
construction.  ``N`` also contributes to the overall boost, so the
overall rapidity is not just ``φₐ``.

One nice feature of this decomposition is that we can compute the
factorization of a given Lorentz transformation ``Λ`` fairly simply,
without any transcendental functions — just basic algebra and one
square root.  The method is described in detail in [the `Quaternionic`
documentation](@extref Quaternionic :std:label:`iwasawa-kan`), and
implemented in [`Quaternionic.KAN`](@extref Quaternionic
:jl:function:`Quaternionic.KAN`).

However, for our purposes, we do not actually need the full
factorization; we just need the ``K`` factor, which is somewhat faster
to compute.  This is implemented in the [`aberration`](@ref
Scri.aberration) function, for reasons explained below.  An important
fact is that

```math
𝐑_{AN} = 𝐑_{φₐ𝐳}\, 𝐑_{\boldsymbol{ξ}}
```

leaves the *direction* ``\boldsymbol{ℓ}`` invariant, but just rescales
the vector as

```math
𝐑_{AN} \boldsymbol{ℓ} \bar{𝐑}_{AN} = e^{φₐ} \boldsymbol{ℓ}.
```

## Iwasawa and Hopf

We can further decompose the ``K`` factor ``\mathrm{Spin}(3)`` via the
Hopf fibration — not into sub*groups*, but into sub*spaces*.
``\mathrm{Spin}(3)`` is homeomorphic to ``𝕊³``, and the Hopf
fibration decomposes it into ``𝕊²`` and ``𝕊¹``.  Here, ``𝕊²``
corresponds to points on the null cone, while ``𝕊¹`` corresponds to
rotations about the null direction.

Specifically, given the choice of ``𝐳``, we have the Hopf map

```math
\begin{aligned}
𝔥 &: \mathrm{Spin}(3) \to 𝕊² \\
𝔥 &: 𝐑 \mapsto k̂ = 𝐑 𝐳 𝐑̄.
\end{aligned}
```

A section is a choice of ``𝐑 ∈ \mathrm{Spin}(3)`` for a given point
in ``𝕊²``, serving as a sort of inverse of ``𝔥``.  We choose the
"smallest" ``𝐑`` that rotates ``𝐳`` to the desired direction, if
defined, and make a particular choice for the one case where it is
not, at ``-𝐳``:

```math
\begin{aligned}
σ &: 𝕊² \to \mathrm{Spin}(3) \\
σ &: k̂ \mapsto \begin{cases}
\frac{1 - k̂𝐳}{\sqrt{2 + 2 k̂⋅𝐳}} & \text{if } k̂ \neq -𝐳, \\
\exp\left[\frac{π}{2} 𝐲𝐳 \right] & \text{if } k̂ = -𝐳.
\end{cases}
\end{aligned}
```

Note that ``𝔥 ∘ σ`` is the identity function on ``𝕊²``.  And we can
calculate the fiber (preimage of ``𝔥``) over any ``k̂ ∈ 𝕊²`` as

```math
𝔥⁻¹(k̂) = \left\{ σ(k̂)\, \exp\left[\frac{γ}{2} 𝐱𝐲 \right] \mathrel{\Big|} γ ∈ ℝ \right\}.
```

In this case, we could actually interpret the fiber ``𝕊¹`` as being
``\mathrm{Spin}(2) ≃ U(1)``.  So along with the KAN decomposition, we have decomposed ``\mathrm{Spin}^+(3,1)`` into a product of the subspace ``𝕊²`` and a series of groups:

```math
\begin{gathered}
\mathrm{Spin}^+(3,1)
≅
\bigg\{ k̂ \bigg\} ×
\bigg\{ \exp\bigg[\frac{γ}{2} 𝐱𝐲 \bigg] \bigg\} ×
\bigg\{ \exp\bigg[\frac{φₐ}{2} \, 𝐭𝐳 \bigg] \bigg\} ×
\bigg\{ \exp\bigg[\frac{1}{2} \boldsymbol{ℓ ξ}\bigg] \bigg\} \\
Λ = 𝐑_n\, 𝐑_γ\, 𝐑_{φₐ}\, 𝐑_{\boldsymbol{ξ}}.
\end{gathered}
```

We can also parametrize the ``𝕊²`` via a vector ``\boldsymbol{ζ}`` in
the ``𝐱-𝐲`` plane.  Then ``R_{\boldsymbol{ζ}} =
\exp{\boldsymbol{ζ}𝐳}`` will rotate ``𝐳`` onto a point on the
sphere.  For example, if we write

```math
\boldsymbol{ζ} = 𝐱\frac{θ}{2}\, \exp\left[𝐱𝐲ϕ\right]
```

the rotor ``R = \exp{\boldsymbol{ζ}𝐳}`` carries ``𝐳`` to the point
given by the spherical coordinates ``(θ, ϕ)``, and transports the
standard dyad at the pole to the dyad at that point.  Note that this
parameterization is similar to *but different from* the stereographic
coordinates ``ζ`` used most commonly in the literature.  The most
noticeable difference is that the stereographic coordinates map the
sphere to the entire complex plane, while ``\boldsymbol{ζ}`` maps the
sphere to a disk within a radius of ``π/2``.  Specifically, the factor
of ``θ/2`` above becomes ``\tan(θ/2)`` for stereographic coordinates.
Note that the entire bound of the disk at radius ``π/2`` corresponds
to the single point ``-𝐳`` on the sphere, and the mapping then wraps
back around the sphere, making this mapping fail to be globally
homeomorphic.

Then, the general element of ``\mathrm{Spin}^+(3,1)`` can be written
(not always uniquely) as

```math
Λ = 𝐑_{\boldsymbol{ζ}}\, 𝐑_γ\, 𝐑_{φₐ}\, 𝐑_{\boldsymbol{ξ}}.
```

This decomposition is particularly nice because the latter three
factors preserve the null direction, whereas ``𝐑_{\boldsymbol{ζ}}``
just rotates the null direction.  In fact, we can commute the rotors
so that ``𝐑_{\boldsymbol{ζ}}`` comes *first*, if we reinterpret each
of the other rotors as picking out the null direction as the special
direction, rather than the original ``𝐳`` direction.  In this case,
the spin and boost terms describe rotating about the null direction
and boosting along the null direction — exactly as the spin and boost
weights are defined.

## Reversing the order

The order of operations given above is not the only order that works.
In fact, we can reverse the order of operations, and thereby come
closer to how we actually think of the transformations relevant to an
arbitrary point in the sky.  First, we just get a little more explicit
about the decomposition and write

```math
Λ = 𝐑_{\boldsymbol{ζ}}\, 𝐑_{γ𝐳}\, 𝐑_{φₐ𝐳}\, 𝐑_{\boldsymbol{ξ}},
```

so that we can see the axes about which these transformations are
applied.  Then we permute the order of operations to write

```math
Λ = 𝐑_{\boldsymbol{ξ}'''}\, 𝐑_{φₐ𝐳''}\, 𝐑_{γ𝐳'}\, 𝐑_{\boldsymbol{ζ}},
```

where this is *exactly the same transformation* as before, but written
with different axes.  We *first* rotate the ``𝐳`` axis onto the
direction in which the *final* null ray will appear, then we rotate
about that null ray, followed by a boost along that null ray, and
finally we apply a null rotation about that same null ray.  The
``\boldsymbol{ζ}`` generator of the rotation is exactly the same as
before, and the other three parameters ``\boldsymbol{ξ}``, ``φₐ``,
``γ`` all have the same values as before (where ``\boldsymbol{ℓ ξ}``
is thought of just in terms of components with respect to whichever
basis is operative at the given stage of the transformation), but the
*interpretations* of these parameters have now changed to being with
respect to the final null direction, rather than the original ``𝐳``
direction.

That is, the factor ``𝐑_K = 𝐑_{\boldsymbol{ζ}}\, 𝐑_{γ𝐳}`` produces
the same null direction as the full Lorentz transformation ``Λ``:

```math
𝐑_K \boldsymbol{ℓ} \bar{𝐑}_K ∝ Λ \boldsymbol{ℓ} \bar{Λ}.
```

(The proportionality factor is just ``e^{φₐ}``.)  That is, ``𝐑_K``
produces exactly the same null direction for the unboosted observer as
the full ``KAN`` produces for the boosted observer.

So if we want to know the values of the fields along that null vector
in the *boosted* frame, we need to know the values of the fields in
the *unboosted* frame along the null vector produced by ``𝐑_K``.

That's not quite the full story, however, because [spin-weighted
fields *require*](@cite Boyle_2016) some information about a chosen
tangent direction.  More generally, a function with spin weight ``s``
and boost weight ``b`` is a function on ``\mathrm{Spin}^+(3,1)`` that
satisfies the right-equivariance condition

```math
f\left(Λ\, 𝐑_{γ𝐳}\, 𝐑_{φₐ𝐳}\right)
= e^{i s γ} e^{b φₐ} f(Λ).
```

The null rotation ``𝐑_{\boldsymbol{ξ}}`` cannot satisfy such a simple
equivariance condition; it must be accounted for by considering mixing
of the tetrad as shown above.  But the effect of ``𝐑_{φₐ𝐳}`` can be
accounted for with a simple multiplication by ``e^{b φₐ}`` — which we
will see is essentially the conformal factor ``κ`` to the appropriate
power.  So information about the direction of propagation and tangent
direction needed for a spin-weighted field is entirely contained in
the ``𝐑_K`` factor.

That is, to transform the fields, we first need to evaluate them on
``𝐑_K``, then include the effects of the null rotation and boost.
This is the reason why the [`aberration`](@ref Scri.aberration)
function just computes the ``K`` factor, rather than the full ``KAN``
decomposition produced by [`Quaternionic.KAN`](@extref Quaternionic
:jl:function:`Quaternionic.KAN`).

## Flagpoles and flagplanes

Penrose's picture of a spinor gives the preceding discussion a
concrete geometric vocabulary.  A spinor determines a *flagpole* — a
null vector, here ``\boldsymbol{ℓ} = (𝐭+𝐳)/\sqrt{2}`` at the pole —
and a *flagplane*, a null half-plane containing the flagpole, here
``\boldsymbol{ℓ} ∧ 𝐱``.  A spin-weighted field is really a function
of both: the flagpole is the direction of propagation, and the
flagplane picks out the tangent orientation whose rotation the spin
weight measures.  In the language of the equivariance condition
above, the flagpole is what ``𝐑_{φₐ𝐳}`` rescales (boost weight) and
the flagplane is what ``𝐑_{γ𝐳}`` rotates (spin weight).

The ``AN`` factors preserve both the flagpole's direction and the
flagplane, so all the *positional* information — which flagpole, and
which flagplane through it — is encoded in the ``K`` factor alone.
That is the geometric content of the statement above that ``𝐑_K``
contains everything a spin-weighted field needs: it is the (unique, up
to the equivariance already accounted for) rotation carrying the
reference flagpole and flagplane at the pole onto those produced by
the full transformation ``Λ``.

### Working in the unprimed frame

One further subtlety arises in the implementation.  The grid rotors
``𝐑'_p`` describing the output pixels are *constructed* with
components relative to the primed (transformed) frame, but every
computation in [`transform!`](@ref) is carried out with components
relative to the unprimed frame.  In effect, the code deliberately
"misinterprets" the components of ``𝐑'_p`` as unprimed-frame
components.  What we *want* is the rotor whose generating bivectors
are the primed basis planes — e.g. ``𝐑'_p{}^x\, 𝐳'𝐲'`` — but
building it from those components with the *unprimed* planes
``𝐳𝐲``, etc., yields a different rotor; the two are related by
conjugation with the frame transformation ``Λ = B\,R`` itself,

```math
𝐑_p = Λ\, 𝐑'_p\, Λ̃.
```

Now, the flagpole this conjugated rotor should act on is the *primed*
one, ``\boldsymbol{ℓ}' = Λ\boldsymbol{ℓ}Λ̃``, and conjugation
telescopes:

```math
𝐑_p\, \boldsymbol{ℓ}'\, 𝐑̃_p
= Λ\, 𝐑'_p\, Λ̃\; Λ\boldsymbol{ℓ}Λ̃\; Λ\, 𝐑̃'_p\, Λ̃
= (Λ\, 𝐑'_p)\, \boldsymbol{ℓ}\, \widetilde{(Λ\, 𝐑'_p)}
```

— and the same holds for the flagplane.  So the flagpole and
flagplane we need are simply those produced by the *product*
``Λ\,𝐑'_p`` acting on the reference pair at the pole, and the
pure-rotation factor carrying the reference pair onto them is the
``K`` factor of that product.  This is exactly what
[`aberration`](@ref Scri.aberration)`(R * R′ₚ, v⃗)` computes.  The
``A`` and ``N`` factors that the extraction discards are exactly the
ones accounted for elsewhere — the conformal factor ``κ`` and the
component mixing.

## Reinterpreting ``i``

TL;DR: ``i ∈ ℂ`` is replaced by ``𝐈₃ = 𝐈𝐭``, the spatial
pseudoscalar.  It actually transforms whenever ``𝐭`` transforms, but
if we just write expressions in terms of ``𝐈₃' = 𝐈𝐭'`` without
explicitly transforming ``𝐭``, that should be fine, because the
interpretation of ``i`` also needs to change.  And at that point, it's
just a bookkeeping device, so we don't need to worry about the fact
that it transforms.  The reason ``𝐈₃`` appears is because it is
central in the spatial subalgebra; it commutes with everything, which
is why it can act like ``i``.

The unit imaginary ``i ∈ ℂ`` is a purely algebraic object that has no
geometric meaning to Newman and Penrose.  In Geometric Algebra, we try
to identify the geometric meaning of all algebraic objects.  But the
replacement for ``𝐦`` is not so clear.  We need something that
ensures ``𝐦𝐦=0``, while also transforming reasonably under null
rotations.  The obvious guess is ``i↦𝐈``, which is invariant under
(proper, orthochronous) Lorentz transformations.  Unfortunately,
``(𝐱+𝐈𝐲)²`` simply does not have zero scalar part.  The next
obvious guess is ``i↦𝐱𝐲``, the pseudoscalar of the "screen" space
that ``𝐦`` represents.  Unfortunately, ``(𝐱+𝐱𝐲𝐲)=2𝐱``, which
also obviously does not square to zero.  Finally, we come to
``𝐈₃=𝐈𝐭``.  This does actually work correctly, with the caveat that
``𝐭`` also transforms; when we transform a quantity involving ``i``,
we have to remember that ``i`` will have new meaning in the new frame.

!!! warning "Transformation of 𝐈₃"

    ``𝐈₃`` itself transforms under null rotations, so we have to
    expect our transformation law for ``𝐦`` to reflect this.
    Specifically, we need to factor as
    ``𝐑 𝐈₃𝐲 𝐑̄  = (𝐑 𝐈₃ 𝐑̄ )\, (𝐑 𝐲 𝐑̄)``.

This is *almost* the null tetrad used in, e.g., the Newman-Penrose
formalism, except our definitions of ``𝐦`` and ``𝐦̄`` do not use the
unit imaginary ``i ∈ ℂ``, but rather the unit pseudoscalar ``𝐈 ∈
𝒢(ℝ^{3,1})``.  In fact, with these definitions, ``𝐦`` and ``𝐦̄``
are not even vectors, but more general multivectors.  This makes
almost no difference to the calculations, but it does allow us to work
entirely within the geometric algebra, without the gratuitous and
geometrically meaningless use of complex numbers in just part of the
tetrad.
