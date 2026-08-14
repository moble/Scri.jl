# The Lorentz Group

Here, our focus is on using Geometric Algebra (GA) to understand and
implement Lorentz transformations — and to decompose those
transformations into conventional pieces.  In the context of spacetime
and Lorentz transformations, GA allows us to work with null tetrads
and the Lorentz group ``ℒ`` in a more intuitive and geometrically
meaningful way.  In fact, it will be easier to work with its double
cover ``\mathrm{Spin}⁺(3,1)``, which is the group of rotors in the GA
of 3+1 spacetime.  This treatment ties in neatly with the classical
approach via stereographic coordinates and complex analysis
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
we can apply any Lorentz transformation to it to get the null tetrad
appropriate to any other observer.  In particular, we can apply an
ordinary rotation to "point" ``\boldsymbol{ℓ}`` in any direction we
like.  This suggests the standard factorization of the Lorentz group,
beginning with rotations tied to the spherical coordinates, followed
by rotations about the radial direction, followed by boosts in the
radial direction, and finally followed by null rotations about the
rotated ``\boldsymbol{ℓ}``.  We will discuss these transformations in
more detail below, but the important point is that we decompose a
general Lorentz transformation *specifically with respect to* either
``\boldsymbol{ℓ}`` or ``𝐧``.

## The Lorentz group ``ℒ``

Before we decompose the (proper orthochronous) Lorentz group, we need
to understand how the Lorentz group shows up in Geometric Algebra.
The answer is already implicit in [the primer](@ref "A Primer on
Geometric Algebra"): every proper orthochronous Lorentz transformation
can be expressed as a rotor ``R ∈ \mathrm{Spin}⁺(3,1)`` — the product
of an even number of unit-norm vectors — acting on a vector ``𝐯`` by
conjugation,

```math
𝐯' = R\, 𝐯\, R̃,
```

with composition of transformations given by multiplication of rotors.
Each rotor is the exponential of a bivector, and the *type* of
transformation is fixed by the sign of that bivector's square:
negative for pure rotations, positive for pure boosts, and zero for
null rotations.

A **rotation** comes from a *spatial* bivector, which squares to
``-1``.  Writing ``𝐢 = 𝐳𝐲``, ``𝐣 = 𝐱𝐳``, ``𝐤 = 𝐲𝐱`` for the
three basis spatial bivectors, a rotation by angle ``θ`` in the plane
of a unit bivector ``𝐁`` is

```math
R = \exp\left[\frac{θ}{2} 𝐁\right] = \cos\frac{θ}{2} + 𝐁 \sin\frac{θ}{2}.
```

These involve no factor of ``𝐭`` and generate the maximal compact
subgroup ``\mathrm{Spin}(3)``.  Note that this is *exactly* [Euler's
formula](https://en.wikipedia.org/wiki/Euler's_formula), with the
bivector playing the role of ``i``.  (In fact, complex analysis is
just GA in two dimensions, with ``i = 𝐱𝐲`` or ``i = 𝐲𝐱``.)

A **boost** comes from a *timelike* bivector, which squares to ``+1``.
A boost along ``𝐳`` with rapidity ``φ`` is

```math
B = \exp\left[\frac{φ}{2} 𝐭𝐳\right] = \cosh\frac{φ}{2} + 𝐭𝐳 \sinh\frac{φ}{2},
```

the trigonometric functions turned hyperbolic by ``(𝐭𝐳)² = +1``.
The speed is ``β = \tanh φ`` and the Lorentz factor is ``γ = \cosh
φ``.  This is *similar to* Euler's formula, but with hyperbolic trig
functions replacing the circular ones because the bivector squares to
``+1`` instead of ``-1``.

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
the null rotation *simultaneously* boosts in a direction orthogonal to
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
``\boldsymbol{ξ}``-``𝐳`` plane.  To see this, we expand the product

```math
\boldsymbol{ℓ ξ} = \frac{𝐭\boldsymbol{ξ} + 𝐳\boldsymbol{ξ}}{\sqrt{2}}.
```

The ``𝐭\boldsymbol{ξ}`` term generates a boost in the
``\boldsymbol{ξ}`` direction, and the ``𝐳\boldsymbol{ξ}`` term
generates the rotation in the ``\boldsymbol{ξ}``-``𝐳`` plane that
counteracts the mixing that would give ``\boldsymbol{ℓ}`` a component
along ``\boldsymbol{ξ}``.

Define the rotor

```math
R = \exp\left[ \frac{1}{2} \boldsymbol{ℓ ξ} \right].
```

As usual, the exponential is evaluated by way of its series expansion,
but in this case the series terminates after the second term because
the generator squares to zero because ``\boldsymbol{ℓ}² = 0``.
Specifically, ``\boldsymbol{ℓ}`` and ``\boldsymbol{ξ}`` are orthogonal
so they anticommute, and we can calculate the square of the generator
that would appear in the second term:

```math
\left(\boldsymbol{ℓ ξ}\right)^2
=\boldsymbol{ℓ ξ ℓ ξ}
=-\boldsymbol{ℓ ℓ ξ ξ}
=-\boldsymbol{ℓ}^2 \boldsymbol{ξ}^2
= 0.
```

Therefore, we have

```math
R = 1 + \frac{1}{2} \boldsymbol{ℓ ξ},
```

which makes calculations particularly simple.  The fact that the
generator of this transformation ``\tfrac{1}{2} \boldsymbol{ℓ ξ}``
squares to zero makes it "nilpotent", which means that the group of
null rotations about ``\boldsymbol{ℓ}`` is a nilpotent Lie group.
This becomes important below, where the nilpotent property gives this
group its distinction as the "N" of the "KAN" decomposition.

It can be helpful to see that the null rotations about a given null
vector form a group.  Let ``\boldsymbol{ξ}_1`` and
``\boldsymbol{ξ}_2`` be two vectors in the ``𝐱``-``𝐲`` plane.  We
have

```math
\begin{aligned}
R_{\boldsymbol{ξ}_2} R_{\boldsymbol{ξ}_1}
&=
\left(1 + \frac{1}{2} \boldsymbol{ℓ ξ}_2\right)
\left(1 + \frac{1}{2} \boldsymbol{ℓ ξ}_1\right)
\\ &=
1 + \frac{1}{2} \boldsymbol{ℓ}
\left( \boldsymbol{ξ}_2 + \boldsymbol{ξ}_1 \right)
+ \frac{1}{4} \boldsymbol{ℓ ξ}_2 \boldsymbol{ℓ ξ}_1
\\ &=
1 + \frac{1}{2} \boldsymbol{ℓ}
\left( \boldsymbol{ξ}_2 + \boldsymbol{ξ}_1 \right)
\\ &=
R_{\boldsymbol{ξ}_2 + \boldsymbol{ξ}_1}.
\end{aligned}
```

That is, the multiplicative group of null rotations about *a given
null vector* is isomorphic to the additive group of vectors in the
orthogonal plane.

The action on the basis ``(𝐭, 𝐱, 𝐲, 𝐳)`` is not enlightening, but
the results on the ``(\boldsymbol{ℓ}, 𝐦, 𝐦̄, 𝐧)`` basis are nicely
systematic:

```math
\begin{aligned}
R \boldsymbol{ℓ} R̄ &= \boldsymbol{ℓ}, \\
R 𝐦 R̄ &= 𝐦 + 𝐱\boldsymbol{ξℓ}, \\
R 𝐦̄ R̄ &= 𝐦̄ - 𝐱\boldsymbol{ξℓ}, \\
R 𝐧 R̄ &= 𝐧 + \boldsymbol{ξ} + \frac{1}{2} ξ² \boldsymbol{ℓ}.
\end{aligned}
```

This is just the usual conformal transformation of the tetrad, though
exhibited in a simpler and more geometric form.  (Note that the
factors of ``𝐱`` in the ``𝐦`` and ``𝐦̄`` expressions essentially
make their terms "complex", like the ``𝐦`` and ``𝐦̄`` objects
themselves — meaning that those terms are sums of vectors and
trivectors.)

!!! details "Conformal representations and the stereographic approach"

    In general, it is possible to represent the geometry of a normed
    vector space ``ℝ^{p,q}`` using the "conformal representation",
    which introduces two new vectors, one timelike and one spacelike.
    Points in ``ℝ^{p,q}`` are represented as null vectors in
    ``ℝ^{p+1,q+1}``.  The value of adding these extra dimensions is
    that *extended* objects like lines, planes, hyperplanes, circles,
    and higher-dimensional spheres can all be represented by objects
    in the geometric algebra of this higher-dimensional space.
    Intersections, tangencies, and transformations of these extended
    objects are then obtained as simple algebraic operations using
    geometric algebra.

    In fact, the field of "Conformal Geometric Algebra" is built on
    this fact and the facility with which GA represents its
    transformations [DoranLasenby_2003](@cite), and so has become a
    standard tool in computer graphics and robotics for describing
    transformations of the plane and 3D space.

    A closely related fact is that the conformal group of ``ℝ^{p,q}``
    is isomorphic to the connected component of the identity of the
    group ``\mathrm{SO}(p+1,q+1)``.  Though generally only true for
    ``p+q > 2``, the case ``(p,q)=(2,0)`` qualifies with some caveats
    [Schottenloher_2008; Theorems 2.9 and 2.11](@cite).  In our case,
    ``ℝ^{2,0}`` is the  ``𝐱``-``𝐲`` plane, and this says that the
    conformal group of the plane is isomorphic to
    ``\mathrm{SO}^+(3,1)``, which is exactly the Lorentz group!  But
    complex algebra is precisely the geometric algebra of the
    ``𝐱``-``𝐲`` plane, so complex analysis is exactly what we need to
    analyze the Lorentz group in this approach.

    In *most* applications of these facts, ``ℝ^{p,q}`` is generally
    the physical space of interest, and the extra two dimensions are
    essentially fictitious — invented to take advantage of this
    mathematical trick.  What's remarkable about the stereographic /
    complex approach pioneered by Penrose is that it *reverses* this
    picture.  Null vectors in the full physical spacetime are the
    objects of interest and Lorentz is the crucial group, whereas the
    complex plane and its conformal transformations are mathematical
    artifice.  And yet, with help from complex analysis, this plane
    becomes the primary analysis tool for the null cone and for the
    Lorentz group in the BMS literature.

    Geometric Algebra actually bridges both pictures; it allows us to
    directly manipulate both general Lorentz transformations and the
    complex representation using the same language and notation.

Nothing in this construction relies on our particular choice of
``\boldsymbol{ℓ}``.  In particular, exchanging ``\boldsymbol{ℓ} ↔ 𝐧``
throughout yields the family of null rotations that fix ``𝐧``
instead, shifting ``𝐦`` along ``𝐧`` and sending ``\boldsymbol{ℓ} ↦
\boldsymbol{ℓ} + \boldsymbol{ξ} + \tfrac{1}{2} ξ²\, 𝐧``.  When the
distinction matters we subscript each family by the null vector it
fixes — ``N_{\boldsymbol{ℓ}}`` versus ``N_{𝐧}``.  The question of
*which* family acts on radiation data at null infinity turns out to be
surprisingly delicate; see ["Which null direction?"](@ref
which_null_direction) below.

## Cartan (polar) decomposition

There are two natural ways to factor a general rotor ``Λ ∈
\mathrm{Spin}⁺(3,1)``, and we will use both.  The first is somewhat
simpler, factoring an arbitrary transformation into a boost followed
by a rotation:

```math
Λ = R\, B.
```

Equivalently, we could factor it as a rotation followed by a boost,
though the boost direction would be rotated by ``R``:

```math
Λ = R\, B = R\, B\, R⁻¹\, R = \left(R\, B\, R⁻¹\right)\, R = B'\, R.
```

This is the (global) *Cartan decomposition*, which — in general —
takes an arbitrary semisimple Lie group ``G`` and factors it into the
product of a compact subgroup ``K`` and non-compact part ``P`` so that
``G = K\, P``.  Here, the compact subgroup is the rotation group
``\mathrm{Spin}(3)``, and the non-compact part is the set of boosts —
which is not a subgroup, because of course it is not closed under
multiplication/composition.

Obviously, the decomposition is not unique, but is determined by a
choice of time axis ``𝐭``, which allows us to specify what is meant
by a "spatial" rotation.  We can understand the decomposition at the
simplest level by working backwards: start with a rotation and a
boost, and see what their product looks like.  Any boost is specified
by a spatial direction ``v̂`` and rapidity  ``φ``.  The corresponding
rotor is

```math
B = \exp\left[\tfrac{φ}{2} 𝐭v̂\right]
= \cosh\frac{φ}{2} + 𝐭v̂ \sinh\frac{φ}{2}.
```

We don't need any specific form for the rotation, except to note that
it is generated exclusively by spatial components so its rotor
contains only scalar, ``𝐲𝐱``, ``𝐱𝐳``, and ``𝐳𝐲`` components —
precisely the space of the quaternions.  We can expand

```math
Λ = R\, B = R \cosh\frac{φ}{2} + R 𝐭 v̂ \sinh\frac{φ}{2}.
```

Note that the ``R 𝐭 v̂`` contains *only* terms with factors of ``𝐭``
since neither ``R`` nor ``v̂`` contains another factor of ``𝐭`` for
it to cancel with.  That is, the second term is a linear combination
of ``𝐭𝐱``, ``𝐭𝐲``, ``𝐭𝐳``, and the pseudoscalar — precisely the
*dual space* of the quaternions.  Therefore, we can simply separate
the terms with a factor of ``𝐭`` from those without, leaving us with
``R \cosh\tfrac{φ}{2}``.  We can then simply normalize that result to
get ``R``.  Finally, we can just multiply on the left by ``R⁻¹ = R̃``
to find ``B = R̃Λ``.

We can be a little more formal about the separation of a Lorentz
spinor into its quaternionic and dual-quaternionic parts by looking at
the effects of conjugation by ``𝐭``:

```math
Θ(Λ) = 𝐭 Λ 𝐭⁻¹.
```

If we examine the effect on the decomposed ``Λ = R B``, this does not
affect the rotation, but it *does* flip the sign of the boost.  In
fact, this operation is important enough that it has a name: the
**global Cartan involution** [Knapp_1996](@cite).  For now, the main
fact that we need is that it is an involution — meaning that ``Θ(Θ(Λ))
= Λ`` — so we can use it to split a Lorentz spinor into the parts that
are invariant and the parts that change sign under ``Θ``, which we
then use to define the "real" and "imaginary" parts:[^1]

```math
\begin{aligned}
\Re_𝐭(Λ) &= \frac{Λ + 𝐭 Λ 𝐭⁻¹}{2}, \\
\Im_𝐭(Λ) &= \frac{Λ - 𝐭 Λ 𝐭⁻¹}{2𝐈}.
\end{aligned}
```

Here, we write the pseudoscalar ``𝐈 = 𝐭𝐱𝐲𝐳``.  (Division by
``𝐈`` is just multiplication by ``𝐈⁻¹ = -𝐈``, which happens to
commute with all spinors, so the fraction notation is unambiguous in
this case; it is used here to emphasize the similarity to the complex
definitions.)  With this notation, the factorization can be written
clearly:

```math
R = \frac{\Re_𝐭\, Λ}{\lVert \Re_𝐭\, Λ \rVert},
\qquad
B = R̃\, Λ.
```

The computation is cheap, with no transcendental functions or SVD as
is seen in matrix methods — just a single real square root and simple
algebraic operations.  This is implemented in the functions
[`Quaternionic.RB`](@extref Quaternionic
:jl:function:`Quaternionic.RB`) and [`Quaternionic.BR`](@extref
Quaternionic :jl:function:`Quaternionic.BR`).

[^1]: This is treated very nicely and using geometric algebra — though
    with a flavor closer to that found in quantum mechanics — by
    [FrancisKosowsky_2005](@citet).  Incidentally, they have another
    paper [FrancisKosowsky_2004](@cite) that isn't specifically about
    spinors, but has some overlap and would likely be quite
    interesting to readers of this page.

!!! details "Complexified quaternions and Quaternionic.jl"

    In the complexified-quaternion representation of spinors (which is
    how [`Quaternionic.Lorentz`](@extref Quaternionic
    :jl:type:`Quaternionic.Lorentz`) implements the representation),
    the choice of ``𝐭`` is implicit in the identification of the real
    and imaginary parts.  The components with a ``𝐭`` factor are
    precisely the components proportional to the unit imaginary ``i``,
    which functions *exactly* as the pseudoscalar ``𝐈`` and so can be
    seen as simply a different notation for it.  However, this also
    means that the real and imaginary parts defined above by ``\Re_𝐭``
    and ``\Im_𝐭`` can be computed by extracting the usual complex real
    and imaginary parts of the four complex components.  Those
    operations are accomplished with
    [`Quaternionic.ℂreal`](@extref
    Quaternionic :jl:function:`Quaternionic.ℂreal`) and
    [`Quaternionic.ℂimag`](@extref
    Quaternionic :jl:function:`Quaternionic.ℂimag`).

We can write the Cartan decomposition as

```math
\mathrm{Spin}^+(3,1) \mathrel{\cong_{\text{Top}}} K \times P,
```

where ``K = \mathrm{Spin}(3)`` is the (maximal compact) subgroup of
rotations, and ``P`` is the set of boosts.  Note that this
relationship is not *group isomorphism*
``\mathrel{\cong_{\text{Grp}}}``, because the ``P`` factor is not a
group; rather, it is just *topological isomorphism* (homeomorphism)
``\mathrel{\cong_{\text{Top}}}``.  But the point is that any Lorentz
transformation can be written as the product of a rotor and a boost.
The decomposition is not unique, but is determined by a choice of time
axis ``𝐭``.

It is possible to go further by *also* choosing a preferred spatial
direction — or equivalently a preferred null direction.  First, we
will proceed with the "Iwasawa" decomposition, which splits the ``P``
factor into two groups — ``A`` representing boosts along the spatial
direction and ``N`` representing null rotations about that null
direction.  The result is a ``KAN`` decomposition, where the ``K``
factor is the same as in the Cartan decomposition.  (However, because
the ``N`` factor is simultaneously a boost and a rotation, the
particular element of ``K`` picked out for a given ``Λ`` will differ
between the two decompositions.)  The selection of a preferred spatial
direction also lets us select a subgroup of rotations that leave that
direction invariant, which is the "Hopf-Levi" decomposition — giving
us what we might call the ``SMAN`` decomposition.

## Iwasawa's ``KAN`` decomposition

Iwasawa showed how to decompose any semisimple Lie group
[Knapp_1996](@cite) (which includes the Lorentz group and its double
cover ``\mathrm{Spin}^+(3,1)``) into three subgroups:

```math
G = KAN.
```

In the case of the Lorentz and ``\mathrm{Spin}^+(3,1)`` groups, this
decomposition is specified by a choice of a preferred time axis and a
preferred spatial or null direction.  For definiteness, we choose the
time axis to be ``𝐭`` and the spatial direction to be ``𝐳``, which
corresponds to the null direction ``\boldsymbol{ℓ}`` — though the
procedure is the same for any other choices of orthogonal timelike and
spacelike unit vectors.

We use those to pick out boosts along that spatial direction, and null
rotations that fix the corresponding null vector.  The decomposition
results in the following subgroups:

- **``K``** is the maximal compact subgroup ("kompakt" in German,
  hence the "K"), which in the case of the Lorentz group is the
  rotation group ``\mathrm{SO}(3)``, and in the case of
  ``\mathrm{Spin}(3,1)`` is ``\mathrm{Spin}(3)`` — which we choose to
  be given by the usual construction via the spatial bivectors,
  exponentiating ``𝐢 = 𝐳𝐲``, ``𝐣 = 𝐱𝐳``, and ``𝐤 = 𝐲𝐱``.
  Note that this subgroup does not involve any factors of ``𝐭``; it
  is purely spatial.
- **``A``** is the (abelian, hence the "A") subgroup of boosts along
  the z-axis, so that we have ``A = \left\{ \exp\left[\tfrac{φₐ}{2} \,
  𝐭𝐳 \right] \mid φₐ ∈ ℝ \right\}``.
- **``N``** is the (nilpotent, hence the "N") subgroup of null
  rotations about the null vector ``\boldsymbol{ℓ}``, so that we have
  ``N = \left\{ \exp\left[\tfrac{1}{2} \boldsymbol{ℓ ξ}\right] \mid
  \boldsymbol{ξ} = ξˣ𝐱 + ξʸ𝐲 \right\}``.  Nilpotency means that the
  generator ``\boldsymbol{ℓ ξ}`` raised to an integer power is zero —
  in this case, ``(\boldsymbol{ℓ ξ})² = 0`` because ``\boldsymbol{ℓ}``
  and ``\boldsymbol{ξ}`` anticommute and ``\boldsymbol{ℓ}² = 0``.

Note that ``φₐ`` here is *not* the rapidity of the overall boost if we
factor a transformation as a boost and a rotation.  Rather, it is the
rapidity of the particular boost picked out as ``A`` by this
construction.  ``N`` also contributes to the overall boost, so the
overall rapidity is not just ``φₐ``.  Again, it is important to note
that the ``N`` factor is not a pure boost; it incorporates some
spatial rotation as well, which is why the particular elements of
``K`` picked out by the Cartan and Iwasawa decompositions are not the
same, even though they are both elements of the same group.

It will be important to note that the choice of ``𝐳`` instead of
``-𝐳`` is essentially what lets us pick out a particular ``N``
factor.  The opposite choice of sign[^2] would have picked out the
null rotations that fix ``Θ(\boldsymbol{ℓ}) = 𝐧`` instead.  The ``K``
and ``A`` factors would remain precisely the same, but the ``N``
factor would be ``\left\{ \exp\left[\tfrac{1}{2} \boldsymbol{𝐧
ξ}\right] \mid \boldsymbol{ξ} = ξˣ𝐱 + ξʸ𝐲 \right\}`` — and again,
the particular values picked out for a given transformation would
differ.  But that is another legitimate choice, which we wil come back
to [below](@ref which_null_direction).

[^2]: Formally, the decomposition is determined by a choice of
    "positivity" of the "restricted roots" of the full Lie algebra
    ``𝔰𝔭𝔦𝔫(3,1)`` relative to ``𝔞₀``.  Such a restricted root is
    a nonzero element of the dual space, ``λ ∈ 𝔞₀^*`` such that there
    exist some nonzero ``b ∈ 𝔰𝔭𝔦𝔫(3,1)`` with

    ```math
    [a, b] = λ(a)b \quad \text{for all} \quad a∈𝔞₀.
    ```

    But since ``𝔞₀`` is one-dimensional, and this equation is linear
    in ``a``, this is equivalent to

    ```math
    [𝐭𝐳, b] = λ(𝐭𝐳)b.
    ```

    Moreover, ``𝔞₀^*`` is also one-dimensional, so we can just think
    of each ``λ`` as a real number in this equation.

    Therefore, this is just an ordinary eigenvalue problem over the
    6-dimensional space of bivectors.  We can readily compute the
    solutions: the eigenbivectors with eigenvalue ``0`` are the
    bivectors that commute with ``𝐭𝐳``, namely ``𝐲𝐱`` and ``𝐭𝐳``
    (which make up the Levi factor); ``+2`` provides
    ``\boldsymbol{ℓ}𝐱`` and ``\boldsymbol{ℓ}𝐲``; and ``-2`` provides
    ``𝐧𝐱`` and ``𝐧𝐲``.  A choice of positivity then picks out
    either of the latter two pairs, and ``N`` is generated by linear
    combinations of those two bivectors.

One nice feature of this decomposition is that we can compute the
factorization of a given Lorentz transformation ``Λ`` fairly simply,
without any transcendental functions.  Much like the Cartan
decomposition, it involves just basic algebra and one square root.
The method is described in detail in [the `Quaternionic`
documentation](@extref Quaternionic :std:label:`iwasawa-kan`), and
implemented in [`Quaternionic.KAN`](@extref Quaternionic
:jl:function:`Quaternionic.KAN`).

However, for some of our specific applications, we do not always need
the full factorization; we just need the ``K`` factor, which is
slightly simpler to compute.  We explain below why this is all we
need, but the crucial detail is that

```math
𝐑_{AN} = 𝐑_{φₐ𝐳}\, 𝐑_{\boldsymbol{ξ}}
```

leaves the *direction* ``\boldsymbol{ℓ}`` invariant, but just rescales
the vector as

```math
𝐑_{AN} \boldsymbol{ℓ} \bar{𝐑}_{AN} = e^{φₐ} \boldsymbol{ℓ}.
```

The basic point is that the transformed field along ``\boldsymbol{ℓ}``
is determined by the *original* field along ``\boldsymbol{ℓ}``, so we
just need the ``K`` factor to determine that null direction.  Of
course, we then have to use that field *along with* the ``A`` and
``N`` factors to fully determine the transformed field, which is why
this explanation is not yet complete.  We will fill in those details
below.

## [Iwasawa × Hopf-Levi](@id IwasawaHopfLevi)

We can further decompose Iwasawa's ``K`` factor ``\mathrm{Spin}(3)``
via the Hopf fibration — not into sub*groups*, but into sub*spaces*.
``\mathrm{Spin}(3)`` is homeomorphic to ``𝕊³``, and the Hopf
fibration decomposes it into a bundle over ``𝕊²`` with ``𝕊¹``
fibers.  Here, ``𝕊²`` corresponds to points on the null cone — that
is directions on the sky — while ``𝕊¹`` corresponds to rotations
*about* that null direction.  This is a nice geometric picture with
simple calculations making it concrete.  But we can also take a more
formal algebraic approach and arrive at the same result using the
"Levi decomposition" [Knapp_1996](@cite).

### The Hopf fibration

Given the choice of ``𝐳``, we have the Hopf map

```math
\begin{aligned}
𝔥 &: \mathrm{Spin}(3) \to 𝕊² \\
𝔥 &: 𝐑 \mapsto k̂ = 𝐑 𝐳 𝐑̄.
\end{aligned}
```

A section is a mapping producing a choice of ``𝐑 ∈ \mathrm{Spin}(3)``
for a given point in ``𝕊²``, serving as a sort of inverse of ``𝔥``.
We choose the "smallest" ``𝐑`` that rotates ``𝐳`` to the desired
direction, if defined, and make a particular choice for the one case
where it is not, at ``-𝐳``:

```math
\begin{aligned}
σ &: 𝕊² \to \mathrm{Spin}(3) \\
σ &: k̂ \mapsto \begin{cases}
\frac{1 - k̂𝐳}{\sqrt{2 + 2 k̂⋅𝐳}} & \text{if } k̂ \neq -𝐳, \\
\exp\left[\frac{π}{2} 𝐳𝐲 \right] & \text{if } k̂ = -𝐳.
\end{cases}
\end{aligned}
```

Note that ``𝔥 ∘ σ`` is the identity function on ``𝕊²`` — which is
the main defining feature of a section ``σ``.  But it is obviously not
continuous at ``-𝐳``, so it is not a *global* section.

We can also calculate the fiber (preimage of ``𝔥``) over any ``k̂ ∈
𝕊²`` as

```math
𝔥⁻¹(k̂) = \left\{ σ(k̂)\, \exp\left[\frac{γ}{2} 𝐲𝐱 \right]
  \mathrel{\Big|} γ ∈ ℝ \right\}.
```

In this case, we can actually interpret the fiber ``𝕊¹`` as being the
circle group ``\mathrm{Spin}(2) \mathrel{\cong_{\text{Grp}}} U(1)``.
So along with the KAN decomposition, we have decomposed
``\mathrm{Spin}^+(3,1)`` into a product of the *subspace* ``𝕊²`` and
three *subgroups*:

```math
\begin{gathered}
\mathrm{Spin}^+(3,1)
\mathrel{\cong_{\text{Top}}}
\bigg\{ k̂ \bigg\} ×
\bigg\{ \exp\bigg[\frac{γ}{2} 𝐲𝐱 \bigg] \bigg\} ×
\bigg\{ \exp\bigg[\frac{φₐ}{2} \, 𝐭𝐳 \bigg] \bigg\} ×
\bigg\{ \exp\bigg[\frac{1}{2} \boldsymbol{ℓ ξ}\bigg] \bigg\} \\[15pt]
Λ = 𝐑_{k̂}\, 𝐑_γ\, 𝐑_{φₐ}\, 𝐑_{\boldsymbol{ξ}}.
\end{gathered}
```

Again, this relationship is not *group isomorphism*
``\mathrel{\cong_{\text{Grp}}}``, because the ``k̂`` factor is not a
group; rather, it is just *topological isomorphism* (homeomorphism)
``\mathrel{\cong_{\text{Top}}}``.

### The Levi factor

The Levi factor is a formal algebraic approach to the same
decomposition.  Though it is a more general construction, here we will
start from the ``KAN`` decomposition.  Recall that the ``A`` factor is
the one-dimensional group of boosts along the ``𝐳`` axis, so its Lie
algebra ``𝔞₀`` is just the one-dimensional vector space consisting of
scalar multiples of the generator ``𝐭𝐳``.  We define the
"centralizer" of ``𝔞₀`` in ``K`` — the set of elements of ``K`` that
commute with all elements of ``𝔞₀``.  Conventionally denoted ``M``,
that is

```math
M = Z_K(𝔞₀) = \left\{ k∈K \mathrel{\big|} k𝐭𝐳k^{-1} = 𝐭𝐳 \right\}.
```

Here, ``k ∈ \mathrm{Spin}(3)`` commutes with ``𝐭`` — essentially by
definition — so this is equivalent to the set of elements of
``\mathrm{Spin}(3)`` that commute with ``𝐳``.  We can express this in
the Lie algebra ``𝔨₀`` (the generators of ``\mathrm{Spin}(3)``) as

```math
𝔪₀ = \left\{ B∈𝔨₀ \mathrel{\big|} B𝐳B^{-1} = 𝐳 \right\},
```

which we can solve in terms of the basis bivectors to find that it is
the one-dimensional vector space spanned by ``\{𝐲𝐱\}``.  That is, we
have

```math
M = \bigg\{ \exp\bigg[\frac{γ}{2} 𝐲𝐱 \bigg] \bigg\},
```

just as we found in the Hopf picture.  The Levi factor *per se* is
the product of ``M`` and ``A``.

### Parameterizing the ``𝕊²`` factor

Note that there are several ways to parametrize the ``𝕊²`` factor,
and the corresponding value of ``γ`` for each ``k̂`` depends on that
choice.  Essentially, these are different choices of sections.  The
most familiar choice is to use the spherical coordinates ``(θ, ϕ)`` of
the point on the sphere, and then we can write

```math
𝐑_{θ,ϕ} = \exp\left[\frac{ϕ}{2}𝐲𝐱\right]
\exp\left[\frac{θ}{2}𝐱𝐳\right].
```

The points with ``θ = 0`` and ``θ = π`` are singular, so there is a
degeneracy between the values of ``ϕ`` and ``γ`` at those points —
though together they are still capable of representing any element of
the group.  Another common choice is to use stereographic coordinates,
which are complex numbers related (in the conventions of
[PenroseRindler_1984](@citet) for the future null sphere) to the
spherical coordinates by

```math
ζ = \cot\frac{θ}{2}\, \exp\left[ϕ 𝐲𝐱\right],
```

which maps very neatly to a rotor as

```math
𝐑_{ζ} = \frac{ζ + 𝐱𝐳}{\sqrt{1 + |ζ|²}}.
```

This parameterization removes the singularity at ``θ = π``, but ``ζ →
∞`` as ``θ → 0``, so it is also not without its problems — as any
mapping of a sphere to a plane must be.[^3]

[^3]: We could also parametrize the ``𝕊²`` via a vector
    ``\boldsymbol{ζ}`` in the ``𝐱-𝐲`` plane.  Then
    ``R_{\boldsymbol{ζ}} = \exp{\boldsymbol{ζ}𝐳}`` will rotate ``𝐳``
    onto a point on the sphere.  For example, if we write

    ```math
    \boldsymbol{ζ} = \frac{θ}{2}\, \exp\left[ϕ𝐲𝐱\right]𝐱
    ```

    the rotor ``R = \exp{\boldsymbol{ζ}𝐳}`` carries ``𝐳`` to the
    point given by the spherical coordinates ``(θ, ϕ)``, and
    transports the standard dyad at the pole to the dyad at that
    point.  Note that this parameterization is similar to *but
    different from* the stereographic coordinates ``ζ`` used most
    commonly in the literature.  The most noticeable difference is
    that the stereographic coordinates map the sphere to the entire
    complex plane, while ``\boldsymbol{ζ}`` maps the sphere to a disk
    within a radius of ``π/2``.  Specifically, the factor of ``θ/2``
    above becomes ``\tan(θ/2)`` for stereographic coordinates.

    Note that the entire bound of the disk at radius ``π/2``
    corresponds to the single point ``-𝐳`` on the sphere, and the
    mapping then wraps back around the sphere, making this mapping
    fail to be globally homeomorphic.  This approach trades the
    familiar conformal structure and complex analysis for a more
    direct relationship to the geometry and group structure.

One important point is that ``R_{θ,ϕ}`` and ``R_{ζ}`` both map ``𝐳``
to the same point on the sphere, but they *do not* map ``𝐱`` or
``𝐲`` to the same tangent directions at that point.  The difference
is accounted for by the ``γ`` factor, which rotates the tangent
directions at that point.  Nonetheless, regardless of the conventions
chosen, we can always find some combination that produces the desired
overall rotation.  So we will just write the ``𝕊²`` generically as
``𝐑_{k̂}``, remembering that it depends on a choice of section, as does
the corresponding value of ``γ``.

To summarize, the Cartan decomposition as ``KP`` splits by Iwasawa's
decomposition as ``KAN``, which can be further extended by the
Hopf-Levi decomposition to ``SMAN``:

```math
\mathrm{Spin}^+(3,1)
\mathrel{\cong_{\text{Top}}}
S \times M \times A \times N,
```

and the general element of ``\mathrm{Spin}^+(3,1)`` can be written
(not uniquely) as

```math
Λ = 𝐑_{k̂}\, 𝐑_{γ𝐲𝐱}\, 𝐑_{φₐ𝐭𝐳}\, 𝐑_{\boldsymbol{ℓ ξ}}.
```

This decomposition is particularly nice because the latter three
factors preserve the null direction, whereas ``𝐑_{k̂}`` just rotates
the reference tetrad to the chosen null direction.  The spin and boost
weights are related to the middle two factors — the effect of which
can be seen more clearly if we commute the factors, reversing the
order entirely.

## Reversing the order

The order of operations given above is not the only order that works.
In fact, we can reverse the order of operations, and thereby come
closer to how we actually think of the transformations relevant to an
arbitrary point in the sky.  In the decomposition of the previous
equation, the subscripts on the three latter factors indicate (twice)
the bivector generators of the transformations, which are
exponentiated to produce the full rotors.  An important fact about the
exponential is that it commutes with conjugation by another element of
the algebra — meaning that if we conjugate an exponential, the result
is the same as if we had conjugated the generator first and then
exponentiated.  For example, we have

```math
\begin{aligned}
𝐑_{k̂} \, 𝐑_{γ𝐲𝐱}
&= 𝐑_{k̂} \, \exp\left[\frac{γ}{2}𝐲𝐱\right]\, 𝐑_{k̂}^{-1}\, 𝐑_{k̂} \\
&= \exp\left[\frac{γ}{2}𝐑_{k̂} \, 𝐲𝐱\, 𝐑_{k̂}^{-1}\right]\, 𝐑_{k̂} \\
&= \exp\left[\frac{γ}{2}\left(𝐑_{k̂} \, 𝐲\, 𝐑_{k̂}^{-1}\right)
  \left(𝐑_{k̂} \, 𝐱\, 𝐑_{k̂}^{-1}\right)\right]\, 𝐑_{k̂} \\
&= \exp\left[\frac{γ}{2}𝐲'𝐱'\right]\, 𝐑_{k̂} \\
&= 𝐑_{γ𝐲'𝐱'}\, 𝐑_{k̂}.
\end{aligned}
```

That is, we can move the factor on the left to the right side if we
also use that factor to act on the generator of the term it is moving
past.  Extending this idea, we can permute the order of *all* of the
operations to write

```math
Λ = 𝐑_{\boldsymbol{ℓ}''' \boldsymbol{ξ}'''}\, 𝐑_{φₐ𝐭''𝐳''}\,
𝐑_{γ𝐲'𝐱'}\, 𝐑_{k̂},
```

where this is *exactly the same transformation* as before, but written
with different axes.  We *first* rotate the ``𝐳`` axis onto the
direction in which the *final* null ray will appear, then we rotate
about that null ray, followed by a boost along that null ray, and
finally we apply a null rotation about that same null ray.

The generator of the ``R_{k̂}`` rotation is exactly the same as
before, and the other three parameters ``\boldsymbol{ξ}``, ``φₐ``,
``γ`` all have the same values as before (where ``\boldsymbol{ℓ ξ}``
is thought of just in terms of components with respect to whichever
basis is operative at the given stage of the transformation), but the
*interpretations* of these parameters have now changed to being with
respect to the final null direction, rather than the original ``𝐳``
direction.  That is, the factor ``𝐑_K = 𝐑_{k̂}\, 𝐑_{γ𝐲𝐱} =
𝐑_{γ𝐲'𝐱'}\, 𝐑_{k̂}`` produces the same null *direction* as the
full transformation ``Λ``:

```math
𝐑_K \boldsymbol{ℓ} \bar{𝐑}_K ∝ Λ \boldsymbol{ℓ} \bar{Λ}.
```

(The proportionality factor is just ``e^{-φₐ}``.)  So ``𝐑_K``
produces exactly the same null direction for the unboosted observer as
the full ``KAN`` produces for the boosted observer.  If we want to
know the values of the fields along that null vector in the *boosted*
frame, we need to know the values of the fields in the *unboosted*
frame along the null vector produced by ``𝐑_K``.

!!! tip "Telescopes and polarizers"

    The natural way to think about this is in terms of pointing a
    telescope: the ``K`` factor rotates the telescope to point in the
    right direction, while the ``AN`` factors adjust the boost along
    that direction, without actually changing the direction of the
    null vector being observed.  That's not quite the full story,
    however, because our fields are polarized, and the polarization is
    defined with respect to a tangent direction that is also affected
    by ``K``.  Specifically, [spin-weighted fields *require*](@cite
    Boyle_2016) some information about a chosen tangent direction.

    In this view of things, the ``𝕊²`` factor of ``K``, ``𝐑_{k̂}``, is
    responsible for pointing the telescope in the right direction,
    while the ``𝕊¹`` factor ``𝐑_{γ𝐲'𝐱'}`` rotates the polarizer.  We
    can formalize this a little more by establishing a "home" position
    for the telescope and its polarizer — conventionally with the
    telescope pointing along the ``𝐳`` axis and the polarizer aligned
    with the ``𝐱`` axis, both in the rest frame of the original
    coordinate system.

    That, of course, is the picture for past null infinity, where the
    null rays are incoming to the telescope.  For future null
    infinity, we think of an emitting source sending out polarized
    signals in the particular direction.  The idea is fundamentally
    identical, but perhaps a little less familiar.

More generally, a function with spin weight ``s`` and boost weight
``b`` is a function on ``\mathrm{Spin}^+(3,1)`` that satisfies the
right-equivariance condition[^4]

```math
f\left(Λ\, 𝐑_{γ𝐲𝐱}\, 𝐑_{φₐ𝐭𝐳}\right)
= e^{-i s γ} e^{b φₐ} f(Λ).
```

The null rotation ``𝐑_{\boldsymbol{ℓ ξ}}`` does not satisfy such a
simple equivariance condition.  In fact no equivariance weight for
null rotations can exist at all;[^5] it must be accounted for by
mixing of the corresponding field components.  Which null rotations do
that mixing — those fixing ``\boldsymbol{ℓ}`` or those fixing ``𝐧`` —
is a sharper question than it looks; see ["Which null
direction?"](@ref which_null_direction) below.  But the effect of
``𝐑_{φₐ𝐭𝐳}`` can be accounted for with a simple multiplication by
``e^{b φₐ}`` — which we will see is essentially the conformal factor
``κ`` to the power ``b``.  So information about the direction of
propagation and tangent direction needed for a spin-weighted field is
entirely contained in the ``𝐑_K`` factor.

[^4]: The exponentials in the equivariance condition are called
    "characters" — homomorphisms from the group to the multiplicative
    group of complex numbers.  (Specifically, these are "general",
    rather than "unitary" characters.)  They can have various
    "weights"; in this case, ``-is`` with half-integer ``s``, and
    ``b`` which could be any complex number.  The spin weight is
    restricted in this case because the corresponding rotations are
    compact, while the boosts are not.  Our applications to field
    components will require that ``s`` and ``b`` are half-integer real
    numbers for spinor fields, and integer real numbers for tensor
    fields.

[^5]: The group ``N`` on its own does actually have characters
    ``e^{λ⋅ξ}`` for ``λ∈ℂ²`` (both components of which must be purely
    imaginary for *unitary* characters).  But considered as a
    *sub*group of ``MAN``, it does not have any nontrivial character.
    The easiest way to see this is to consider the requirement that a
    character must be a homomorphism.  If we have two group elements
    ``Λ₁,\, Λ₂ ∈ MAN``, the character of their group commutator must
    satisfy ``χ(Λ₁Λ₂Λ₁⁻¹Λ₂⁻¹) = χ(Λ₁)χ(Λ₂)χ(Λ₁)⁻¹χ(Λ₂)⁻¹ = 1``.  Now
    take ``Λ₁ = 𝐑_{φₐ}`` and ``Λ₂ = 𝐑_{\boldsymbol{ℓ ξ}}``.  Then
    their group commutator is ``1 + \tfrac{1}{2}(e^{φₐ} − 1)
    \boldsymbol{ℓ ξ}``, so choosing ``φₐ = \ln 2`` the group
    commutator is actually just ``𝐑_{\boldsymbol{ℓ ξ}}``.  Therefore,
    the homomorphism condition shows that ``χ(𝐑_{\boldsymbol{ℓ ξ}}) =
    1`` for all ``\boldsymbol{ξ}`` — the trivial character.  In fact,
    ``N`` is exactly the "commutator subgroup" [Aluffi_2009](@cite) of
    ``MAN``, which implies that the *only* characters of ``MAN`` are
    those of ``MA``.

That is, to transform the fields, we first need to evaluate them on
``𝐑_K``, then include the effects of the null rotation and boost.
This is the reason why the [`aberration`](@ref Scri.aberration)
function just computes the ``K`` factor, rather than the full ``KAN``
decomposition produced by [`Quaternionic.KAN`](@extref Quaternionic
:jl:function:`Quaternionic.KAN`).  The spin rotation is included in
the ``K`` factor, and is dealt with by treating the spin-weighted
spherical harmonics as functions on ``\mathrm{Spin}(3)``
[Boyle_2016](@cite); the boost weight is included separately simply by
multiplication; and the null rotation is included by [mixing
components](@ref Mixing:-the-peeling-tower).

## Flagpoles and flagplanes

It might help to step back for a moment and remember the broader goal.
We are dealing with polarized fields propagating radially along null
directions, and finding how those fields transform under Lorentz
transformations.  So in each direction, and for each observer, we need
to be able to describe the null direction and the polarization.  The
literature has two geometric objects that are used to describe these
two pieces of information: the "flagpole" and the "flagplane".

We first define the null vector

```math
𝐊
≔ Λ \boldsymbol{ℓ} \bar{Λ}
∝ 𝐑_K \boldsymbol{ℓ} 𝐑_K
= 𝐑_{k̂} \boldsymbol{ℓ} 𝐑_{k̂}.
```

This is the direction of propagation of the field, as a full
four-vector of the spacetime.  It's important to note that ``𝐑_{k̂}``
picks out the direction, and that *direction* is invariant under
``𝐑_{γ𝐲'𝐱'}``, ``𝐑_{φₐ𝐭''𝐳''}``, and ``𝐑_{\boldsymbol{ℓ}'''
\boldsymbol{ξ}'''}``.  This is precisely the object often referred to
as the "flagpole" in the spinor literature [PenroseRindler_1984; p.
37](@cite).

Now, to specify the polarization of the field components, we need to
specify a fiducial direction in the plane orthogonal to ``𝐊``
[Boyle_2016](@cite).  Conventionally, when dealing only with
rotations, this fiducial direction is just the ``𝐱`` basis vector
after rotation.  We can define[^6]

```math
𝐋 ≔ Λ 𝐱 \bar{Λ}.
```

We would *like* to use this as the fiducial direction for the
polarization, but as long as the ``N`` factor is non-trivial, it is no
longer be spatial.  It is not immediately clear which value of the
field should be extracted from the original frame to allow us to
compute the transformed field.  Somehow, we have to get "as close as
we can" to some sort of concordance between the frames.

[^6]: Note that we often include ``𝐲`` by constructing the "complex
    vector" ``𝐦 = 𝐑_K \tfrac{𝐱 + i 𝐲}{\sqrt{2}} 𝐑̄_K``.  This is,
    in fact, how spin-weighted spherical harmonics are defined
    [Boyle_2016](@cite).

Writing out the components ``\boldsymbol{ξ} = ξˣ𝐱 + ξʸ𝐲``, a
straightforward calculation shows that the final transformed direction
of ``𝐱`` is

```math
x''''
= Λ 𝐱 \bar{Λ}
=  𝐑_K\, \left(𝐱 + ξˣe^{φₐ}\boldsymbol{ℓ}\right)\, 𝐑̄_K
= 𝐋 + ξˣe^{φₐ}𝐊.
```

Now, the set

```math
\{c𝐊 + 𝐑_K\, 𝐱\, 𝐑̄_K \mathrel{|} c∈ℝ\}
```

is the set of all possible transformations of ``𝐱`` for the given
flagpole.  So we might take this equivalence class as the object that
determines the fiducial direction for the polarization; any given
observer will just need to choose the particular member that is
spatial relative to their own frame.

It might help to make contact with the rest of the literature here.
Conventionally, we take all *positive* scalings of that set to be the
"flagplane":

```math
\{a𝐊 + b𝐋 \mathrel{|} a,b∈ℝ; b>0\}.
```

[PenroseRindler_1984](@Citet), for example, take great pains to
explain that a flagplane is a *half*-plane — the purpose being to
ensure that the flagplane must rotate through a full ``2π`` to return
to itself, rather than just the ``π`` that would return a full plane
to itself.  This is exactly the geometry encoded by an *oriented*
plane, which is what a simple bivector represents.[^7]  Therefore, we
can also represent the flagplane by the bivector

```math
𝐅 = 𝐊 ∧ 𝐋 = Λ \left(\boldsymbol{ℓ} ∧ 𝐱\right) \bar{Λ}.
```

This bivector scales as ``e^{φₐ}`` under the ``A`` factor, but is
invariant under the ``N`` factor.  This could be a more elegant
representation of the flagplane than the conventional one, because it
is a geometric object that also participates in the algebra of the
spacetime, and retains information about the magnitude.

[^7]: A "simple" bivector is one that can be written as the wedge
    product of two vectors — a distinction that matters in four or
    more dimensions.  For example, the bivector ``𝐭∧𝐱`` is simple, but
    ``𝐭∧𝐱+𝐲∧𝐳`` is not.  Only a *simple* bivector represents an
    oriented plane.  This should be familiar from any of the standard
    treatments of differential geometry [MisnerThorneWheeler_1973,
    DoranLasenby_2003, Frankel_2011, Lee_2019](@cite).

!!! details "Comparison to Penrose-Rindler"

    [PenroseRindler_1984](@Citet) use the two-spinor formalism
    extensively.  This is closely related to the quaternion formalism:
    if ``Q = a𝟏 + b𝐢 + c𝐣 + d𝐤`` is the quaternion, the corresponding
    two-spinor is ``(ξ, η) = (a + id, c + ib)``.  Alternatively, given
    a two-spinor ``(ξ, η)``, we re-imagine the unit imaginary ``i`` as
    the bivector ``𝐤`` and we have ``Q = ξ + 𝐣 η``.  We can then
    check, for example, that the ``𝐋`` and ``𝐊`` vectors defined above
    are precisely equal to the ones given by Eqs. (1.4.14) and
    (1.4.16), respectively, of [PenroseRindler_1984](@cite).

Our field components are really functions of both the flagpole picking
out the direction of propagation, and the flagplane picking out the
tangent orientation that allows us to describe the polarization.  In
the language of the equivariance condition above, the flagpole is what
``𝐑_{φₐ𝐭𝐳}`` rescales (boost weight) and the flagplane is what
``𝐑_{γ𝐲𝐱}`` rotates (spin weight).  The ``AN`` factors preserve
both the flagpole's direction and the flagplane, so all the
*positional* information — which flagpole, and which flagplane through
it — is encoded in the ``K`` factor alone.

## Working in the unprimed frame

At this point, we are in a good position to understand a subtlety in
the implementation of the transformation.  The field is expressed in
either frame on a series of time slices, and on each time slice the
field components are decomposed in spin-weighted spherical harmonics.
During the transformation, these have to be evaluated at a series of
"pixels" — rotors from ``\mathrm{Spin}(3)`` [Boyle_2016](@cite), which
correspond as discussed above to flagpoles and flagplanes.
Specifically, there is a natural choice of pixels in the *transformed*
frame which will allow us to compute the mode weights.  We will denote
those pixels as ``𝐑'_p``.  But the field at each such pixel has to be
computed from the original data.  So we need to evaluate it at some
corresponding ``𝐑_p``.  The previous section tells us that the
flagpole and flagplane should be the same for both:

```math
\begin{aligned}
𝐑'ₚ \boldsymbol{ℓ}' 𝐑̄'ₚ
&∝ 𝐑ₚ \boldsymbol{ℓ} 𝐑̄ₚ,
\\
𝐑'ₚ \left(𝐱' ∧ \boldsymbol{ℓ}'\right) 𝐑̄'ₚ
&∝ 𝐑ₚ \left(𝐱 ∧ \boldsymbol{ℓ}\right) 𝐑̄ₚ.
\end{aligned}
```

The reference vectors being acted on here are related by the
transformation:

```math
\boldsymbol{ℓ}' = Λ \boldsymbol{ℓ} \bar{Λ},
\qquad
𝐱' = Λ 𝐱 \bar{Λ}.
```

Therefore, the expression above can be rewritten as

```math
\begin{aligned}
𝐑'ₚ Λ \boldsymbol{ℓ} \bar{Λ} 𝐑̄'ₚ
&∝ 𝐑ₚ \boldsymbol{ℓ} 𝐑̄ₚ,
\\
𝐑'ₚ Λ \left(𝐱 ∧ \boldsymbol{ℓ}\right) \bar{Λ} 𝐑̄'ₚ
&∝ 𝐑ₚ \left(𝐱 ∧ \boldsymbol{ℓ}\right) 𝐑̄ₚ,
\end{aligned}
```

so we might expect that the correct ``𝐑ₚ`` to use would be

```math
𝐑ₚ = 𝐑'ₚ Λ.
```

Or better yet, the previous section argued that all we actually need
is the ``K`` factor of that product, so we might expect that the right
choice would be

```math
𝐑ₚ = K\left(𝐑'ₚ Λ\right),
```

where ``K`` is the function that extracts the ``K`` factor from a
``KAN`` decomposition.

*However*, the subtlety is that the `𝐑'ₚ` actually generated by the
code is just an array of four numbers.  And throughout the code, the
frame of those four numbers is never explicitly specified.  Therefore,
such an array is inherently always in the same frame.  Since the rest
of the calculation assumes the frame is the unprimed one, the code is
effectively treating the array as being in the *unprimed* frame.  On
the other hand, the code generating that array doesn't know anything
about that distinction, and our information is relative to the
*primed* frame.  That is, in principal we expect the code to generate

```math
𝐑'ₚ = R'ₚʷ + R'ₚˣ 𝐢' + R'ₚʸ 𝐣' + R'ₚᶻ 𝐤',
```

but what it actually stores is the array of components

```julia
𝐑'ₚ = [R'ₚʷ,  R'ₚˣ,  R'ₚʸ,  R'ₚᶻ]
```

(Note that the previous equation was in "math" face, while this one is
in "code" face, to distinguish what *should* happen analytically from
what *does* happen practically.)  To resolve the mismatch, we just
have to conjugate by ``Λ`` to get the correct rotor relative to the
primed frame, which we can do implicitly and arrive at

```julia
𝐑ₚ = K(Λ * 𝐑'ₚ)
```

as the correct rotor to use in the unprimed frame.  This is exactly
what [`aberration`](@ref Scri.aberration)`(R′ₚ, Λ)` computes.  Again,
the ``A`` and ``N`` factors that the extraction discards are exactly
the ones accounted for elsewhere: the conformal factor ``κ``, and the
component mixing, which re-enters in as a null rotation about the
generator ``𝐧``, as discussed next.

## [Which null direction?](@id which_null_direction)

!!! note "TO DO"

    I need to rewrite this section entirely.  I should use
    ``ΘN_{\boldsymbol{ℓ}} = N_n``.

    I guess the point is that the same Λ can be decomposed in two
    different ways, but they still mean the same thing.  So for
    understanding the coordinates on the null cone, we decompose with
    respect to the null ray corresponding to the given point on the
    sphere.  But to understand the transformation of the field, we
    decompose with respect to the generator of the null hypersurface
    carrying the field.  Different observers generate different
    tetrads, but they are both (pseudo-)orthonormal tetrads, so they
    are related by a Lorentz transformation.  Moreover, the two
    tetrads are both adapted to the same null hypersurface, so they
    are related by a transformation in the stabilizer of the
    generator.

Everything above treated ``\boldsymbol{ℓ}`` as the special null
direction: the rotation of ``M`` was chosen to fix ``\boldsymbol{ℓ}``,
the boost of ``A`` was chosen to be along that direction, and the null
rotation of ``N`` was chosen to fix that vector.  But nothing in the
construction requires this choice of null direction.  For example,
fixing the conjugate null vector ``𝐧`` works just as well, and gives
a second Iwasawa decomposition:

```math
\mathrm{Spin}^+(3,1) = K\,A\,N_{\boldsymbol{ℓ}} = K\,A\,N_{𝐧},
```

with the *same* ``K`` and ``A``, and the nilpotent subgroups
subscripted by the null vector they fix.  The two decompositions of a
single rotor have different factors, so we must face a question the
group theory alone cannot answer: when transforming radiation data,
which null direction is the right one to build the decomposition
around?

Start with how much the two choices *share*.  Write

```math
M = \left\{ \exp\left[\frac{γ}{2} 𝐲𝐱\right] \right\}
\mathrel{\cong_{\text{Grp}}}
\mathrm{Spin}(2)
```

for the rotations about ``𝐳`` — the Hopf fiber inside ``K`` — so that
``MA`` is the combined *spin–boost* subgroup.  Every element of ``MA``
preserves *both* null rays at once, merely rescaling them oppositely:

```math
𝐑\, \boldsymbol{ℓ}\, 𝐑̃ = e^{φₐ}\, \boldsymbol{ℓ},
\qquad
𝐑\, 𝐧\, 𝐑̃ = e^{-φₐ}\, 𝐧,
\qquad
𝐑 = 𝐑_{γ𝐳}\, 𝐑_{φₐ𝐳} ∈ MA.
```

The full stabilizer of the ray ``[\boldsymbol{ℓ}]`` is the
four-parameter subgroup ``MAN_{\boldsymbol{ℓ}}``, and the stabilizer
of ``[𝐧]`` is ``MAN_{𝐧}``.  (Such stabilizers of null rays are
called *parabolic subgroups*, and the piece ``MA`` common to both is
called their *Levi factor*.)  The two stabilizers differ *only* in
their nilpotent factors.  Now recall [the equivariance condition](@ref
"Reversing the order"): spin weight and boost weight are the
multipliers ``e^{isγ}`` and ``e^{bφₐ}`` picked up under
right-multiplication by ``MA`` — they are characters of ``MA`` alone.
Since ``MA`` is blind to the difference between ``\boldsymbol{ℓ}`` and
``𝐧``, so is the entire weight apparatus.  Spin weight, boost weight,
and the conformal rescaling are perfectly symmetric between the two
null directions, and cannot distinguish them even in principle.

Could a *null-rotation weight* break the tie?  No — and the reason is
instructive.  Conjugating a null rotation by the boost dilates its
generator,

```math
𝐑_{φₐ𝐳}\, \boldsymbol{ℓξ}\, 𝐑̃_{φₐ𝐳} = e^{φₐ}\, \boldsymbol{ℓξ},
```

since the boost rescales ``\boldsymbol{ℓ}`` and leaves the transverse
``\boldsymbol{ξ}`` alone.  A multiplier ``χ`` for null rotations,
consistent with the group structure, would therefore have to satisfy
``χ(𝐑_{\boldsymbol{ξ}}) = χ(𝐑_{e^{φₐ}\boldsymbol{ξ}})`` for every
``φₐ``; taking ``φₐ → -∞`` and using continuity forces ``χ ≡ 1``.
There is no such thing as null-rotation weight.  A field component can
at best be *invariant* under a family of null rotations; otherwise it
must *mix* with other components — which is exactly what happens, on
[the tetrad page](@ref "BMS Action on the Tetrad").

So the choice cannot come from the group, nor from the weights.  It
comes from the *surface carrying the data*.  A null hypersurface
singles out, at each of its points, exactly one null direction — its
generator, the degenerate direction of its induced metric — and a
transformation relating two frames adapted to the surface must
preserve that direction.  The frame transitions therefore live in the
generator's stabilizer.  For radiation at ``ℐ⁺`` the generator is
``𝐧``, not ``\boldsymbol{ℓ}`` — even though the radiation propagates
along ``\boldsymbol{ℓ}`` — because ``ℐ⁺`` is ruled like an
*absorption* cone rather than an emission cone; this is explained at
[From the null cone to ``ℐ⁺``](@ref from_cone_to_scri).  The
transitions at ``ℐ⁺`` thus lie in ``MAN_{𝐧}``, and components mix by
null rotations about the generator.  At ``ℐ⁻``, and on the outgoing
null cone of an emitter, the generator is ``\boldsymbol{ℓ}`` and the
roles revert to ``MAN_{\boldsymbol{ℓ}}``.

None of this demotes the ``\boldsymbol{ℓ}``-adapted decomposition used
above; it just delimits its two jobs.  Its ``K`` factor is the *base
map* — which generator, which flagpole and flagplane — and because
``𝐑_K`` is a *rotation*, it carries the whole null pair
``(\boldsymbol{ℓ}, 𝐧)`` at the pole onto the pair at the new point
simultaneously, so the same ``K`` serves no matter which stabilizer
the data surface selects.  Its ``A`` factor carries the boost weight,
and ``A`` is common to both parabolics anyway.  The one factor that is
*not* shared, ``N_{\boldsymbol{ℓ}}``, is also the one factor that
never acts on the data by a weight; its physical effect re-enters in
mirrored form, as the null rotation about the generator ``𝐧`` whose
parameter is computed directly from the coordinate transformation on
[the tetrad page](@ref "BMS Action on the Tetrad").

## Reinterpreting ``i``

TL;DR: ``i ∈ ℂ`` is replaced by ``𝐈₃ = 𝐈𝐭``, the spatial
pseudoscalar.  It actually transforms whenever ``𝐭`` transforms, but
if we just write expressions in terms of ``𝐈₃' = 𝐈𝐭'`` without
explicitly transforming ``𝐭``, that should be fine, because the
interpretation of ``i`` also needs to change.  And at that point, it's
just a bookkeeping device, so we don't need to worry about the fact
that it transforms.  The reason ``𝐈₃`` appears is because it is
central in the spatial subalgebra; it commutes with everything *in
that subalgebra*, which is why it can act like ``i``.

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
