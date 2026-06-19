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
to understand how the Lorentz group shows up in Geometric Algebra.

## Null rotations

A null rotation is a very particular type of Lorentz transformation
(sometimes called a [parabolic
transformation](https://en.wikipedia.org/wiki/Lorentz_group#Parabolic))
that leaves a chosen null vector invariant.  For simplicity, let us
choose the null vector

```math
\boldsymbol{ℓ} = \frac{𝐭+𝐳}{\sqrt{2}}
```

as the invariant null vector.  Its complementary null vector is

```math
𝐧 = \frac{𝐭-𝐳}{\sqrt{2}},
```

which will be important.  Now, for any (not necessarily unit) vector
``\boldsymbol{ξ}`` in the ``𝐱``-``𝐲`` plane, the bivector
``\boldsymbol{ℓ ξ} = -\boldsymbol{ξ ℓ}`` generates a null
rotation.  More specifically, this bivector generates a boost in the
``\boldsymbol{ξ}`` direction, and *simultaneously* a rotation in the
``\boldsymbol{ξ}``-``𝐳`` plane.  Define the spinor

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

## Iwasawa's ``KAN`` decomposition

[Knapp_1996](@cite) describes the Iwasawa decomposition of the Lorentz
group, which is a factorization of the group into three subgroups:

```math
G = KAN.
```

This decomposition is most useful when we have a preferred time axis
and preferred spatial direction.  We use those to pick out boosts
along that spatial direction, and null rotations that fix the
corresponding null vector.  Specifically, we have the following
subgroups:

- **``K``** is the maximal compact subgroup, which in the case of the
  Lorentz group is isomorphic to the rotation group
  ``\mathrm{SO}(3)``.  This subgroup consists of all rotations that
  preserve the spatial part of spacetime.
- **``A``** is the abelian subgroup of boosts in a fixed direction.
  We take this direction to be the z-axis, so that we have ``A =
  \left\{ \exp\left[\tfrac{φₐ}{2} \, 𝐭𝐳 \right] \mid φₐ ∈ ℝ
  \right\}``.
- **``N``** is the nilpotent subgroup of null rotations, which
  consists of all Lorentz transformations that can be represented as
  null rotations about a fixed null vector.  We take this vector to be
  the null vector ``\boldsymbol{ℓ}``, so that we have ``N = \left\{
  \exp\left[\tfrac{1}{2} \boldsymbol{ℓ ξ}\right] \mid \boldsymbol{ξ} =
  ξˣ𝐱 + ξʸ𝐲 \right\}``.

Note that ``φₐ`` here is *not* the rapidity of the overall boost if we
factor a transformation as a boost and a rotation.  Rather, it is the
rapidity of the particular boost picked out as ``A`` by this
construction.  ``N`` also contributes to the overall boost, so the
overall rapidity is not just ``φₐ``.

One nice feature of this decomposition is that we can compute the
factorization of a given Lorentz transformation ``Λ`` fairly simply.
We first introduce the idempotent

```math
u₊ = \frac{1}{2} (1 + 𝐭𝐳).
```

It is easy to verify that ``u₊𝐭𝐳 = 𝐭𝐳u₊ = u₊`` and hence the
idempotent identity: ``u₊² = u₊``.  Next, we write

```math
Λ = Rₖ Rₐ Rₙ,
```

where ``Rₖ ∈ K``, ``Rₐ ∈ A``, and ``Rₙ ∈ N``.  Starting with the
rightmost factor, we have

```math
Rₙ = \exp\left[\tfrac{1}{2} \boldsymbol{ℓ ξ}\right] = 1 + \frac{1}{2} \boldsymbol{ℓ ξ},
```

because ``\boldsymbol{ℓ}`` and ``\boldsymbol{ξ}`` anticommute, and
because ``\boldsymbol{ℓ}² = 0``, their product ``\boldsymbol{ℓ ξ}``
itself is also nilpotent.  (Remember, that's what the "N" stands for.)
We can simply compute the product ``Rₙ u₊`` by using the fact that
``\boldsymbol{ξ} u₊ =  u₊ \boldsymbol{ξ}`` and expanding terms in our
basis, then find

```math
Rₙ u₊ %&= \left(1 + \frac{1}{2} \boldsymbol{ℓ ξ}\right) u₊, \\
%  &= u₊ + \frac{1}{2} \boldsymbol{ℓ ξ} u₊, \\
%  &= u₊ + \frac{1}{2} \boldsymbol{ℓ} u₊ \boldsymbol{ξ}, \\
%  &= u₊ + \frac{1}{4\sqrt{2}} (𝐭+𝐳) (1+𝐭𝐳) \boldsymbol{ξ}, \\
%  &= u₊ + \frac{1}{4\sqrt{2}} (𝐭+𝐳+𝐭𝐭𝐳+𝐳𝐭𝐳) \boldsymbol{ξ}, \\
%  &= u₊ + \frac{1}{4\sqrt{2}} (𝐭+𝐳-𝐳-𝐭) \boldsymbol{ξ}, \\
= u₊.
```

Now, with

```math
Rₐ = \exp\left[\tfrac{φₐ}{2} \, 𝐭𝐳 \right]
= \cosh\left(\frac{φₐ}{2}\right) + \sinh\left(\frac{φₐ}{2}\right) 𝐭𝐳,
```

and the fact that ``𝐭𝐳u₊ = u₊``, we have

```math
\begin{aligned}
Rₐ u₊ &= \left[\cosh\left(\frac{φₐ}{2}\right) + \sinh\left(\frac{φₐ}{2}\right)\right] u₊,
 &= e^{φₐ/2} u₊.
\end{aligned}
```

Putting these together, we have

```math
Λ u₊ = e^{φₐ/2} Rₖ u₊.
```

Recall that ``e^{φₐ/2}`` is a strictly positive real number, and
``Rₖ`` is a pure rotation so it is "ℂ-real" — meaning that it is a
linear combination of only basis elements that do not have a factor of
``𝐭`` — whereas ``u₊`` is just a combination of 1 and a "ℂ-imaginary"
part.  Therefore, we can take the real part of this expression to find

```math
ℂ\Re\{Λ u₊\} = \frac{e^{φₐ/2}}{2} Rₖ.
```

This is just a strictly positive real number times ``Rₖ``, which we
know must have unit magnitude, so we can find ``Rₖ`` by normalizing
this real part:

```math
Rₖ = \mathrm{normalize}\left(ℂ\Re\{Λ u₊\}\right).
```

For the purposes of this package, this factor is actually all we need
to find.  But if desired, we can also separate the remaining factors
``Rₐ`` and ``Rₙ``.  The above also shows us that

```math
φₐ = 2 \ln \left( 2\left| ℂ\Re\{Λ u₊\} \right| \right),
```

And we can immediately plug this into the expression for ``Rₐ`` to
find that factor.  Finally, we can find the remaining factor by
rearranging the original factorization:

```math
Rₙ = Rₐ^{-1} Rₖ^{-1} Λ.
```

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
