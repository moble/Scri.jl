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
\boldsymbol{\ell} &= \frac{𝐭+𝐳}{\sqrt{2}}, \\
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
\boldsymbol{\ell}² = 𝐦² = 𝐦̄² = 𝐧² = 0,
```

and the only nonzero inner products are

```math
\boldsymbol{\ell} \cdot 𝐧 = -1, \qquad 𝐦 \cdot 𝐦̄ = 1.
```

This simple tetrad is aligned with the axes of the standard basis, but
we can apply any Lorentz transformation to it to get a more general
null tetrad.  In particular, we can apply an ordinary rotation to
"point" ``\boldsymbol{\ell}`` in any direction we like.  This suggests
the standard factorization of the Lorentz group, beginning with
rotations tied to the spherical coordinates, followed by rotations
about the radial direction, followed by boosts in the radial
direction, and finally followed by null rotations about the rotated
``\boldsymbol{\ell}``.  We will discuss these transformations in more
detail below, but the important point is that we decompose a general
Lorentz transformation *specifically with respect to* either
``\boldsymbol{\ell}`` or ``𝐧``.

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
\boldsymbol{\ell} = \frac{𝐭+𝐳}{\sqrt{2}}
```

as the invariant null vector.  Its complementary null vector is

```math
𝐧 = \frac{𝐭-𝐳}{\sqrt{2}},
```

which will be important.  Now, for any (not necessarily unit) vector
``\boldsymbol{ξ}`` in the ``𝐱``-``𝐲`` plane, the bivector
``\boldsymbol{\ell ξ} = -\boldsymbol{ξ \ell}`` generates a null
rotation.  More specifically, this bivector generates a boost in the
``\boldsymbol{ξ}`` direction, and *simultaneously* a rotation in the
``\boldsymbol{ξ}``-``𝐳`` plane.  Define the spinor

```math
𝐑 = \exp\left[ \frac{1}{2} \boldsymbol{\ell ξ} \right].
```

Because ``\boldsymbol{\ell}² = 0``, the exponential series terminates
after the second term, and we have

```math
𝐑 = 1 + \frac{1}{2} \boldsymbol{\ell ξ},
```

which makes calculations particularly simple.  The results on the
basis ``(𝐭, 𝐱, 𝐲, 𝐳)`` are not enlightening, but the results on
the ``(\boldsymbol{\ell}, 𝐧, 𝐦, 𝐦̄)`` basis are very interesting:

```math
\begin{aligned}
𝐑 \boldsymbol{\ell} \bar{𝐑} &= \boldsymbol{\ell}, \\
𝐑 𝐦 \bar{𝐑} &= 𝐦 + 𝐱\boldsymbol{ξ\ell}, \\
𝐑 𝐦̄ \bar{𝐑} &= 𝐦̄ - 𝐱\boldsymbol{ξ\ell}, \\
𝐑 𝐧 \bar{𝐑} &= 𝐧 + \boldsymbol{ξ} + \frac{1}{2} ξ² \boldsymbol{\ell}.
\end{aligned}
```

```math
\begin{aligned}
𝐑 \boldsymbol{\ell} \bar{𝐑}
&= \left(1 + \frac{1}{2} \boldsymbol{\ell ξ} \right) \boldsymbol{\ell} \left(1 - \frac{1}{2} \boldsymbol{\ell ξ} \right) \\
&= \boldsymbol{\ell} + \frac{1}{2} \boldsymbol{\ell ξ \ell} - \frac{1}{2} \boldsymbol{\ell \ell ξ} - \frac{1}{4} \boldsymbol{\ell ξ \ell \ell ξ} \\
&= \boldsymbol{\ell},
\end{aligned}
```

```math
\begin{aligned}
𝐑 𝐧 \bar{𝐑}
&= \left(1 + \frac{1}{2} \boldsymbol{\ell ξ} \right) 𝐧 \left(1 - \frac{1}{2} \boldsymbol{\ell ξ} \right) \\
&= 𝐧 + \frac{1}{2} \boldsymbol{\ell ξ} 𝐧 - \frac{1}{2} 𝐧 \boldsymbol{\ell ξ} - \frac{1}{4} \boldsymbol{\ell ξ 𝐧 \ell ξ} \\
&= 𝐧 - \frac{1}{2} \boldsymbol{(\ell 𝐧 + 𝐧 \ell) ξ} - \frac{1}{4} \boldsymbol{\ell 𝐧 \ell ξ²} \\
&= 𝐧 + \boldsymbol{ξ} + \frac{1}{2} \boldsymbol{ξ}² \boldsymbol{\ell},
\end{aligned}
```

```math
\begin{aligned}
𝐑 𝐱 \bar{𝐑}
&= \left(1 + \frac{1}{2} \boldsymbol{\ell ξ} \right) 𝐱 \left(1 - \frac{1}{2} \boldsymbol{\ell ξ} \right) \\
&= 𝐱 + \frac{1}{2} \boldsymbol{\ell ξ} 𝐱 - \frac{1}{2} 𝐱 \boldsymbol{\ell ξ} - \frac{1}{4} \boldsymbol{\ell ξ 𝐱 \ell ξ} \\
&= 𝐱 + \frac{1}{2} \boldsymbol{\ell ξ} 𝐱 + \frac{1}{2} \boldsymbol{\ell 𝐱 ξ} - \frac{1}{4} \boldsymbol{\ell² ξ 𝐱 ξ} \\
&= 𝐱 + \boldsymbol{\ell (ξ \cdot 𝐱)}
\end{aligned}
```

```math
\begin{aligned}
𝐑 𝐲 \bar{𝐑} &= 𝐲 + \boldsymbol{\ell (ξ \cdot 𝐲)}
\end{aligned}
```

```math
\begin{aligned}
𝐑 𝐈₃𝐲 \bar{𝐑}
&= \left(1 + \frac{1}{2} \boldsymbol{\ell ξ} \right) 𝐈₃𝐲 \left(1 - \frac{1}{2} \boldsymbol{\ell ξ} \right) \\
&= 𝐈₃𝐲 + \frac{1}{2} \boldsymbol{\ell ξ} 𝐈₃𝐲 - \frac{1}{2} 𝐈₃𝐲 \boldsymbol{\ell ξ} - \frac{1}{4} \boldsymbol{\ell ξ 𝐈₃𝐲 \ell ξ} \\
&= 𝐈₃𝐲 + \frac{1}{2} \boldsymbol{\ell ξ} 𝐈₃𝐲 - \frac{1}{2} 𝐈₃𝐲 \boldsymbol{\ell ξ} - \frac{1}{4} \boldsymbol{\ell ξ 𝐈₃𝐲 \ell ξ} \\
\end{aligned}
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
