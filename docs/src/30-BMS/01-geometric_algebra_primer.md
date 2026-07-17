# A Primer on Geometric Algebra

!!! tip "The short version"
    Geometric Algebra is — among many other things — a very simple and
    intuitive yet powerful way to implement rotations and boosts.  The
    core idea is to multiply vectors together.  This "geometric product"
    is not generally commutative, but it has several key properties:

      1. it is associative,
      2. it is distributive,
      3. it commutes with scalar multiplication, and
      4. the product of a vector with itself is just the inner (dot)
         product of that vector with itself.

    The last property is what makes the geometric product so interesting,
    because it means that parallel vectors commute, while orthogonal
    vectors anticommute.  This makes it particularly easy to work with
    orthogonal basis elements.  More generally, because of distributivity,
    we can separate a vector into parts that are parallel and orthogonal
    to another vector:
    ```math
    𝐚𝐛 = 𝐚\left(𝐛_∥ + 𝐛_⟂\right) = \left(𝐛_∥ - 𝐛_⟂\right)𝐚.
    ```
    This last result begins to show that reflections can be implemented
    by multiplication, and reflections give rise to rotations and boosts.
    (See the discussion of the Cartan-Dieudonné theorem for details.)

    While the idea of a rotation about an axis doesn't generalize, the
    idea of "rotation in a plane" works in any number of dimensions,
    and the generator of a rotation in a plane is simply the product of
    two orthogonal vectors in that plane.  Exponentiating half the
    "angle" times the generator produces a *rotor*,
    ```math
    R = \exp \left[\frac{θ}{2}𝐱𝐲\right]
    = \cos \frac{θ}{2} + 𝐱𝐲 \sin \frac{θ}{2},
    ```
    which transforms any vector by conjugation, ``𝐯' = R𝐯R̃`` (the
    tilde is the *reverse*, defined below), so that composition of
    transformations is just multiplication of rotors.  The character of
    the transformation is fixed by the square of the generator:
    ``(𝐱𝐲)² = -1`` produces the trigonometric form above — an ordinary
    rotation; ``(𝐭𝐳)² = +1`` makes the functions hyperbolic — a boost;
    and a *null* generator truncates the series — a "null rotation".
    All three cases are worked out on [The Lorentz Group](@ref) page.

Geometric Algebra (GA) is a powerful mathematical framework that
unifies various algebraic systems, including complex numbers,
quaternions, vector calculus, and differential forms.  It provides a
geometric interpretation of algebraic operations, making it
particularly useful in physics, engineering, and computer graphics.
Nonetheless, it is simple enough to be accessible to anyone familiar
with basic algebra.  GA is fundamentally identical to Clifford
Algebra; the different names reflect different emphases and traditions
in the literature.  GA is typically developed over ``ℝ``, rather than
``ℂ`` — the complex structures we usually encounter in physics
appearing naturally within the algebra, rather than being introduced
_ad hoc_.  GA also emphasizes the geometric interpretation of the
algebraic structures, over what is usually an algebraic emphasis in
literature using the name Clifford Algebra.

Before getting into the details of null tetrads and Lorentz
transformations, we need to review some of the basics of Geometric
Algebra.  This is not meant to be a comprehensive introduction to GA,
but rather a quick primer on the key concepts and operations that we
will need for our purposes.  For a more comprehensive introduction,
see [DoranLasenby_2003](@cite).

## The geometric product

Geometric Algebra starts with a _real_ vector space ``𝕍``, equipped
with an inner product taking a pair of vectors ``𝐯, 𝐰 ∈ 𝕍`` to
``𝐯⋅𝐰 = 𝐰⋅𝐯 ∈ ℝ``.  We then introduce a product called the
_geometric product_, which is associative and distributive, but not
necessarily commutative.  We express the geometric product between
vectors ``𝐯`` and ``𝐰`` simply as juxtaposition: ``𝐯𝐰``.  The
geometric product is essentially the tensor product, subject to the
identification that ``𝐯𝐯 = 𝐯⋅𝐯``.  We also have compatibility with
scalar multiplication, so that we have

```math
s(𝐯𝐰) = (s𝐯)𝐰 = 𝐯(s𝐰) = (𝐯𝐰)s
\qquad \text{for any }s ∈ ℝ.
```

These rules are enough to extend the algebra to arbitrary dimensions,
using inner products of arbitrary (even degenerate) signature.

An example is very helpful in clarifying.  Consider the vector space
``ℝ²``, with the standard basis vectors ``𝐱`` and ``𝐲``.  Consider
the sum ``𝐱+𝐲``.  The product of this vector with itself is
identified with the inner product

```math
(𝐱+𝐲)(𝐱+𝐲) = (𝐱+𝐲)⋅(𝐱+𝐲) = 2.
```

On the other hand, we can use the distributive property to derive

```math
\begin{aligned}
(𝐱+𝐲)(𝐱+𝐲)
&= 𝐱𝐱+𝐱𝐲+𝐲𝐱+𝐲𝐲 \\
&= (𝐱𝐱+𝐲𝐲) + (𝐱𝐲+𝐲𝐱) \\
&= (𝐱⋅𝐱+𝐲⋅𝐲) + (𝐱𝐲+𝐲𝐱) \\
&= 2 + (𝐱𝐲+𝐲𝐱).
\end{aligned}
```

Comparing the two expressions, we see that the last term,
``𝐱𝐲+𝐲𝐱`` must vanish:

```math
𝐱𝐲 = -𝐲𝐱.
```

That is, these orthogonal vectors anticommute under the geometric
product.  This has an important consequence:

```math
(𝐱𝐲)(𝐱𝐲) = 𝐱(𝐲𝐱)𝐲 = -𝐱(𝐱𝐲)𝐲 = -(𝐱𝐱)(𝐲𝐲) = -1.
```

That is, ``(𝐱𝐲)² = -1``; the product ``𝐱𝐲`` _is the unit
imaginary_ associated to the ``𝐱``-``𝐲`` plane.

Obviously, parallel vectors commute, since they can be expressed as
scalar multiples of each other and scalars commute with all vectors.
These are the two critical features of the geometric product: parallel
vectors commute, while orthogonal vectors anticommute.  Combined with
associativity and distributivity, these properties allow us to
calculate quite general geometric products in arbitrary dimensions.

Specifically, we can decompose the geometric product of two vectors
into symmetric and antisymmetric parts:[^1]

```math
\begin{aligned}
𝐯𝐰 &= \frac{1}{2}(𝐯𝐰 + 𝐰𝐯) + \frac{1}{2}(𝐯𝐰 - 𝐰𝐯) \\
&= (𝐯⋅𝐰) + (𝐯 ∧ 𝐰),
\end{aligned}
```

where ``𝐯 ∧ 𝐰`` is called the exterior or [wedge
product](https://en.wikipedia.org/wiki/Wedge_product), producing a
[bivector](https://en.wikipedia.org/wiki/Bivector).  When ``𝐯²≥0``
and ``𝐰²≥0``, this result is a _general complex number_ associated
with the plane spanned by ``𝐯`` and ``𝐰``, with ``𝐯⋅𝐰`` being the
real part and ``𝐯 ∧ 𝐰`` being the imaginary part which squares to a
negative number.  Note that the wedge product corresponds to the usual
cross product in three dimensions, but generalizes to arbitrary
dimensions and signatures.

[^1]: It is remarkable that this formula actually has a scalar being
    added to the wedge product of two vectors — which is a rank-2
    tensor.  In Physics, we are frequently taught that scalars must
    never be added to vectors — never mind tensors!  This is typically
    good for helping students catch elementary mistakes, but not
    actually necessary.  Mathematicians routinely define the tensor
    space to allow for adding arbitrary ranks together.

!!! important "Geometric Algebra generates Complex Algebra"

    The geometric product of two vectors is precisely a complex number,
    with the real part being the inner product of the vectors, and the
    imaginary part being a bivector representing the plane spanned by
    those vectors.  This is a generalization of the fact that the
    product of two orthogonal vectors is a bivector that squares to -1,
    and thus can be identified with the unit imaginary.

## Reflections and rotations

One of the reasons Geometric Algebra is so powerful is that it
provides a very natural way to represent reflections — which, in turn,
give rise to orthogonal and conformal transformations in all
dimensions and (nondegenerate) signatures.

First, note the fact that the geometric product of two vectors tracks
both the "dot product" and "cross product" means that we can often
find _inverses_ of vectors, which simplifies many calculations.
Specifically, if the norm of a vector ``𝐧`` is nonzero, then we can
define its inverse as

```math
𝐧⁻¹ = \frac{𝐧}{𝐧²},
```

where the denominator is just a scalar.  Obviously, we then have
``𝐧 𝐧⁻¹ = 𝐧 𝐧 / 𝐧² = 1``.

Now, choose any invertible vector ``𝐧``.  Any other vector ``𝐯``
decomposes into a part that commutes with ``𝐧`` and a part that
anticommutes with ``𝐧`` — which we denote as ``𝐯_∥`` and ``𝐯_⟂``,
respectively.  Given the properties of the geometric product shown
above, we have

```math
-𝐧 𝐯 𝐧⁻¹ = -𝐧 𝐯_∥ 𝐧⁻¹ - 𝐧 𝐯_⟂ 𝐧⁻¹ = -𝐧𝐧⁻¹ 𝐯_∥ + 𝐧 𝐧⁻¹ 𝐯_⟂
= -𝐯_∥ + 𝐯_⟂.
```

That is, this negative conjugation by ``𝐧`` reflects the vector
``𝐯`` along the line defined by ``𝐧``; reflections are represented
as simple conjugations in the algebra.  We can compose reflections,
just by applying this transformation repeatedly, which is equivalent
to negative conjugation by the product of the vectors defining the
reflections.  For any such ``𝐧``, the negative conjugation results in
a reflection, and any reflection can be represented in this way for
some choice of ``𝐧``.  But note that the choice of ``𝐧`` is not
unique; ``-𝐧`` will achieve exactly the same reflection.  Thus, the
multiplicative group of unit vectors is a double cover of the group of
reflections.

This may seem like a trivial curiosity, but it has profound
implications because of the Cartan-Dieudonné theorem
[Garling_2011](@cite):

> If ``T`` is an isometry of a regular quadratic space ``(E, Q)``,
> ``T`` is the product of at most ``\mathrm{dim}\, E`` simple
> reflections.

For our purposes, the isometries are just the orthogonal group, and a
"regular quadratic space" is just a vector space with a
_nondegenerate_ inner product (corresponding to the quadratic form
``Q`` above).  We are only interested in _real_ vector spaces, so our
vector space is essentially just ``ℝ^{p,q}``, where ``p`` is the
number of positive terms in the signature, and ``q`` is the number of
negative terms.  (E.g., we will take Minkowski space as ``ℝ^{3,1}``.)
Thus, we might rephrase the theorem more simply as

> Any orthogonal transformation of ``ℝ^{p,q}`` can be expressed as the
> product of at most ``p+q`` simple reflections.

There is a corollary that is also important for our purposes:

> Any _special_ orthogonal transformation of ``ℝ^{p,q}`` can be
> expressed as the product of _an even number_ of at most ``p+q``
> simple reflections.

Composing two reflections is the same as conjugation by a product of two unit
vectors — a _rotor_ ``R`` — acting as ``𝐯 ↦ R𝐯R̃``.  The tilde here
denotes the _reverse_ operation, which swaps the order of the vectors
in any product; for a product of unit vectors this is the inverse, so
``RR̃ = 1``.  For a rotor ``R = \exp[θ𝐱𝐲/2]`` generated by the plane
``𝐱``-``𝐲``, any vector orthogonal to that plane anticommutes with
both ``𝐱`` and ``𝐲``, hence _commutes_ with ``𝐱𝐲``, and is left
unchanged by the conjugation; the components in the plane are rotated
by the full angle ``θ`` — the same halving of the exponent familiar
from quaternions and spinors.

In fact, the even products of unit vectors form the _Spin_ group, and
it is a _double_ cover of the special orthogonal group: just as with a
single reflection, ``±R`` give the same transformation ``𝐯 ↦ R𝐯R̃``,
but now the two sheets are connected, so a continuous rotation by
``2π`` returns the transformation to the identity while returning ``R``
to ``-R``.  This is the key to understanding how spinors arise in
physics.  A spinor is an object the Spin group acts on _directly_ — by a
single rotor ``R``, rather than by the two-sided conjugation ``R\,(⋅)\,R̃``
— and so it detects the very sign that the orthogonal transformation
forgets.  The rotors we build in the following pages are exactly such
objects, and the spin-weighted fields we ultimately transform are built
from spinors — the "square roots" of the null tetrad legs — which is why
the Spin group, not merely the rotation group, is the right language
throughout.

## Higher-dimensional products

The bivector is part of a broader pattern.  The product of two
orthogonal vectors — equivalently, the antisymmetric part of the
geometric product of any two vectors — is a bivector representing the
plane they span, and it records the
[attitude](https://en.wikipedia.org/wiki/Orientation_(geometry)) of
that plane, the
[orientation](https://en.wikipedia.org/wiki/Orientability) of the
bivector, and a magnitude.  The same construction climbs through the
dimensions: the product of three orthogonal vectors is a trivector
representing the volume they span, and so on.  In particular, the
product of ``d`` orthogonal vectors in a ``d``-dimensional space is a
_pseudoscalar_ representing the oriented volume of the whole space, and
it squares to either ``+1`` or ``-1`` depending on the signature.  The
pseudoscalar is usually denoted ``𝐈``, and it plays a central role in
the algebra.  It is essentially the [volume
form](https://en.wikipedia.org/wiki/Volume_form) of the space, and
provides the [Hodge
dual](https://en.wikipedia.org/wiki/Hodge_star_operator) by simple
multiplication; in each of the cases we use below — ``𝐈₂ = 𝐱𝐲``,
``𝐈₃ = 𝐱𝐲𝐳``, and ``𝐈₄ = 𝐭𝐱𝐲𝐳`` — it happens to square to ``-1``,
making it a complex structure for its space.

These graded products have a simple accounting.  There is — by
definition — exactly one scalar, and linear dependence forces the wedge
product of more than ``d`` vectors to vanish.  In between, the space of
products of ``k`` independent vectors has dimension ``\binom{d}{k}``, so
the whole algebra has dimension

```math
\sum_{k=0}^d \binom{d}{k} = 2^d.
```

There is likewise just one pseudoscalar, the product of all ``d`` basis
vectors.  In ``d = 2`` there is a single bivector, which _is_ the
pseudoscalar — the unit imaginary we met above.  That complex numbers
are then a linear combination of the scalar and that one bivector, and
so two-dimensional, is a pure coincidence of ``d = 2``.  Hamilton was
misled into expecting the same structure at ``d = 3``, where the
coincidence is instead that the number of vectors ``\binom{3}{1} = 3``
equals the number of bivectors ``\binom{3}{2} = 3``.  Chasing these
seemingly magical coincidences caused decades of confusion that only the
development of Geometric Algebra has fully resolved.

## Structure of the algebra

Collecting the grades, the algebra over a ``d``-dimensional space is a
direct sum of its grade-``k`` subspaces, each of dimension
``\binom{d}{k}`` — one row of Pascal's triangle.  The three cases we use
are ``d = 2``, ``d = 3``, and ``d = 4`` (Minkowski):

| algebra          | scalar | vector | bivector | trivector | pseudoscalar | total       |
|:-----------------|:------:|:------:|:--------:|:---------:|:------------:|:-----------:|
| ``𝒢(ℝ²)``        |   1    |   2    |    1     |           |              | ``2² = 4``  |
| ``𝒢(ℝ³)``        |   1    |   3    |    3     |     1     |              | ``2³ = 8``  |
| ``𝒢(ℝ^{3,1})``   |   1    |   4    |    6     |     4     |      1       | ``2⁴ = 16`` |

The grades we lean on are the vectors (grade 1) and the bivectors (grade
2), and in each case it is the _even_ grades — scalar, bivector, and (in
four dimensions) pseudoscalar — that form the rotors.

In ``𝒢(ℝ²)`` the single bivector ``𝐱𝐲`` is the pseudoscalar and the
unit imaginary, so the even part — scalar plus bivector — is a copy of
``ℂ``.  In ``𝒢(ℝ³)`` the three bivectors ``𝐢 = 𝐳𝐲``, ``𝐣 = 𝐱𝐳``,
``𝐤 = 𝐲𝐱`` each square to ``-1`` and multiply like Hamilton's
quaternions, so the even part — scalar plus those three bivectors — is a
copy of ``ℍ``, whose unit-norm elements are ``\mathrm{Spin}(3)``, the
rotors of spatial rotation.  In ``𝒢(ℝ^{3,1})`` the six bivectors split
into three _spatial_ ones (``𝐲𝐳, 𝐳𝐱, 𝐱𝐲``, squaring to ``-1``,
generating rotations) and three _timelike_ ones (``𝐭𝐱, 𝐭𝐲, 𝐭𝐳``,
squaring to ``+1``, generating boosts).  The eight-dimensional even part
— scalar, six bivectors, and pseudoscalar — contains
``\mathrm{Spin}⁺(3,1)`` as its unit-norm elements: the rotors we use
throughout to represent Lorentz transformations.  Finally, the spatial
pseudoscalar ``𝐈₃ = 𝐱𝐲𝐳`` and the spacetime pseudoscalar ``𝐈₄ =
𝐭𝐱𝐲𝐳`` are the objects that will [stand in for the unit
imaginary](@ref "Reinterpreting ``i``") in the null tetrad.
