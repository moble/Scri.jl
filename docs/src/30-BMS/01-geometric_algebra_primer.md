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
    idea of "rotation in a plane" works in any number of dimensions.
    For example, in three dimensions, a rotation "about" ``𝐳`` should
    really be thought of as a rotation in the ``𝐱``-``𝐲`` plane.  The
    generator of this rotation is simply ``𝐱𝐲`` — the pseudoscalar for
    that plane — times *half* the angle of rotation.  We simply
    exponentiate this to obtain the full rotation operator, where the
    exponential can be defined by the usual power series expansion:
    ```math
    R = \exp \left[\frac{θ}{2}𝐱𝐲\right]
    = \cos \frac{θ}{2} + 𝐱𝐲 \sin \frac{θ}{2}.
    ```
    Note the resemblance to Euler's formula for complex numbers, which
    is not a coincidence — the derivation relies only on the fact that
    ```math
    (𝐱𝐲)² = 𝐱(𝐲𝐱)𝐲 = -𝐱(𝐱𝐲)𝐲 = -(𝐱𝐱)(𝐲𝐲) = -1.
    ```
    This ``R`` rotates a vector ``𝐯`` as
    ```math
    𝐯' = R𝐯R̃,
    ```
    where the tilde denotes the *reverse* operation — which swaps the
    order of the vectors in any product.  In the case of ``R``, this
    just flips the sign of the ``𝐱𝐲`` term.  Note that ``RR̃=1``.  Any
    vector orthogonal to the ``𝐱``-``𝐲`` plane anticommutes with both
    ``𝐱`` and ``𝐲``, and therefore *commutes* with ``𝐱𝐲``, and is
    therefore unaffected by the rotation.  On the other hand, it is easy
    to verify that the components of ``𝐯`` in the ``𝐱``-``𝐲`` plane are
    rotated by the angle ``θ``.

    Importantly, this result holds in any number of dimensions, but it
    does depend on the signature of ``𝐱`` and ``𝐲``.  If we choose, for
    example, the ``𝐭`` and ``𝐳`` basis vectors of Minkowski, then
    ``(𝐭𝐳)² = 1`` so the trigonometric functions in the expression for
    ``R`` become hyperbolic:
    ```math
    R = \exp \left[\frac{φ}{2}𝐭𝐳\right]
    = \cosh \frac{φ}{2} + 𝐭𝐳 \sinh \frac{φ}{2}.
    ```
    This "rotation" in the ``𝐭``-``𝐳`` plane is just a boost in the
    ``𝐳`` direction with rapidity ``φ``.  The final case to consider is
    when the square of the generator is zero.  This happens when we
    multiply a null vector ``𝐧`` by another vector ``𝐯`` that is
    orthogonal to it, so that
    ```math
    (𝐧𝐯)² = 𝐧(𝐯𝐧)𝐯 = -𝐧(𝐧𝐯)𝐯 = -(𝐧𝐧)(𝐯𝐯) = 0.
    ```
    In this case, the exponential truncates after the linear term:
    ```math
    R = \exp \left[\frac{1}{2}𝐧𝐯\right]
    = 1 + \frac{1}{2}𝐧𝐯.
    ```
    This is called a ["null
    rotation"](https://en.wikipedia.org/wiki/Lorentz_group#Parabolic),
    because it does not affect the null vector ``𝐧``, though it does
    affect other directions.  All of these ``R`` objects are called
    *rotors*, even when they do not represent spatial rotations.
    Because they act on vectors by conjugation, we can express the
    composition of transformations simply by multiplying their rotors.

    Another very important type of element is the product of all vectors
    in an orthonormal basis, called the *pseudoscalar*.  In two
    dimensions, using the usual basis ``(𝐱, 𝐲)``, then the product is
    ``𝐈₂ = 𝐱𝐲``.  In three dimensions, we include ``𝐳`` and get
    ``𝐈₃ = 𝐱𝐲𝐳``.  In Minkowski, we include ``𝐭`` and get ``𝐈₄ = 𝐭𝐱𝐲𝐳``.
    These all happen to square to ``-1``, making them complex structures
    for their respective spaces.  And they are essentially the [volume
    forms](https://en.wikipedia.org/wiki/Volume_form) of their spaces,
    and can provide the [Hodge
    dual](https://en.wikipedia.org/wiki/Hodge_star_operator) by simple
    multiplication.

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

In fact, the group of special orthogonal transformations is a double
cover of the group of even products of unit vectors, which is called
the _Spin_ group.  This is the key to understanding how spinors arise
in physics, and how they are related to

## Higher-dimensional products

[#TODO: convert the "details" into a proper section.]

!!! details "Details on the geometric product"

    This is part of a broader pattern: the product of two orthogonal
    vectors — or just the antisymmetric part of the geometric product of
    any two vectors — is a bivector representing the plane spanned by
    those vectors, carrying information about the
    [attitude](https://en.wikipedia.org/wiki/Orientation_(geometry)) of
    that plane, [orientation](https://en.wikipedia.org/wiki/Orientability)
    of the bivector, and magnitude.  That extends to higher dimensions, as
    well, with the product of three orthogonal vectors being a trivector
    representing the volume spanned by those vectors, and so on.  In
    particular, the product of ``d`` orthogonal vectors in a
    ``d``-dimensional space is a pseudoscalar representing the oriented
    volume of the entire space, and squaring to either +1 or -1 depending
    on the signature of the inner product.  The pseudoscalar is often
    denoted as ``𝐈``, and it plays a critical role in the algebra, as we
    will see below.

    These products have interesting properties.  There is — by definition
    — just one scalar.  And linear dependence shows that the wedge product
    of more than ``d`` vectors must vanish.  In between, the space of
    products of ``k`` independent vectors has dimension ``\binom{d}{k}``,
    and the entire algebra has dimension
    ```math
    \sum_{k=0}^d \binom{d}{k} = 2^d.
    ```
    Note that there is just one pseudoscalar, which is the product of all
    the basis vectors.  In ``d=2``, there is just one bivector, which is
    the pseudoscalar, which we've seen is the complex unit imaginary.  The
    fact that complex numbers are a linear combination of the scalar and
    the pseudoscalar for two spatial dimensions, and thus itself has two
    dimensions, is pure coincidence.  Hamilton was misled into believing
    that there must be a similar structure for ``d=3``; here the
    coincidence is that the number of vectors ``\binom{3}{1} = 3`` happens
    to equal the number of bivectors ``\binom{3}{2} = 3``.  These
    seemingly magical coincidences led to confusion that has only recently
    been resolved by the development of Geometric Algebra.

## Structure of the algebra

[Show diagrams of ``Cl(2)``, ``Cl(3)``, and ``Cl(3,1)`` in terms of
their bases, noting the numbers as binomial coefficients.]
