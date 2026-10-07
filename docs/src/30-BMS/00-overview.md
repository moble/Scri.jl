# [Overview](@id bms-overview)

## Motivation

This section provides a detailed description of the
Bondi-Metzner-Sachs (BMS) group, its action on the standard coordinate
tetrad, and its action on the fields defined with respect to that
tetrad — including the strain ``h`` and Weyl components such as
``Ψ₄``.

The BMS group is an *asymptotic* symmetry group of asymptotically flat
spacetimes in general relativity.  It consists of two parts: the
Lorentz transformations and supertranslations.  The supertranslations
generalize ordinary spacetime translations, and are essentially
direction-dependent time translations at null infinity — they allow us
to "re-slice" the asymptotic limits of spacetime with different
notions of constant-time slices.  While they are conceptually less
familiar, they turn out to be fairly simple to deal with in practice.
On the other hand, the Lorentz transformations are much more familiar,
but we need to examine them carefully to really understand how they
act on the tetrad and the fields.

The key fact about the Lorentz group that we will use *repeatedly* is
that by choosing a preferred timelike vector ``𝐭`` we can decompose
the Lorentz group into a product of four factors:

* ``𝕊²``: The sphere of spacelike directions ``k̂`` along which a
    future-directed null vector ``𝐤 = 𝐭 + k̂`` can point.
* ``M``: The one-dimensional subgroup of *rotations about* ``k̂``.
* ``A``: The one-dimensional subgroup of *boosts along* ``k̂``.
* ``N``: The two-dimensional subgroup of *null rotations* that leave
    ``𝐤`` invariant.

The ``M`` and ``A`` factors are directly related to the spin- and
boost-weights of the fields.  This decomposition is very easy to
calculate, and lets us do two very important things.

First, this decomposition allows us to understand the aberration that
simply relates *"where"* two observers are measuring a field — not
just the direction of the null ray representing the direction of
propagation, but also the rotation about that direction that
represents the polarization basis.  Our main goal will be to compute
the field measured by the second observer given the field measured by
the first observer, but to compute that we first need to know which
direction the second observer is looking in, and which polarization
basis they are using.  Essentially, these two pieces of information
are given by the ``𝕊²`` and ``M`` factors; the ``A`` and ``N``
factors drop out.

That information is necessary to compute the field, but it is not
sufficient because there is a major complication:

> Two different observers will actually use *two different tetrads* at
> each point.

Their tetrads are *not related* by the Lorentz transformation that
relates the observers.  The importance of this subtlety cannot be
understated, and is hidden by sloppy language frequently used in the
literature — words like "invariant" and "scalar" used to refer to
numbers that are actually nothing more than components of a tensor
with respect to a particular tetrad, which is itself defined naively
in terms of arbitrary and changeable coordinate systems.

The second capability the decomposition gives us is that it allows us
to relate the two tetrads used at each point.  Specifically, the
tetrads are each chosen to be pseudo-orthonormal, so we know that
there *exists* a Lorentz transformation relating them — though it is
distinct from the Lorentz transformation relating the observers.  To
distinguish them, we will refer to the Lorentz transformation relating
the observers as the **global Lorentz transformation**, and the
Lorentz transformation relating the tetrads at a given point as the
**local Lorentz transformation**.  The other important detail is that
the tetrads *at null infinity* are conventionally chosen to respect
the one invariant geometric feature at each point of null infinity:
the null direction along the generator of the surface.  Therefore, we
can decompose the local Lorentz transformation with respect to that
null direction — which shows us the *form* of one tetrad relative to
the other.  This makes the calculation of that relationship simple.

## Outline

All of this is made much easier if we use a formalism that is adapted
to the geometry while still being capable of expressing the algebra.
This, of course, is supplied by Geometric Algebra (GA).  The GA
formalism is surprisingly simple and elegant.  The algebraic
manipulations are just as easy as manipulating matrices in linear
algebra, with one simple but powerful addition: parallel vectors
commute and perpendicular vectors anti-commute.  And because the
algebra uses vectors directly, the geometric meaning of GA expressions
is transparent and directly tied to the physics of the situation.  The
[next page](@ref "A Primer on Geometric Algebra") provides a brief
introduction to the GA formalism.

We then move on to our discussion of [the Lorentz group](@ref "The
Lorentz Group").  The GA formalism allows us to easily manipulate
elements of the group, use them to operate on vectors and other
objects, and understand their structure.  This is also where we
introduce the decomposition of the Lorentz group described above.

Next, [we introduce the BMS group](@ref bms_group) itself.  One
important feature is the conformal factor ``κ``, which is not *part*
of the group per se, but is an important part of understanding how the
group acts on coordinates, tetrads, and fields.  We discuss the group
structure in detail, including its composition law.

We are [then](@ref "BMS Action on the Tetrad") in a position to
understand how observers related by a BMS transformation construct
their tetrads.  As noted above, this is a subtle but critical part of
understanding how the BMS group behaves.

Finally, we discuss [how the BMS group act on fields](@ref "BMS Action
on Fields") used to describe the physical phenomena of interest.  This
is a relatively straightforward and simple application of the
coordinate- and tetrad-transformation laws.
