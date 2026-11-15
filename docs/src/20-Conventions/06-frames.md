# [Frames and Components](@id frames_and_components)

Throughout this package we exclusively use passive transformations:
the spacetime events stay fixed, and the coordinates used to describe
them change.  Similarly, physical vectors stay fixed, and the basis
used to decompose them changes.  For simplicity, we will discuss the
conventions here as if we were dealing with a pseudo-Euclidean space,
assuming the naive relationship between coordinates and basis vectors:
the coordinates of an event are just the components of its position
vector relative to the origin, so they change opposite to the frame.

Perhaps most importantly, we will be using "rotors" from Geometric
Algebra (described in [the primer](@ref geometric_algebra_primer)) to
represent Lorentz transformations, so we also need to establish the
conventions for how these apply to vectors.  This page fixes what the
rotors ``𝐑`` and ``\boldsymbol{Λ}`` mean, and how they relate to the
velocity ``v⃗`` and conformal factor ``κ``.

Some of this page will seem elementary, belaboring basic facts.  But
the details of all these expressions are absolutely crucial for
getting things right.  The essential difficulty — and the reason we
need to treat this so carefully — is twofold.  First, the code doesn't
have access to the simple and elegant abstract-index notation; we need
to compute actual values, meaning that we need specific components in
a particular basis.  Second, the code stores a simple numerical array,
without any clear understanding of the basis in which those numbers
are expressed.  We need to be very careful to distinguish between the
possibilities in the code, and to make sure that the user can input
and read the values in a way that is consistent with the
documentation.

## Rotating vectors and rotating frames

We first assume an initial coordinate system and corresponding basis
frame ``(𝐭, 𝐱, 𝐲, 𝐳)``, which we denote more compactly as
``𝐞_μ``.  Here, the subscript ``μ ∈ (0,1,2,3)`` is just a label for
the basis vector, not a component of the vector.  We can then
introduce a second coordinate system and basis frame ``(𝐭', 𝐱', 𝐲',
𝐳')``, denoted ``𝐞'_μ``.  To relate the two frames, we use a Lorentz
transformation, represented by the rotor ``\boldsymbol{Λ}``, which
acts via the "sandwich" operation:

```math
𝐞'_μ = \boldsymbol{Λ}\,𝐞_μ\,\boldsymbol{Λ̃},
```

where the tilde denotes the reverse of the rotor, and juxtaposition
denotes the geometric product.  This equation is true for each choice
of ``μ``.

It is important to express a fixed vector in both frames — which
necessarily requires using components.  Index notation has a
convenient way of recording the frame in which the components of a
vector or more general tensor are expressed: by placing primes (or
sometimes hats or similar modifiers) on the indices themselves.[^1]
For example, we might have a "geometric" vector — meaning that it
represents a concept that is independent of any particular coordinate
system or basis — that we label ``𝐯``.  In this basis, ``𝐯`` has
components ``v^μ`` such that

```math
𝐯 = v^μ\,𝐞_μ
= v^0\,𝐭 + v^1\,𝐱 + v^2\,𝐲 + v^3\,𝐳.
```

We can then introduce a second basis ``𝐞'_{ν}``, which gives rise to
a second set of components ``v^{ν'}`` for the same vector:

```math
𝐯 = v^{ν'}\,𝐞'_{ν}
= v^{0'}\,𝐭' + v^{1'}\,𝐱' + v^{2'}\,𝐲' + v^{3'}\,𝐳'.
```

Note that we choose to put a prime directly on the ``𝐞'`` label *and*
on the index ``ν'``, but there is no prime directly on the ``v`` in
``v^{ν'}`` as that would suggest that the components are for a
distinct physical vector ``𝐯'``.[^2]  The components of our fixed
vector in the transformed basis follow logically:

```math
𝐯 = v^{μ'}\,𝐞'_μ = v^{μ'}\,\boldsymbol{Λ}\,𝐞_μ\,\boldsymbol{Λ̃}
\quad⟹\quad
\boldsymbol{Λ̃}𝐯\boldsymbol{Λ} = v^{μ'}\,𝐞_μ
\quad⟹\quad
v^{μ'} = \bigl(\boldsymbol{Λ̃}𝐯\boldsymbol{Λ}\bigr)^μ.
```

That is, the components relative to the new basis can be computed as
the components of the sandwich using the *reversed* rotor as
calculated entirely in the original frame.  Crucially, the last
equation is *not a geometric expression*; it is a purely numerical
statement.

[^1]: This was first introduced by Schouten [Schouten_1954](@cite) as
    the "kernel-index method", which is generally used with slight
    modifications throughout the rest of the literature.  It has also
    been supplemented throughout the literature with the
    abstract-index notation pioneered by Penrose; see Sec. 2.2 of
    [PenroseRindler_1984](@cite) for a thorough treatment.  Of course,
    abstract indices do not refer to any particular basis, and so
    never need primes to distinguish between bases.  Unfortunately,
    because we are doing explicit computations, we need to use the
    kernel-index method, and so we need to be careful to distinguish
    between the two bases.

[^2]: [Wald_1984](@citet) chose to include the prime on the
    *component* label ``v`` as well, presumably because ``ν'`` is a
    "dummy" (aka "bound") index, and Wald assumes that we must be free
    to relabel dummies arbitrarily.  [Schouten_1954](@citet)
    specifically advised against it — even going so far as to
    recommend a *prefixed* prime if needed to distinguish between
    physically distinct vectors.  We might thread the needle by
    suggesting that the index is just ``ν``, and the prime is a
    separate notation indicating the natural basis of the components,
    in which case it might be more correct to write
    ``v^{ν'}\,𝐞'_{ν}``.

[MisnerThorneWheeler_1973](@citet) (MTW) use similar notation for
Lorentz transformations (§2.9).  A single "kernel" letter ``Λ`` serves
for both transformations in both directions, and the positions of the
primes indicate which direction is meant: ``Λ^{μ'}{}_ν`` takes the old
components of a fixed vector to its new ones, ``v^{μ'} =
Λ^{μ'}{}_ν\,v^ν``, while the inverse matrix ``Λ^ν{}_{μ'}`` takes the
new components back to the old ones, and builds the new basis vectors
from the old, ``𝐞'_μ = 𝐞_ν\,Λ^ν{}_{μ'}``.  Comparing with ``𝐞'_μ =
\boldsymbol{Λ}\,𝐞_μ\,\boldsymbol{Λ̃}``, we see that our
``\boldsymbol{Λ}`` with the sandwich action is the same as their
action with ``Λ^ν{}_{μ'}``, and our ``\boldsymbol{Λ̃}`` is the same as
their ``Λ^{μ'}{}_ν``.  Rotor notation has no indices to record the
direction, so we must state once and for all which rotor each symbol
denotes.  In this package ``\boldsymbol{Λ}`` is always the rotor
applied to the basis vectors, as in ``𝐞'_μ =
\boldsymbol{Λ}\,𝐞_μ\,\boldsymbol{Λ̃}``, and we never use the symbol
``\boldsymbol{Λ}`` for its reverse.  The [comparison with other
sources](@ref frames_literature) below gives MTW's matrices for a
rotation and a boost, and relates both rotors to the Möbius
transformations used in the BMS literature.

## Spinors and rotors

In geometric algebra a general spinor ``𝐒`` acts on a vector by the
"sandwich" operation ``𝐯 ↦ 𝐒𝐯𝐒̃``, and the result is necessarily
another vector.  A *rotor* ``𝐑`` is a particular kind of spinor that
satisfies ``𝐑𝐑̃ = 1``, and therefore ``𝐑^{-1} = 𝐑̃`` — in which
case the sandwich operation preserves the norm of every vector, and so
represents an element of ``\mathrm{SO}⁺(3,1)``, while the rotors
themselves form ``\mathrm{Spin}⁺(3,1)``.  The sandwich always produces
a new vector from the one it is given.  Whether a calculation
describes an active or a passive transformation is determined only by
what the sandwich is applied to:

- Applied to a physical vector, with the basis held fixed, it moves
  the vector.  This is an *active* transformation.
- Applied to the basis vectors, with every physical vector held fixed,
  it builds a new basis ``𝐞'_μ = 𝐑𝐞_μ𝐑̃``.  This is a *passive*
  transformation, or change of frame, and ``𝐑`` is then called the
  *frame rotor*.

We saw above that the new components of a fixed vector are then those
of ``𝐑̃𝐯𝐑``.  That is, the sandwich by the *reverse* of the frame
rotor takes the old components of a vector to its new ones — and,
since the coordinates of an event are the components of its position
vector, the old coordinates of an event to its new ones.  The
numerical array of components does not record the basis, so the same
numbers admit more than one reading.  In particular, an active
transformation by ``𝐑̃`` and a change of frame whose frame rotor is
``𝐑`` produce identical components:

| Reading of ``𝐑``      | What ``𝐑`` does                        | New components from old ``𝐯`` |
|:----------------------|:---------------------------------------|:------------------------------|
| active transformation | moves physical vectors; basis fixed    | ``𝐑𝐯𝐑̃``                       |
| frame rotor           | moves the basis vectors; vectors fixed | ``𝐑̃𝐯𝐑``                       |

## The frame rotor

The following statements fix the conventions.  Everything else on the
convention and BMS pages is derived from them.

1. The signature is ``{-}{+}{+}{+}``, and the geometric algebra is
   taken over ``ℝ``.
2. A new frame is related to the old one by the frame rotor
   ``\boldsymbol{Λ}``, applied to the legs as ``𝐞'_μ =
   \boldsymbol{Λ}\,𝐞_μ\,\boldsymbol{Λ̃}``.
3. The new frame is specified by a boost velocity ``v⃗`` and a spatial
   rotation ``𝐑``, with

   ```math
   \boldsymbol{Λ} = 𝐁(v⃗)\,𝐑,
   ```

   where ``𝐁(v⃗)`` is the boost rotor `Boost(v⃗)`.  ``v⃗`` is the
   velocity of the new observer relative to the old one, in the old
   frame's components, with the two origins coinciding at ``t = t' =
   0``.  In particular, ``\boldsymbol{Λ} = 𝐑`` when ``v⃗ = 0``.
4. The new components of a fixed vector ``𝐯`` are the old components
   of ``\boldsymbol{Λ̃}𝐯\boldsymbol{Λ}``, and the new coordinates of
   an event are the old components of
   ``\boldsymbol{Λ̃}𝐱\boldsymbol{Λ}``, where ``𝐱`` is its position
   vector from the common origin.

The order of the factors in item 3 shows that ``v⃗`` is a vector in
the old frame's components.  The new time leg is

```math
𝐭' = \boldsymbol{Λ}\,𝐭\,\boldsymbol{Λ̃}
= 𝐁(v⃗)\,𝐑\,𝐭\,𝐑̃\,𝐁̃(v⃗)
= 𝐁(v⃗)\,𝐭\,𝐁̃(v⃗)
= γ(𝐭 + v⃗),
```

because ``𝐑`` leaves ``𝐭`` unchanged.  The other order would give

```math
𝐑\,𝐁(v⃗)\,𝐭 \,𝐁̃(v⃗)\,𝐑̃ = γ(𝐭 + 𝐑v⃗𝐑̃).
```

Equivalently, the legs are first rotated by ``𝐑`` and then boosted by
``v⃗``, with both specified in the old frame's components.

A `BMS` element stores the frame rotor ``\boldsymbol{Λ}``; the
accessor [`lorentz`](@ref Scri.lorentz) returns that full spacetime
rotor, and [`boost_velocity`](@ref Scri.boost_velocity) and
[`frame_rotation`](@ref Scri.frame_rotation) recover ``v⃗`` and
``𝐑``.  The function [`transform!`](@ref Scri.transform!) can
directly take ``v⃗`` and ``𝐑`` defined in this way (among other
arguments).

### Some crucial applications

For a pure rotation the frame rotor is just ``\boldsymbol{Λ} = 𝐑``,
and the new axes are the old axes turned by ``𝐑``.  For example,
`Rotor(√2, 0, 0, √2)` is a quarter turn about ``𝐳``, and puts ``𝐱'``
along the old ``𝐲`` and ``𝐲'`` along the old ``-𝐱``.  The event
with old coordinates ``(t, x, y, z) = (0, 1, 0, 0)``, whose position
vector is ``𝐱 = -𝐲'``, therefore has new coordinates ``(0, 0, -1,
0)``: the coordinates turn the opposite way to the axes, as item 4
says.

For a pure boost, ``\boldsymbol{Λ} = 𝐁(v⃗) = \cosh(φ/2) +
𝐭v̂\,\sinh(φ/2)``, with ``v⃗ = \tanh φ\,v̂`` and ``γ = \cosh φ``.
The bivector ``𝐭v̂`` squares to ``+1`` and anticommutes with both
``𝐭`` and ``v̂``.  Then item 2 gives the new time leg,

```math
𝐭' = \boldsymbol{Λ}\,𝐭\,\boldsymbol{Λ̃} = \boldsymbol{Λ}^2\,𝐭
= (\cosh φ + 𝐭v̂\,\sinh φ)\,𝐭 = γ(𝐭 + v⃗),
```

which is the statement of item 3 that the new observer moves at
``v⃗``.  Item 4 gives the new time component of the null vector ``𝐤 =
𝐭 + k̂``.  The part of ``k̂`` perpendicular to ``v̂`` commutes with
``\boldsymbol{Λ}`` and has no time component, while ``𝐭`` and ``v̂``
anticommute with ``𝐭v̂``, so

```math
\bigl(\boldsymbol{Λ̃}\,𝐤\,\boldsymbol{Λ}\bigr)^0
= \bigl(\boldsymbol{Λ̃}^2\,𝐭\bigr)^0
  + (v̂⋅k̂)\,\bigl(\boldsymbol{Λ̃}^2\,v̂\bigr)^0
= \cosh φ - (v̂⋅k̂)\,\sinh φ
= γ\,(1 - v⃗⋅k̂).
```

This is the inverse of the conformal factor ``κ`` defined below, and
computed as `κ⁻¹ᵢ` in the implementation of `transform!`.

## [Comparison with other sources](@id frames_literature)

Most sources describe Lorentz transformations by matrices, or — in the
BMS literature — by Möbius transformations of the sphere, rather than
by rotors.  To translate, we need to know two things about each
source: whether its transformations are active or passive, and in
which direction each of its maps goes.  Fortunately, the sources we
rely on are explicit about both.

### The matrices of MTW

MTW treat Lorentz transformations as passive.  Of their two matrices,
``Λ^{μ'}{}_ν`` and ``Λ^ν{}_{μ'}``, "the first transforms coordinates
from an unprimed frame to a primed frame, while the second goes from
primed to unprimed" (§2.9), and the examples in their Box 2.4 are the
"Rotation of Frame of Reference" and the "Boost of Frame of
Reference."  These are the two applications above.  Box 2.4 marks the
new coordinates with bars; we write primes instead.

For a rotation by an angle ``θ`` in the ``x``-``y`` plane, Box 2.4
gives

```math
x = x' \cos θ - y' \sin θ,
\qquad
y = x' \sin θ + y' \cos θ,
```

so the new ``x'`` axis lies at the angle ``θ`` from the old ``x``
axis, toward ``y``.  This is the frame rotor ``\boldsymbol{Λ} = 𝐑``
of a turn by ``θ`` about ``𝐳``, `Rotor(cos(θ/2), 0, 0, sin(θ/2))`,
and these equations give the old coordinates of an event from its new
ones.  With ``θ = π/2`` they reproduce the quarter turn above.

For a boost along ``z``, the signs in Box 2.4 are chosen so that an
"object at rest at ``z' = 0`` in rocket frame moves in direction of
increasing ``z`` in lab frame."  Writing our rapidity ``φ`` in place
of their velocity parameter, they give

```math
\bigl\|Λ^μ{}_{ν'}\bigr\| =
\begin{pmatrix}
\cosh φ & 0 & 0 & \sinh φ \\
0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0 \\
\sinh φ & 0 & 0 & \cosh φ
\end{pmatrix},
```

and ``\bigl\|Λ^{μ'}{}_ν\bigr\|`` is the same matrix with ``-\sinh φ``
in place of ``\sinh φ``.  Since ``𝐞'_μ = 𝐞_ν\,Λ^ν{}_{μ'}``, the
columns of ``\|Λ^μ{}_{ν'}\|`` are the old components of the new legs.
The first column is the new time leg ``\cosh φ\,𝐭 + \sinh φ\,𝐳 =
γ(𝐭 + v⃗)``, which we found above for the frame rotor
``\boldsymbol{Λ} = 𝐁(v⃗)`` with ``v⃗ = \tanh φ\,𝐳``.  So
``\|Λ^μ{}_{ν'}\|`` is the matrix of our ``𝐁(v⃗)``, and
``\|Λ^{μ'}{}_ν\|`` is the matrix of ``𝐁̃(v⃗)``.

### Möbius transformations in the BMS literature

At ``ℐ⁺``, each point of the sphere is the direction ``k̂`` of an
outgoing null ray ``𝐤 = 𝐭 + k̂``, and the BMS literature usually
labels it by the stereographic coordinate ``ζ = e^{iϕ} \cot(θ/2)``.  A
Lorentz transformation then acts on these labels as a Möbius
transformation,

```math
ζ ↦ \frac{aζ + b}{cζ + d},
\qquad
ad - bc = 1.
```

As with the matrices, the direction of such a map depends on whether
it is read actively or passively.

[Sachs_1962b](@citet) (§II), [NewmanPenrose_1966](@citet) (§II), and
[Boyle_2015](@citet) (§II A) read it passively.  Sachs "subjects the
coordinates … to a Lorentz transformation," and finds that the angles
at infinity undergo a conformal transformation, which he writes as a
Möbius transformation of ``z = \cot(θ/2)\,e^{iϕ}``.  Newman and
Penrose define the BMS group by "transformations on the ``θ``, ``ϕ``,
``u`` coordinates."  And in Boyle, the stereographic coordinates
``ζ_A`` and ``ζ_B`` of the same point in frames ``A`` and ``B`` are
related by ``ζ_B = (aζ_A + b)/(cζ_A + d)``.  In each case the map
takes the old label of a fixed null direction to its new label.  So it
is the action of ``\boldsymbol{Λ̃}``: the new label of the ray ``𝐤``
is the direction of the spatial part of
``\boldsymbol{Λ̃}\,𝐤\,\boldsymbol{Λ}``.

[PenroseRindler_1984](@citet) read it actively: "If not further
qualified, Lorentz and Poincaré transformations in this book will be
understood to be active" (§1.1).  Their Möbius transformation (their
Eq. 1.2.17) moves each null direction to its image, with both labeled
in a single tetrad.  They also show that "the active Lorentz
transformation ``L`` that carries ``g_{\mathbf{i}}`` into
``g_{\hat{\mathbf{i}}}`` is *formally* equivalent, in its effect on
the coordinates of a vector, to the passive Lorentz transformation
``G^{-1}``," where ``G`` is the passive transformation induced by
passing from the old tetrad ``g_{\mathbf{i}}`` to the new tetrad
``g_{\hat{\mathbf{i}}}``.  This is the duality in the table above.
The active transformation that carries our old tetrad into the new one
is ``\boldsymbol{Λ}`` itself, so for a given change of frame, their
Möbius transformation is the inverse of the passive one in the other
three sources.  It takes the new label of a null direction to its old
label, which is exactly what `transform!` needs in order to find where
to read the old data.

This map is implemented as part of the action of a [`BMS`](@ref
Scri.BMS) element on the points of ``ℐ``.  The element `g` is a
functor: calling it as [`g(t′, k̂′)`](@ref Scri.BMS(::Real,
::Scri.QuatVec) where {T<:Real,A,I}) returns the old coordinates ``(t,
k̂)`` of the point whose new coordinates are ``(t', k̂')``.  The
direction is given by Penrose and Rindler's Möbius transformation, as
the old ray ``𝐤 ∝ \boldsymbol{Λ}\,𝐤'\,\boldsymbol{Λ̃}`` with ``𝐤' =
𝐭 + k̂'``, and the time by inverting the time law [below](@ref "The
time coordinate"), as ``t = t'/κ + c_α\,α(k̂)``.  The inverse map,
from old coordinates to new, is computed by the unexported function
[`Scri.new_coordinates`](@ref), and its angular part is the passive
Möbius transformation of the other three sources.

Two examples will help fix the direction of the passive map.  In the
quarter turn above, the azimuth of every fixed direction decreases by
``π/2``, so ``ζ' = -i\,ζ``; the old ``𝐱`` direction, with ``ζ = 1``,
has ``ζ' = -i``, which is the direction of ``-𝐲'``.  For a boost with
rapidity ``φ`` along ``𝐳``, the angle between a fixed ray and the
direction of motion grows, with ``\tan(θ'/2) = e^{φ} \tan(θ/2)`` as in
Eq. (1) of [Boyle_2015](@citet), so that ``ζ' = e^{-φ}\,ζ``.

In summary, for a given change of frame:

| Our rotor          | Coordinates | MTW            | Penrose and Rindler | Möbius transformation of Sachs, Newman-Penrose, and Boyle |
|:-------------------|:------------|:---------------|:--------------------|:----------------------------------------------------------|
| ``\boldsymbol{Λ}`` | new to old  | ``Λ^ν{}_{μ'}`` | active ``L``        | its inverse, ``ζ' ↦ ζ``                                   |
| ``\boldsymbol{Λ̃}`` | old to new  | ``Λ^{μ'}{}_ν`` | passive ``G``       | ``ζ ↦ ζ'``                                                |

### The active transformations in `Quaternionic` and `SphericalFunctions`

The documentation in this package consistently uses passive
transformations.  However, the two main packages that `Scri` builds
on, `Quaternionic` and `SphericalFunctions`, generally use active
transformations in their documentation.  Numerically, the new
components of a fixed vector are given by the *active* transformation
by ``\boldsymbol{Λ̃}`` that `Quaternionic` and `SphericalFunctions`
apply to arrays of components — as, for example, `Λ̃(V)` and `D(Λ̃,
ℓₘₐₓ) * f` respectively.  That is why the documentation of
`Quaternionic` calls `inv(Boost(v⃗))` "the passive transformation" for
an observer moving with velocity ``v⃗``: it is ``\boldsymbol{Λ̃}`` for
``𝐑 = 1``.

## Composition and translating geometry into code

Above, we glossed over one important subtlety: the basis in which
*rotors* are expressed.  Implicitly, expressions above all use the
original frame's basis to express the components of the rotors.  But
when we apply a second change of frame, we need to be careful about
the basis in which the second rotor's components are expressed.

We assume that there are three vector-basis frames: ``𝐞_μ``,
``𝐞'_μ``, and ``𝐞''_μ``.  We further assume that they are related by

```math
𝐞'_μ = \boldsymbol{Λ}₁\,𝐞_μ\,\boldsymbol{Λ̃}₁
\quad \text{and} \quad
𝐞''_μ = \boldsymbol{Λ}₂\,𝐞'_μ\,\boldsymbol{Λ̃}₂.
```

We also see that

```math
𝐞''_μ
= \boldsymbol{Λ}₂\, 𝐞'_μ\, \boldsymbol{Λ̃}₂
= \boldsymbol{Λ}₂\, \boldsymbol{Λ}₁\, 𝐞_μ\, \boldsymbol{Λ̃}₁\,
  \boldsymbol{Λ̃}₂.
```

Thus, the rotor that carries the original frame to the final frame is
``\boldsymbol{Λ}₂\boldsymbol{Λ}₁``.

There are some important subtleties when we try to code up the
elegant, coordinate-independent geometric expressions above.  For
example, a vector is represented in the code as — essentially — a
wrapper around an array of four numbers, and a rotor as a wrapper
around an array of eight (real or four complex) numbers.  But the code
does not keep track of the basis in which those numbers are expressed.
So we need to design some translation rules.  Essentially the
particular basis is irrelevant; the code just assumes a *consistent*
basis for any given calculation.  Thus, our first rule is for the
simple case: *if* each factor in an expression is encoded in the same
basis, then we essentially copy the geometric expression into the
code, and the result will be correct.  For example, if our geometric
expression is

```math
𝐞''_μ
= \boldsymbol{Λ}₂\, \boldsymbol{Λ}₁\, 𝐞_μ\, \boldsymbol{Λ̃}₁\,
  \boldsymbol{Λ̃}₂,
```

and we have corresponding objects in the code `Λ₂` and `Λ₁` expressed
in the same basis, then we can compute this as

```julia
𝐞′′ = Λ₂ * Λ₁ * 𝐞 * conj(Λ₁) * conj(Λ₂).
```

However, if — as is frequently the case — `Λ₂` is known only in the
basis of the primed frame, then we need to rewrite its components to
be valid in the unprimed frame.  Because the code assumes a single
frame for the calculation, it will interpret the components `Λ₂` as if
they were appropriate for the unprimed frame, which they are not.  We
have to explicitly rotate those incorrectly interpreted components
back to the primed frame: the correct code expression in this case is
just `Λ₁ * Λ₂ * conj(Λ₁)`.  We can then compute

```julia
𝐞′′ = Λ₁ * Λ₂ * 𝐞 * conj(Λ₂) * conj(Λ₁).
```

Essentially, we have reversed the order of the factors relative to the
previous expression.  So we simply have to be careful about the basis
in which each factor is expressed.

Generally, in this code, we assume that successive transformations are
expressed in *successive* frames.  For example, if `g₁` and `g₂` are
[`BMS`](@ref Scri.BMS) elements, then we assume that `g₁` is expressed
in the original frame, and `g₂` is expressed in the frame produced by
`g₁`, so that their composition should be written as `g₁ * g₂`.
Similarly, [`compose`](@ref Scri.compose)`(g₁, g₂)` applies `g₁` first
and then `g₂`, and `lorentz(g₁ * g₂) ≈ lorentz(g₁) * lorentz(g₂)`.
The main motivation is that the [`transform!`](@ref Scri.transform!)
function behaves in this way: a call using `g₁` followed by a call
using `g₂` must necessarily interpret the second call's argument as
expressed in the frame produced by the first transformation.  The
result is equivalent to a single call using `g₁ * g₂`.

The `BMS` type is meant to represent the *abstract* BMS group.  The
group law is written in this "moving-frame" reading because it is the
composition of the changes of frame: ``g₂`` acts on the frame that
``g₁`` produces.  This is clearest in the composed supertranslation,
``α₁(𝐤) + α₂(\boldsymbol{Λ̃}₁𝐤\boldsymbol{Λ}₁)/κ(\boldsymbol{Λ}₁,
𝐤)``, where ``α₂`` is evaluated at frame 1's label of the ray, and
the factor ``1/κ`` converts the shift of frame 1's time into a shift
of the original time.  The same order makes the action of an element
on the points of ``ℐ`` compose like a function: `g(t′, k̂′)` returns
the old coordinates of the point whose new coordinates are ``(t',
k̂')``, and `(g₁ * g₂)(t′, k̂′) == g₁(g₂(t′, k̂′)...)`.

## [The conformal factor ``κ``](@id conformal_factor_convention)

A point of ``ℐ`` is labeled by the null ray ``𝐤 = (1, ℐ\,k̂)``, where
``ℐ = ±1`` selects [future or past null infinity](@ref
scri_pm_conventions).  The conformal factor is defined as the ratio of
the old to the new time components of a selected null ray,

```math
κ(\boldsymbol{Λ}, 𝐤) = \frac{k^0}{k^{0'}}
= \frac{k^0}{(\boldsymbol{Λ̃}\,𝐤\,\boldsymbol{Λ})^0}
= \frac{1}{γ\,(1 - ℐ\,v⃗⋅k̂)},
```

where the last form follows from the application calculated above, and
is independent of ``𝐑``.  We name the frame rotor in
``κ(\boldsymbol{Λ}, 𝐤)`` when it matters, and this is the number that
[`conformal_factor`](@ref Scri.conformal_factor)`(Λ, k̂)` returns.
For a fixed ray, intervals of time scale as ``Δt' = κ\,Δt``, so ``κ``
is a redshift factor: frequencies scale as ``ω' = ω/κ``, and ``κ = 1 +
z`` for a source at rest in the old frame.  Near ``ℐ⁺``, an observer
moving along the direction of propagation sees the radiation
redshifted, ``κ > 1``.  In general the radiation is redshifted in
exactly the directions with ``v⃗⋅k̂ > 1 - 1/γ``, because time dilation
alone blueshifts the directions transverse to the motion.  The new
label of the ray is read from ``𝐤' =
\boldsymbol{Λ̃}\,𝐤\,\boldsymbol{Λ}/(\boldsymbol{Λ̃}\,𝐤\,\boldsymbol{Λ})^0``,
and for two changes of frame in succession

```math
κ(Λ₁Λ₂, 𝐤) = κ(Λ₂, Λ̃₁𝐤Λ₁)\,κ(Λ₁, 𝐤).
```

[Sachs_1962a](@citet) — the creator of the full BMS group — and
[NewmanPenrose_1966](@citet) use the symbol ``K`` for this factor; we
follow later references in using ``κ`` instead, to avoid confusion
with the ``K`` of the [Iwasawa decomposition](@ref "Iwasawa's ``KAN``
decomposition") used throughout these docs.

## The time coordinate

A BMS element is represented as a pair ``(Λ, α)``.  The
supertranslation ``α`` is applied first, on the old sphere, and the
Lorentz transformation second, about the new origin set by the
supertranslation rather than the original one:

```math
t' = κ\,\bigl[t - c_α\,α(k̂)\bigr],
```

where ``α`` is a function of the old direction and ``c_α`` is the
[sign convention](@ref conventions-type) of the supertranslation.  The
same law holds at ``ℐ⁺``, where ``t`` represents the retarded time
``u``, and at ``ℐ⁻``, where it represents the advanced time ``v``.
With products written in the order of application, two elements
compose as

```math
(Λ₁, α₁)\,(Λ₂, α₂) = \bigl(Λ₁Λ₂,\ α₁ + \mathrm{Ad}(Λ₁)(α₂)\bigr),
\qquad
\mathrm{Ad}(Λ)(α)(𝐤) = \frac{α(Λ̃\,𝐤\,Λ)}{κ(Λ, 𝐤)},
```

as derived on [the group page](@ref bms_group), and the factor ``c_α``
cancels from the law of composition.

### Comparison with Sachs

[Sachs_1962a](@citet) (§3b) orders the two parts the other way.  With
primes on the new coordinates, his equations (3.5), (3.7), and (3.12)
amount to

```math
r = K\,r' + J,
\qquad
K^2\,dΩ^2 = dΩ'^2,
\qquad
u = K^{-1}\,\bigl[u' + α_{\text{Sachs}}(θ', ϕ')\bigr].
```

The first two show that his ``K`` is our ``κ``, not its reciprocal.
The third reads ``u' = K\,u - α_{\text{Sachs}}``: the Lorentz part is
applied first, and the supertranslation is then a function on the
*new* sphere, with no factor of ``K``.  Both orderings parameterize
the same group.  Comparing ``u' = κu - c_ακα`` with Sachs's form at a
ray whose old label is ``k̂`` and new label is ``k̂'``,

```math
α_{\text{Sachs}}(k̂') = c_α\,κ(k̂)\,α(k̂),
```

so Sachs's element is our ``\bigl(Λ,\ c_α\,(α_{\text{Sachs}}∘Λ̃)/κ
\bigr)``, where ``α_{\text{Sachs}}∘Λ̃`` means ``α_{\text{Sachs}}``
evaluated at the new label of the ray.  Sachs's definition is unusual
in the literature, having been quickly replaced by the form used here.
For example, [NewmanPenrose_1966](@citet) (their Eq. 2.1) write ``u' =
K\,(u - α)``, with both ``K`` and ``α`` functions of the old angles,
which is our law with ``c_α = +1``.

## Fields on the sphere

`SphericalFunctions` treats a spin-weighted function as a function of
a *point frame*: a rotor ``Q`` that takes the reference basis ``(𝐱,
𝐲, 𝐳)`` to a frame at a point of the sphere.  Under a pure rotation,
``\boldsymbol{Λ} = 𝐑``, and a field ``f`` is written in the new frame
as

```math
f'(Q) = f(ΛQ),
\qquad
f' = 𝔇(Λ̃)\cdot f,
```

pointwise and on its mode weights, respectively.  That is, the point
frame whose components in the new frame are ``Q`` has components
``ΛQ`` in the old frame.  The mode law is
`SphericalFunctions`'s [rotation law](@extref SphericalFunctions
:std:label:`spin_rotation_law`) with the rotor ``\boldsymbol{Λ̃}``,
which is again the duality between an active rotation by
``\boldsymbol{Λ̃}`` and a change of frame by ``\boldsymbol{Λ}``.

In general the new field at a new pixel ``R'_p`` is read from the old
field at the point frame

```math
R_p = K\bigl(Λ\,R'_p\bigr),
```

where ``K(\,⋅\,)`` is the rotation factor of the Iwasawa
decomposition, and the value is then multiplied by suitable powers of
``κ`` and mixed with other components as described in the following
sections.  The first part of this process is performed by
`aberration`, which is described in detail on [its own page](@ref
"Aberration of Gravitational Waves").

At ``ℐ⁻`` the labeling of the sphere is antipodal.  The implementation
then uses the parity conjugate of the frame rotor, ``\boldsymbol{Λ}_ℐ
= 𝐁(-v⃗)\,𝐑``, which is `ℂconj(Λ)` and realizes the frame rotor on
the antipodal labeling; the code calls it `Λℐ`, which is just `Λ`
itself at ``ℐ⁺``.

## Symbols

| Symbol               | Meaning                                                                    | In the code                              |
|:---------------------|:---------------------------------------------------------------------------|:-----------------------------------------|
| ``v⃗``                | velocity of the new observer, in old components                            | `v⃗`, `boost_velocity(g)`                 |
| ``𝐑``                | spatial rotation of the frame                                              | `R`, `frame_rotation(g)`                 |
| ``\boldsymbol{Λ}``   | frame rotor ``𝐁(v⃗)𝐑``, with ``𝐞'_μ = \boldsymbol{Λ}𝐞_μ\boldsymbol{Λ̃}``     | `g.Λ`, `lorentz(g)`                      |
| ``\boldsymbol{Λ̃}``   | its reverse; new components are those of ``\boldsymbol{Λ̃}𝐯\boldsymbol{Λ}`` | `inv(g.Λ)`                               |
| ``\boldsymbol{Λ}_ℐ`` | ``\boldsymbol{Λ}`` at ``ℐ⁺``, and `ℂconj(Λ)` at ``ℐ⁻``                     | local `Λℐ`                               |
| ``κ``                | conformal factor ``k^0/k^{0'}``, a redshift                                | `conformal_factor`; `κ⁻¹ᵢ` holds ``1/κ`` |
| ``α``                | supertranslation, on the old sphere                                        | `g.α`, `supertranslation(g)`             |
| ``c_α``              | sign convention of the supertranslation                                    | `c_α`                                    |
| ``ℐ``                | ``+1`` at ``ℐ⁺``, ``-1`` at ``ℐ⁻``                                         | `ℐ`                                      |
| ``R'_p``, ``R_p``    | point frame of a pixel, in new and old components                          | `R′ₚ`, `Rₚ`                              |
| ``K``, ``A``, ``N``  | factors of the Iwasawa decomposition                                       | `Quaternionic.KAN`                       |
