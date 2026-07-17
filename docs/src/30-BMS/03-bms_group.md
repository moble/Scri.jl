# [The BMS Group](@id bms_group)

## Bondi-Sachs coordinates, metric, and gauge

The natural setting for outgoing radiation is a coordinate system built
from the light cones themselves.  We foliate the neighborhood of ``ℐ⁺``
by outgoing null cones, label each cone by a *retarded time* ``u`` that
is constant across it, and coordinatize the cross-sections by a radial
coordinate ``r`` together with two angles ``xᴬ = (θ, φ)`` that stay
fixed along each outgoing null ray.  That ray — the generator of the
cone — is the null vector ``ℓ`` of [the standard tetrad](@ref
"Tetrad").  These are the *Bondi–Sachs coordinates* [Sachs_1962a](@cite).

*Bondi gauge* is the handful of conditions that pin them down.  We take
``r`` to be an areal radius, so that the determinant of the angular
metric is ``r⁴`` times that of the round unit sphere, and we require the
metric to approach Minkowski as ``r → ∞``, with its angular part
approaching ``r² dΩ²`` — the round unit-sphere metric.  The full
asymptotic expansion of the metric is standard [Sachs_1962a](@cite) and
thoroughly coordinate-dependent; we will not need it.  We need exactly
one fact drawn from it: *asymptotically, the angular part of the metric
is ``r²`` times a round unit-sphere metric.*

That single condition is what gives the conformal factor its meaning.
An asymptotic symmetry is free to rescale ``r`` and to change the
round metric on the sphere, but Bondi gauge ties the two together: the
sphere metric may change only by an overall, angle-dependent
rescaling, and the rescaling of the radius against the rescaling of
the sphere is precisely Sachs's conformal factor ``κ``.  The rest of
this page makes that statement precise by dropping the coordinates and
keeping the geometry — the cross-sections of ``ℐ⁺`` assemble into a
*celestial sphere*, and the asymptotic symmetries act on it
conformally.

## The null cone and the celestial sphere

It is cleanest to build the picture first in flat spacetime and lift it
to ``ℐ⁺`` afterward.  Start, then, with the future null cone ``𝒩⁺`` of
the origin — the future-pointing null vectors ``𝐤``, with the metric
the flat ``ηᵦᵧ``.  What matters for radiation is not an individual null
vector but its *direction*: two vectors differing only by a positive
rescaling ``𝐤 ↦ λ𝐤`` describe the same outgoing ray.  The space of
rays is the **celestial sphere**,

```math
𝕊² = 𝒩⁺ / ℝ⁺,
```

the quotient of the cone by that rescaling.  Equivalently, ``𝒩⁺`` is a
fiber bundle over ``𝕊²`` with fiber ``ℝ⁺``, whose projection

```math
π(𝐤) = [𝐤] = \{ λ𝐤 \mid λ > 0 \}
```

sends each null vector to its ray.  This bundle is the whole of the
geometric content; everything that follows is a matter of how we
*represent* it concretely.

A concrete sphere needs a *section* — a chosen representative null
vector for each ray — and a unit timelike vector ``𝐭`` supplies one.
Normalizing ``g(𝐭, 𝐭) = -1``, the projection

```math
π_𝐭(𝐤) = \frac{𝐤}{-g(𝐤, 𝐭)}
```

rescales each null vector ``𝐤`` to have unit time component in the
``𝐭`` frame, and the spatial unit vector

```math
k̂(𝐤) = π_𝐭(𝐤) - 𝐭
```

then gives its direction on the sphere.  Introducing coordinates
``xᴬ`` on ``𝕊²`` — the usual ``(θ, φ)``, say — the section is

```math
σᵝ(xᴬ) = (1, k̂(xᴬ))ᵝ,
```

the future null ray with unit time component in the ``𝐭`` frame.
This is the ``ℓ`` of [the standard tetrad](@ref "Tetrad") (up to
conventions):

```math
ℓᵝ = -ηᵝᵞ\frac{(dt-dr)ᵧ}{\sqrt{2}} \sim -gᵝᵞ\frac{(du)ᵧ}{\sqrt{2}}.
```

The section drags the round metric of the sphere along with it.  Pulling
the flat metric back through ``σ`` gives the induced metric on ``𝕊²``,

```math
dΩ² = σ^* η = ηᵦᵧ \frac{∂σᵝ}{∂xᴬ} \frac{∂σᵞ}{∂xᴮ} dxᴬ dxᴮ,
```

or, in components,

```math
g_{AB} = ηᵦᵧ \frac{∂σᵝ}{∂xᴬ} \frac{∂σᵞ}{∂xᴮ}.
```

For the standard section, this is just the round metric.  Writing it
as a pullback is what lets us read off, in the next section, exactly
what a Lorentz transformation does to it.

## Lorentz transformations ``ℒ`` and the conformal factor ``κ``

A Lorentz transformation ``Λ`` preserves the null cone ``𝒩⁺`` — null
vectors map to null vectors — but it does *not* preserve the section:
the image ``Λσ`` of a unit-time-component null vector generally has some
other time component, so it no longer satisfies ``ℓ⁰ = 1``.  Restoring
the preferred section takes a rescaling, and the factor that does it is
the conformal factor.  Define

```math
κ(Λ, 𝐤) = \frac{k⁰}{{Λ⁰}ᵦkᵝ} = \frac{1}{γ(1-v⃗⋅k̂)},
```

the ratio of the old time component to the new one; the second form is
the familiar Doppler factor of a pure boost with velocity ``v⃗``,
recovered explicitly on the [aberration page](@ref "Aberration of
Gravitational Waves").  A Lorentz transformation acts on the section by

```math
σ' = κ\, (Λ ∘ σ),
```

the factor of ``κ`` undoing the change in time component, so that
``σ'`` again has ``ℓ'^{0'} = 1``.

Now, by the same logic as in the previous section, we can show the
pullback in the new coordinate system:

```math
g'_{AB} = ηᵦᵧ \frac{∂{σ'}ᵝ}{∂xᴬ} \frac{∂{σ'}ᵞ}{∂xᴮ}.
```

Applying ``∂/∂xᴬ`` to ``{σ'}ᵝ = κ\, {Λᵝ}ᵧσᵞ`` produces two terms — one
differentiating ``κ`` and one differentiating ``σ``:

```math
\frac{∂{σ'}ᵝ}{∂xᴬ} = \frac{∂κ}{∂xᴬ} {Λᵝ}ᵧσᵞ + κ {Λᵝ}ᵧ \frac{∂σᵞ}{∂xᴬ}.
```

Three facts collapse the expansion.  First, ``Λ`` preserves the metric:

```math
ηᵦᵧ {Λᵝ}ᵤ {Λᵞ}ᵥ = ηᵤᵥ.
```

Second, the section is null:

```math
ηᵦᵧ σᵝ σᵞ = 0.
```

And third, differentiating that relation shows that the section is
orthogonal to its own derivative:

```math
ηᵦᵧ σᵝ \frac{∂σᵞ}{∂xᴬ} = 0.
```

Every ``∂κ/∂xᴬ`` term multiplies either ``σ·σ`` or ``σ·∂σ`` and so
drops out, leaving only the original metric scaled by ``κ²``:

```math
g'_{AB} = κ² ηᵦᵧ \frac{∂{σ}ᵝ}{∂xᴬ} \frac{∂{σ}ᵞ}{∂xᴮ} = κ² g_{AB}.
```

That is, ``{dΩ'}² = κ² dΩ²``.  This is the key result: **a Lorentz
transformation acts on the celestial sphere by a conformal
rescaling**, with conformal factor ``κ²``.

This brings us back to Bondi gauge, which fixes the angular metric to
asymptotically have the form ``r² dΩ²`` in any frame; if the
transformation is to be an asymptotic isometry of the one physical
spacetime, the two frames must agree there: ``r² dΩ² \sim {r'}²
{dΩ'}²``.  With ``{dΩ'}² = κ² dΩ²`` this requires that

```math
r \sim κ\, r',
```

which is exactly Sachs's *definition* of ``κ``.  The geometric
conformal factor of the sphere and the coordinate rescaling of the
radius are one and the same.

!!! details "A coordinate derivation of the conformal factor"

    Above, we simply defined ``κ``, and then showed that it is the
    conformal factor of the sphere metric under a Lorentz
    transformation. The same result can be obtained more mechanically,
    in explicit coordinates.  It is less illuminating but more
    familiar, so it is worth recording.

    Take the boost to be along ``𝐳``, so that the polar angle obeys
    the aberration formula

    ```math
    \cos θ' = \frac{\cos θ - β}{1 - β \cos θ},
    ```

    and ``κ`` as defined above is expressed as

    ```math
    κ = \frac{\sqrt{1-β²}}{1-β\cos θ}.
    ```

    Differentiating the aberration formula gives ``\sin θ'\, dθ' = κ²
    \sin θ\, dθ``, while ``\sin² θ' = 1 - \cos² θ'`` applied to the
    aberration formula gives ``\sin θ' = κ \sin θ``; together these
    yield ``dθ' = κ\, dθ``.  The unit-sphere metric ``dθ² + \sin²θ\,
    dφ²`` therefore rescales as ``{dΩ'}² = κ²\, dΩ²``.  Rotations
    leave the unit-sphere metric alone, so a general Lorentz
    transformation — a boost composed with rotations — rescales it by
    ``κ²`` with the very same ``κ``.  And, of course, we can use that
    rotational invariance to align the boost with any direction we
    like, so the result is general.

### [From the null cone to ``ℐ⁺``](@id from_cone_to_scri)

The construction so far was embedded in Minkowski space, with ``𝒩⁺``
the null cone of a single point and the fiber ``ℝ⁺`` acting by
dilation.  That picture is scaffolding.  A general asymptotically flat
spacetime has no preferred origin, its interior null congruences focus
and caustic, and it has no global null cone — yet the celestial sphere
survives intact, because the object we actually need is not the null
cone of a point but **future null infinity ``ℐ⁺`` itself**.

In the conformally compactified spacetime, ``ℐ⁺`` is a smooth null
hypersurface with topology ``ℝ × 𝕊²``.  Its null generators are the
integral curves of its degenerate direction, the retarded time ``u``
runs along each generator, and the space of generators is the celestial
sphere,

```math
π : ℐ⁺ → 𝕊² = ℐ⁺ / (\text{generators}).
```

This is the same bundle projection as before, but the change runs
deeper than swapping one fiber for another.  ``ℐ⁺`` is itself a null
cone: in the compactified spacetime its generators all converge on
future timelike infinity, so ``ℐ⁺`` is the *past* light cone of
``i⁺`` — an *absorption* cone — whereas the model we built was the
*future* light cone of a bulk point, an *emission* cone.  And the two
kinds of cone are ruled by *opposite* legs of the null pair.  For a
bulk vertex the comparison is elementary.  The future cone of the
origin has generator ``\{ λ\,(1, k̂) \mid λ > 0 \}`` at sky direction
``k̂``: its tangent is the outgoing leg ``ℓ ∝ (1, k̂)``.  The past
cone of the origin consists of points ``-s\,(1, p̂)``, which sit at
sky direction ``k̂ = -p̂`` and have future-pointing tangent ``(1, p̂)
= (1, -k̂) ∝ 𝐧``: the *ingoing* leg at that sky position.  This is
exactly the situation at ``ℐ⁺``: the outgoing ray with tangent ``ℓ ∝
(1, k̂)`` *ends* there, striking ``ℐ⁺`` transversally, while the
generator through the same point runs along the other leg, ``𝐧 ∝ (1,
-k̂)`` — the ``ñ = \sqrt{2}\,∂ᵤ`` of [the standard tetrad](@ref
"Tetrad"), with the retarded time ``u`` flowing along it.

The null-cone model is therefore faithful as a bundle over the sphere
— in Minkowski the generators of ``ℐ⁺`` are labeled by exactly the
directions ``k̂`` of the null rays through the origin — but *not*
leg-by-leg.  Following any generator of the origin's cone out to
infinity, the retarded time ``u = t - r = 0`` is constant, so the
entire cone arrives at the single cut ``u = 0`` of ``ℐ⁺``: the
dilation direction collapses at infinity, and the internal direction
that replaces it — the one the fiber ``ℝ`` of retarded time runs along
— is the *other* member of the null pair.  The two fibers are
abstractly the same one-dimensional abelian group (``ℝ⁺ ≅ ℝ`` under
``\log``), but they are attached along opposite null legs.  Passing
from the null cone to ``ℐ⁺`` swaps the degenerate direction of the
data surface from ``ℓ`` to ``𝐧``, while keeping the same celestial
sphere with the same labels.

!!! details "Isn't ``ℐ⁺`` also the light cone of ``i⁰``?"

    In compactified Minkowski the same surface can equally be
    described as the future light cone of spatial infinity ``i⁰`` —
    the cone of ``i⁰`` refocuses at ``i⁺``.  But the bulk analogy that
    correctly predicts the tangent structure is the past-cone one: at
    sky position ``k̂`` a bulk past cone is ruled by ``(1, -k̂) ∝
    𝐧``, matching ``ñ = \sqrt{2}\,∂ᵤ``, whereas a bulk future cone is
    ruled by ``ℓ``, which at ``ℐ⁺`` is transverse.  When in doubt,
    trust the tetrad: the generator is the leg the time coordinate
    flows along.

This swap is not bookkeeping; it decides which transformations act
*within* the surface.  Every BMS transformation maps generators of
``ℐ⁺`` to generators — the fibration is preserved — but no foliation
by cross-sections is preserved: there is no invariant family of cuts,
and that failure *is* the supertranslation freedom introduced below.
A transformation relating two frames adapted to ``ℐ⁺`` must therefore
fix the generator direction ``𝐧`` at each point, while nothing fixes
the transverse ``ℓ``.  In the language of [the Lorentz page](@ref
which_null_direction), tetrad transitions at ``ℐ⁺`` lie in the
stabilizer ``MAN_{𝐧}`` of the generator, so field components mix by
null rotations about ``𝐧`` — worked out on [the tetrad page](@ref
"BMS Action on the Tetrad") — whereas characteristic data on an
outgoing null cone, or radiation at ``ℐ⁻`` (the future light cone of
``i⁻``, ruled like an emission cone), selects ``N_{ℓ}`` instead.

Two features of ``ℐ⁺`` are worth isolating, because the rest of the page
rests on them:

* **The cross-section metric is only a conformal class, not a metric.**
    The conformal factor ``Ω`` used to reach ``ℐ⁺`` is fixed only up to
    ``Ω → ωΩ``, which rescales the induced sphere metric ``q → ω²q``.  So
    ``ℐ⁺`` hands us a *conformal class* of round metrics; the specific
    round ``dΩ²`` above is the representative singled out by a choice of
    inertial frame, not a canonical object.  This is exactly why the
    Lorentz action below is conformal rather than isometric, and why
    ``κ`` is a conformal factor.
* **``ℐ⁺`` is smooth by assumption.**  Completeness and regularity of
    the generators are part of the *definition* of asymptotic flatness,
    not something re-established per spacetime.  The caustics and the
    missing origin are features of the bulk, and the asymptotic symmetry
    group never looks at the bulk.

The fiber translations ``u ↦ u − α`` will turn out to be the
supertranslations, and the conformal maps of the base ``𝕊²`` the Lorentz
group; assembling the two is the business of the rest of this page.

## Supertranslations ``𝒮``

We stick with the approach implicit in the original papers, which is
to define our transformations as acting *passively*.  And in this
case, the supertranslations are defined as transformations of the time
coordinate.  (Again, we can refer to either the retarded time ``u`` or
the advanced time ``v``, so we will use the more neutral notation
``t``.)  So a supertranslation relates the time coordinate of a point
``p`` in one frame to the time coordinate of *the same point* ``p`` in
another frame as

```math
t'(p) = t(p) - c_α α(p),
```

where ``p`` refers to a point on ``ℐ``, and ``α`` is a real-valued
function that is independent of ``t``.  The ``c_α`` is just a
bookkeeping device to keep track of the differing conventions for the
sign in the formula above.

## BMS

We now have two groups: the (proper orthochronous) Lorentz group ``ℒ``
and the supertranslations ``𝒮``.  The supertranslation ``α ∈ 𝒮``
acts solely on the time coordinate as

```math
t'(t, 𝐤) = t - c_α α(𝐤).
```

The Lorentz transformation ``Λ ∈ ℒ`` acts on the spatial coordinates
as well as the time coordinate, the latter as

```math
t'(t, 𝐤) = κ(Λ, 𝐤) t.
```

We can combine these two transformations, applying ``α`` first and
then ``Λ``, to get the general transformation of the time coordinate:

```math
t'(t, 𝐤) = κ(Λ, 𝐤) [t - c_α α(𝐤)].
```

Now, we can investigate how repeated transformations compose.  Take
``α₁, α₂ ∈ 𝒮`` and ``Λ₁, Λ₂ ∈ ℒ``, and we apply ``α₁`` followed by
``Λ₁``, and then ``α₂`` followed by ``Λ₂``.  First, note that the null
ray transforms as

```math
𝐤' = Λ₁ 𝐤.
```

Also, using the definition of ``κ``, we can compute[^1]

```math
\begin{aligned}
κ(Λ₂Λ₁, 𝐤)
&= \frac{k⁰}{{(Λ₂Λ₁)⁰}ᵦkᵝ} \\
&= \frac{k⁰}{{(Λ₂)⁰}ᵧ {(Λ₁)ᵞ}ᵦkᵝ} \\
&= \frac{{(Λ₁)⁰}ᵦkᵝ}{{(Λ₂)⁰}ᵧ {(Λ₁)ᵞ}ᵦkᵝ} \frac{k⁰}{{(Λ₁)⁰}ᵦkᵝ} \\
&= κ(Λ₂, Λ₁𝐤)\,κ(Λ₁, 𝐤).
\end{aligned}
```

These allow us to write

```math
\begin{aligned}
t''(t', 𝐤')
&= κ(Λ₂, 𝐤') [t' - c_α α₂(𝐤')] \\
&= κ(Λ₂, Λ₁𝐤) \left\{κ(Λ₁, 𝐤) [t - c_α α₁(𝐤)] - c_α α₂(Λ₁𝐤)\right\} \\
&= κ(Λ₂Λ₁, 𝐤) \left\{t - c_α [α₁(𝐤) + α₂(Λ₁𝐤) / κ(Λ₁, 𝐤)]\right\} \\
\end{aligned}
```

That is, the result of a sequence of two transformations is another
transformation of the same form:

```math
t''(t, 𝐤) = κ(Λ, 𝐤) [t - c_α α(𝐤)],
```

where

```math
\begin{aligned}
Λ &= Λ₂ Λ₁, \\
α(𝐤) &= α₁(𝐤) + α₂(Λ₁ 𝐤)/κ(Λ₁, 𝐤).
\end{aligned}
```

Note that ``α₂`` is to be evaluated at the same geometric point on
``ℐ`` as ``α₁`` and ``κ`` in this expression, but is presumably
expressed with respect to the *transformed* frame — it is given as a
function of the null rays in those transformed components.

[^1]: The equality between the first and last lines of this equation
    is the defining feature of a "crossed homomorphism" — or more
    specifically a "1-cocycle" — of the Lorentz group with values in
    the multiplicative group of positive real numbers (the ``κ``
    function) [Brown_1982, nlab:crossed_homomorphism](@cite).  Though
    we explicitly compute it here, that equality is a general
    consequence of the definition of ``κ`` as the factor needed to
    restore the preferred section after a Lorentz transformation.

We can formalize this a little more directly by forming the (outer)
[semidirect
product](https://en.wikipedia.org/wiki/Semidirect_product#Outer_semidirect_product)
[Aluffi_2009; page 230](@cite) of these two groups:

```math
\text{BMS} = ℒ ⋉_φ 𝒮,
```

where the homomorphism ``φ: ℒ → \text{Aut}(𝒮)`` is given by

```math
φ(Λ)(α) = α/κ.
```

The identity element of this group is the pair ``(1, 0)``, and a
general element is

```math
(Λ, α) = (Λ, 0) ∘ (1, α).
```

The group operation is defined as

```math
(Λ₂, α₂) (Λ₁, α₁) = (Λ₂ Λ₁, α₁ + α₂ ∘ Λ₁ / κ₁).
```

It's easy to see that multiplication by ``1/κ = γ(1-v⃗⋅k̂)`` preserves
the defining properties of the supertranslations, so ``φ(Λ)`` is
indeed an automorphism.  And repeated multiplication by ``1/κᵢ`` is
consistent with the composition of Lorentz transformations, so ``φ``
is a homomorphism.

The groups ``ℒ`` and ``𝒮`` are isomorphic to subgroups of
``\text{BMS}`` in the obvious way: we can identify any Lorentz
transformation with the element ``(Λ, 0)``, and any supertranslation
with the element ``(1, α)``.  We will abuse notation slightly and
refer to the elements ``Λ`` and ``α`` themselves as elements of the
BMS group.  Note that, by this construction, ``𝒮`` is a [normal
subgroup](https://en.wikipedia.org/wiki/Normal_subgroup) of
``\text{BMS}``, but ``ℒ`` is not.  That is, for ``Λ ∈ ℒ`` and ``α ∈
𝒮``, the element ``Λ α Λ⁻¹`` is still an element of ``𝒮``, but ``α Λ
α⁻¹`` is not an element of ``ℒ``.  This fact will be useful later.

We have essentially *constructed* the BMS group, but
[Sachs_1962a](@citet) actually *derived* it from the asymptotic metric
conditions, and analyzed the group structure after the fact
[Sachs_1962b](@cite).  [Sachs_1962a](@citet) was the first to describe
the BMS group.  (He referred to it as the "Generalized Bondi-Metzner
group", or GBM group; later authors renamed it the Bondi-Metzner-Sachs
group to honor his contribution.)  His slightly later paper
[Sachs_1962b](@cite) was more specifically about the BMS group itself,
and proved some important properties, including:

  1. The supertranslations form an abelian normal subgroup ``N`` of
     the BMS group; the factor group is isomorphic to the Poincaré
     group.
  2. The translations form a normal four-dimensional subgroup of the
     proper BMS group.
  3. If ``N'`` is a four dimensional normal subgroup of the proper BMS
     group then ``N'`` is contained in the supertranslation group
     ``N``.
  4. The only normal four dimensional subgroup of the BMS group is the
     translation group.

In particular, the Poincaré group — which Sachs refers to by its
longer name: the inhomogeneous (proper) orthochronous Lorentz group —
is a subgroup of the BMS group, though it is not a normal subgroup
because conjugating a spatial translation by a Lorentz transformation
results in a general supertranslation.

!!! note "Which way around is κ?"

    It is easy to confuse ``κ`` with ``1/κ`` in these formulas — not
    least because different references actually define it both ways.  A
    discerning consistency check is *Poincaré closure*.  That is, we
    perform a transformation that comes from the Poincaré subgroup of
    BMS and check that the BMS result is the same as the simpler
    Poincaré result.  Specifically, the time translation ``α = δt``
    corresponds to a spacetime translation by the 4-vector
    ``\boldsymbol{δ} = (δt, 0, 0, 0)``.  We can readily compute
    ``Λ^{-1} \boldsymbol{δ} Λ = γ\,δt\,(1; v⃗)`` for a pure boost
    with velocity v⃗, which is equivalent to the BMS
    supertranslation ``α' = γ\,δt\,(1 - v⃗⋅k̂) = δt/κ``.

    Our group multiplication law is consistent, because we take
    ``Λ₂=Λ⁻¹``, ``Λ₁=Λ``, ``α₁=0``, and ``α₂=δt`` to find that
    ```math
    (Λ₂, α₂) (Λ₁, α₁) = Λ⁻¹ ∘ δt ∘ Λ = (1, α ∘ Λ / κ) = (1, δt / κ).
    ```
    The last equality is because ``α`` is constant over the sphere,
    so evaluating in the direction ``Λ𝐤`` is the same as evaluating
    in the direction ``𝐤``, and that value is ``δt``.  That factor
    of  ``1/κ`` is given by our composition law, and is exactly the
    same as the result of the Poincaré transformation.

## [Representations of the supertranslation](@id bms_representations)

Strictly speaking, there are *two* BMS groups: ``\text{BMS}⁺`` acting
on ``ℐ⁺`` and ``\text{BMS}⁻`` acting on ``ℐ⁻``.  They are isomorphic —
canonically so, under the antipodal matching of [Future and past null
infinity](@ref scri_pm_conventions) [Strominger_2014,
Strominger_2017](@cite) — so nothing about the group structure
distinguishes them; what differs is purely the *representation* of
elements.  The supertranslation ``α`` is stored as mode weights of a
function of the labels ``k̂`` on the celestial sphere, and the two
ends of null infinity are labeled antipodally: at ``ℐ⁺`` a generator
is labeled by its outgoing propagation direction, while at ``ℐ⁻`` it
is labeled by the direction in which an observer *sees* it — opposite
to the propagation.  The same abstract supertranslation therefore has
two mode representations, related by composition with the antipodal
map ``A: k̂ ↦ -k̂``:

```math
α ↦ α ∘ A,
\qquad
(α ∘ A)_{ℓ,m} = (-1)^ℓ\, α_{ℓ,m},
```

since ``Y_{ℓ,m}(-k̂) = (-1)^ℓ\, Y_{ℓ,m}(k̂)``.  The Lorentz part needs
no such choice — a rotor is representation-neutral — but its
*realization* as a map of labels does depend on the labeling: on the
antipodal labels, ``Λ`` acts through the conjugated map ``A ∘ Λ ∘ A``,
which flips the boost part (``v⃗ ↦ -v⃗``) while leaving rotations
alone — exactly why the conformal factor and the direction map involve
the ``ℐ`` sign.

The [`BMS`](@ref) type therefore retains that choice as data — the
`I` type parameter — declaring which labeling the stored modes refer
to, alongside the analogous supertranslation-sign convention `A` (the
sign with which ``α`` enters the time law ``t' = κ(t - c_αα)``).  Both
default to `+1`: future null infinity, and the plus sign.  The
conversion constructor `BMS(g; c_α, ℐ)` (see [`BMS`](@ref)) converts
an element between representations — exactly, since both conversions
are pure sign flips.

Because the twist in the twisted group law,

```math
α(𝐤) = α₁(𝐤) + α₂(Λ₁𝐤)/κ(Λ₁, 𝐤),
```

is built from the label map ``Λ₁𝐤`` and the conformal factor — both
realized on the *labeled* sphere — the coordinates of a composite
element depend on the representation even though the abstract group
does not.  To make that precise, write ``∘⁺`` and ``∘⁻`` for the two
composition laws: each is the formula above, with the label map and
conformal factor computed from the corresponding section ``𝐤 = (1,
±k̂)``.  Note that ``κ`` itself needs no such decoration: as defined
[above](@ref "Lorentz transformations ``ℒ`` and the conformal factor
``κ``"), it is a single function of ``Λ`` and a null vector — indeed,
being scale-invariant in ``𝐤``, a function of the null *ray* — and
the two labelings merely feed it antipodal rays, ``κ(Λ, (1, -k̂))`` on
``ℐ⁻`` whereas it is ``κ(Λ, (1, +k̂))`` on ``ℐ⁺``.  The same is true
of the label map.  Now promote the relabeling of ``α`` to a map on
whole elements,

```math
\begin{aligned}
P \colon \text{BMS}⁺ &→ \text{BMS}⁻, \\
(Λ, α) &↦ (Λ, α ∘ A),
\end{aligned}
```

which sends the coordinates of a ``\text{BMS}⁺`` element to the
coordinates of the matched ``\text{BMS}⁻`` element — and, being an
involution, back again.  Then the two composition laws are related by

```math
P(g₂) ∘⁻ P(g₁) = P(g₂ ∘⁺ g₁).
```

That is, composing the matched elements at ``ℐ⁻`` yields the matched
composite — precisely the statement that ``P`` is a group isomorphism
``\text{BMS}⁺ → \text{BMS}⁻``.  Pointwise, the matched element acts
antipodally: if ``g`` maps ``(t, k̂) ↦ (t', k̂')`` on ``ℐ⁺``, then
``P(g)`` maps ``(t, -k̂) ↦ (t', -k̂')`` on ``ℐ⁻``.

In the code, all this bookkeeping is handled automatically.  Operations
that combine two elements ([`compose`](@ref Scri.compose), `==`,
`isapprox`) accept any mix of conventions, accounting for the
differences exactly; each convention the two inputs share is kept for
the result, and each convention on which they disagree is resolved to
the default `+1`.  Single-element operations — the action functor,
[`conformal_factor`](@ref Scri.conformal_factor), `inv` — use the
element's own conventions, and the [`transform!`](@ref) bridge
re-represents the element to match the null infinity declared by its
[`DataComponents`](@ref Scri.DataComponents) argument.

Note that the ``ℐ`` relevant to `DataComponents` is conceptually
distinct from the ``ℐ`` relevant to the BMS element.  It may be
quite reasonable, for example, to always consider the effect a BMS
element has on future null infinity, but still be concerned with the
transformations of data on past null infinity.  This is the reason why
both `BMS` and `DataComponents` each have their own `ℐ` parameter.

## Decomposition of BMS

Essentially by our definition, we have already decomposed
``\text{BMS}`` into the supertranslations and the Lorentz
transformations.  However, it can also be useful to further decompose
these parts.  Once we have chosen a frame ``𝐭, 𝐱, 𝐲, 𝐳`` to work
with, we can conventionally decompose any element of the BMS group as
shown in the diagram below.

```@raw html
<div class="composition-diagram" style="text-align: center; margin: 1.5em 0;">
<!--
  Node layout — to move a node, change its translate(cx, cy).
  Row y-centers:  0=30, 1=107, 2=214, 3=335
  Edges use those same centers; boxes & text are local to each <g>.
-->
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1080 390"
     width="100%" style="max-width: 740px;">
  <style>
    .bms-tree text {
      fill: currentColor;
      text-anchor: middle;
    }
    .bms-tree .box {
      fill: var(--bg, #888);
      stroke: currentColor;
      stroke-width: 1.4;
    }
    .bms-tree .edge {
      stroke: currentColor;
      stroke-width: 1.2;
    }
  </style>

  <g class="bms-tree">

    <!-- ===== Edges (drawn first, behind nodes) ===== -->

    <!-- BMS → Supertranslations, Lorentz -->
    <line class="edge" x1="530" y1="30" x2="260" y2="107"/>
    <line class="edge" x1="530" y1="30" x2="800" y2="107"/>

    <!-- Supertranslations → Proper, Spacetime -->
    <line class="edge" x1="260" y1="107" x2="100" y2="214"/>
    <line class="edge" x1="260" y1="107" x2="370" y2="214"/>

    <!-- Lorentz → Rotations, Boosts -->
    <line class="edge" x1="800" y1="107" x2="680" y2="214"/>
    <line class="edge" x1="800" y1="107" x2="920" y2="214"/>

    <!-- Spacetime → Time, Space -->
    <line class="edge" x1="370" y1="214" x2="300" y2="335"/>
    <line class="edge" x1="370" y1="214" x2="440" y2="335"/>

    <!-- Rotations → Spherical, Spin angle -->
    <line class="edge" x1="680" y1="214" x2="630" y2="335"/>
    <line class="edge" x1="680" y1="214" x2="750" y2="335"/>

    <!-- Boosts → 1-D boosts, Null rotations -->
    <line class="edge" x1="920" y1="214" x2="880" y2="335"/>
    <line class="edge" x1="920" y1="214" x2="1000" y2="335"/>

    <!-- ===== Nodes ===== -->

    <!-- Level 0: BMS -->
    <g transform="translate(530, 30)">
      <rect class="box" x="-40" y="-17" width="80" height="34" rx="6"/>
      <text y="6" class="katex"><tspan class="mord">BMS</tspan></text>
    </g>

    <!-- Level 1: Supertranslations 𝒮 -->
    <g transform="translate(260, 107)">
      <rect class="box" x="-100" y="-25" width="200" height="50" rx="6"/>
      <text y="-5">Supertranslations <tspan class="katex"><tspan class="mord mathscr">S</tspan></tspan></text>
      <text y="14" class="katex"><tspan class="mord mathnormal">α</tspan></text>
    </g>

    <!-- Level 1: Lorentz ℒ -->
    <g transform="translate(800, 107)">
      <rect class="box" x="-90" y="-25" width="180" height="50" rx="6"/>
      <text y="-5">Lorentz <tspan class="katex"><tspan class="mord mathscr">L</tspan></tspan></text>
      <text y="14" class="katex"><tspan class="mord">Λ</tspan></text>
    </g>

    <!-- Level 2: Proper supertranslations -->
    <g transform="translate(100, 214)">
      <rect class="box" x="-95" y="-25" width="190" height="50" rx="6"/>
      <text y="-5">Proper</text>
      <text y="14">supertranslations</text>
    </g>

    <!-- Level 2: Spacetime translations -->
    <g transform="translate(370, 214)">
      <rect class="box" x="-85" y="-25" width="170" height="50" rx="6"/>
      <text y="-5">Spacetime</text>
      <text y="14">translations</text>
    </g>

    <!-- Level 2: Rotations Spin(3) -->
    <g transform="translate(680, 214)">
      <rect class="box" x="-65" y="-25" width="130" height="50" rx="6"/>
      <text y="-5">Rotations</text>
      <text y="14" class="katex"><tspan class="mord mathnormal">R</tspan></text>
    </g>

    <!-- Level 2: Boosts -->
    <g transform="translate(920, 214)">
      <rect class="box" x="-50" y="-25" width="100" height="50" rx="6"/>
      <text y="-5">Boosts</text>
      <text y="14" class="katex"><tspan class="mord mathnormal">v⃗</tspan></text>
    </g>

    <!-- Level 3: Time translations δt -->
    <g transform="translate(300, 335)">
      <rect class="box" x="-60" y="-31" width="120" height="62" rx="6"/>
      <text y="-13">Time</text>
      <text y="4">translations</text>
      <text y="21" class="katex"><tspan class="mord mathnormal">δt</tspan></text>
    </g>

    <!-- Level 3: Space translations δx⃗ -->
    <g transform="translate(440, 335)">
      <rect class="box" x="-65" y="-31" width="130" height="62" rx="6"/>
      <text y="-13">Space</text>
      <text y="4">translations</text>
      <text y="21" class="katex"><tspan class="mord mathnormal">δx⃗</tspan></text>
    </g>

    <!-- Level 3: Spherical coordinates (θ, φ) -->
    <g transform="translate(630, 335)">
      <rect class="box" x="-55" y="-31" width="110" height="62" rx="6"/>
      <text y="-13">Spherical</text>
      <text y="4">coordinates</text>
      <text y="21" class="katex"><tspan class="mord mathnormal">(<tspan class="mord mathnormal">θ</tspan>, <tspan class="mord mathnormal">φ</tspan>)</tspan></text>
    </g>

    <!-- Level 3: Spin angle ψ -->
    <g transform="translate(750, 335)">
      <rect class="box" x="-50" y="-25" width="100" height="50" rx="6"/>
      <text y="-5">Spin angle</text>
      <text y="14" class="katex"><tspan class="mord mathnormal">ψ</tspan></text>
    </g>

    <!-- Level 3: 1-D boosts η -->
    <g transform="translate(880, 335)">
      <rect class="box" x="-55" y="-25" width="110" height="50" rx="6"/>
      <text y="-5">1-D boosts</text>
      <text y="14" class="katex"><tspan class="mord mathnormal">η</tspan></text>
    </g>

    <!-- Level 3: Null rotations -->
    <g transform="translate(1000, 335)">
      <rect class="box" x="-50" y="-31" width="100" height="62" rx="6"/>
      <text y="-13">Null</text>
      <text y="4">rotations</text>
      <text y="21" class="katex"><tspan class="mord mathnormal">ζ</tspan></text>
    </g>
  </g>
</svg>
<script>
(function() {
  var bg = getComputedStyle(document.documentElement).backgroundColor;
  document.querySelector('.bms-tree').style.setProperty('--bg', bg);
})();
</script>
</div>
```

The supertranslations decompose naturally into "proper"
supertranslations and ordinary spacetime translations — the latter of
which further decompose into time translations ``δt`` and space
translations ``δx⃗``.  In general, the supertranslations are functions
on the sphere, which we decompose into spherical harmonics.  If the
supertranslation is a pure spacetime translation, then the only
nonzero spherical harmonic modes are the ``ℓ=0`` mode for time
translations and the ``ℓ=1`` modes for space translations.
Equivalently, we can think of such translations in terms of their
physical parameters ``δt`` and ``δx⃗``.

The Lorentz transformations decompose naturally into rotations ``R``
and boosts ``v⃗``.  It is already convenient to represent the rotation
as a quaternion, which can be constructed in numerous ways, including
the generator of the rotation, the axis and angle of the rotation, the
Euler angles, or the spherical coordinates and spin angle (which is
essentially a different version of Euler angles).  Finally, a general
boost may be decomposed into a 1-D boost of rapidity ``η`` in the
``(θ, φ)`` direction, and a null rotation of parameter ``ζ`` about
that direction.  This is not often an intuitive decomposition, but it
is important in discussions of boost weight.
