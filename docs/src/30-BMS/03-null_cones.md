# Null cones

!!! note "TO DO"

    Introduction.  The classical stereographic and complex-analytic
    approach of Penrose is recovered from the geometry, rather than
    adopted as an artificial device.

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
``xᴬ`` on ``𝕊²`` — the usual ``(θ, ϕ)``, say — the section is

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

## The plane of null rotations is the conformal plane

!!! note "TO DO"

    Fill this in

## Lorentz transformations ``ℒ`` and the conformal factor ``κ``

A Lorentz transformation ``Λ`` preserves the null cone ``𝒩⁺`` — null
vectors map to null vectors — but it does *not* preserve the section:
the image ``Λσ`` of a unit-time-component null vector generally has
some other time component, so it no longer satisfies ``k⁰ = 1``.  The
preferred section is restored by a simple rescaling, and it turns out
that the factor that does it is the conformal factor by which the
metric transforms.  Define

```math
κ(Λ, 𝐤) = \frac{k⁰}{{Λ⁰}ᵦkᵝ} = \frac{1}{γ(1-v⃗⋅k̂)},
```

the ratio of the old time component to the new one.  The second form
is the familiar Doppler factor of a pure boost with velocity ``v⃗``,
recovered explicitly on the [aberration page](@ref "Aberration of
Gravitational Waves").  A Lorentz transformation acts on the section
by

```math
σ' = κ\, (Λ ∘ σ),
```

the factor of ``κ`` undoing the change in time component, so that
``σ'`` again has ``k'^{0'} = 1``.

!!! note "The conformal factor and KAN"

    We saw [previously](@ref "Iwasawa's ``KAN`` decomposition") that
    the KAN decomposition gives us a very useful way of understanding
    Lorentz transformations.  The conformal factor ``κ`` will also be
    extremely important in the BMS group, so it will be useful to
    relate the two.  Every KAN decomposition must be taken *with
    respect to* some null direction.  If we choose ``𝐤``, then we saw
    that the effect of the ``AN`` factors is to rescale the null
    vector by a factor ``e^{φₐ}``.  The rotation ``K`` does not change
    the time component, so the conformal factor is simply

    ```math
    κ = e^{-φₐ}
    ```

    Where we emphasize that this ``φₐ`` is the one that results from
    decomposition with respect to ``𝐤``.  Alternatively, we can note
    that if ``R`` is any rotation that takes the standard reference
    vector ``\boldsymbol{ℓ} = (𝐭 + 𝐳) / \sqrt{2}`` to ``𝐤``, then
    ``κ(Λ, 𝐤) = κ(ΛR, \boldsymbol{ℓ})``, so ``φₐ`` also results from
    the KAN decomposition of ``ΛR`` with respect to
    ``\boldsymbol{ℓ}``.

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

We can simplify this expansion using three facts.  First, ``Λ``
preserves the metric:

```math
ηᵦᵧ {Λᵝ}ᵤ {Λᵞ}ᵥ = ηᵤᵥ.
```

Second, the section is null:

```math
ηᵦᵧ σᵝ σᵞ = 0.
```

And third, differentiating that null condition shows that the
section's derivatives are orthogonal to the section itself:

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
    dϕ²`` therefore rescales as ``{dΩ'}² = κ²\, dΩ²``.  Rotations
    leave the unit-sphere metric alone, so a general Lorentz
    transformation — a boost composed with rotations — rescales it by
    ``κ²`` with the very same ``κ``.  And, of course, we can use that
    rotational invariance to align the boost with any direction we
    like, so the result is general.

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

## [From the null cone to ``ℐ⁺``](@id from_cone_to_scri)

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
abstractly the same one-dimensional abelian group (``ℝ⁺ \cong ℝ``
under ``\log``), but they are attached along opposite null legs.
Passing from the null cone to ``ℐ⁺`` swaps the degenerate direction of
the data surface from ``ℓ`` to ``𝐧``, while keeping the same
celestial sphere with the same labels.

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
\cong
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
