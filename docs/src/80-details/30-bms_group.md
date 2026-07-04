# [The BMS Group](@id bms_group)

## Bondi-Sachs coordinates, metric, and gauge

## The null cone and the celestial sphere

* Future null cone and celestial sphere
  * ``𝒩⁺`` with metric asymptotically Minkowski ``ηᵦᵧ``
  * celestial sphere ``𝕊²`` is the space of future null *rays* ``𝒩⁺/ℝ⁺``
    * choose some ``t`` vector
    * each null vector decomposes as a term proportional to ``t`` plus
      a spatial vector
    * each null vector then gives rise to a unique spatial unit vector
      ``n̂`` by normalizing the spatial part
    * consider ``𝕊²`` to be represented by these spatial unit vectors
    * also use coordinates ``xᴬ``, which may be ``(θ, φ)``
  * section ``σ:𝕊²→𝒩⁺`` via ``σ:n̂↦ℓᵝ=(1,n̂)``.  That is,
    ``σᵝ(xᴬ) = (1, n̂(xᴬ))``
  * (note that ``ℓᵝ = -ηᵝᵞ(dt-dr)ᵧ \sim -ηᵝᵞ(du)ᵧ``)
  * induced metric on ``𝕊²`` is the pullback

    ```math
    dΩ² = σ^* η = ηᵦᵧ \frac{∂σᵝ}{∂xᴬ} \frac{∂σᵞ}{∂xᴮ} dxᴬ dxᴮ
    ```

    or

    ```math
    g_{AB} = ηᵦᵧ \frac{∂σᵝ}{∂xᴬ} \frac{∂σᵞ}{∂xᴮ}.
    ```

## Lorentz transformations ``ℒ`` and the conformal factor ``κ``

* Lorentz transformations preserve ``𝒩⁺`` but not the section with
  ``ℓ⁰=1``.
* Define ``κ`` in terms of time component

  ```math
  κ(Λ, 𝐤) = \frac{k⁰}{{Λ⁰}ᵦkᵝ} = \frac{1}{γ(1-v⃗⋅n̂)}.
  ```

* Lorentz transformation ``Λ`` induces a transformation of the section
  as ``σ' = κ (Λ ∘ σ)`` — where we have to rescale by ``κ`` to get
  back to the preferred section with ``ℓ'^{0'}=1``.
* Differentiating ``{σ'}ᵝ = κ {Λᵝ}ᵧσᵞ`` by ``xᴬ`` gives us two terms,
  the first differentiating ``κ`` and the second differentiating
  ``σ``:

  ```math
  \frac{∂{σ'}ᵝ}{∂xᴬ} = \frac{∂κ}{∂xᴬ} {Λᵝ}ᵧσᵞ + κ {Λᵝ}ᵧ \frac{∂σᵞ}{∂xᴬ}.
  ```

  That derivative is then contracted twice with the Minkowski metric
  to give us the new metric in the primed frame:

  ```math
  g'_{AB} = ηᵦᵧ \frac{∂{σ'}ᵝ}{∂xᴬ} \frac{∂{σ'}ᵞ}{∂xᴮ}.
  ```

  We will expand this expression, but we need a few preliminary
  results first.  First, note that in all resulting terms, we have the
  ``Λ``s contracted with the metric, which they preserve:

  ```math
  ηᵦᵧ {Λᵝ}ᵤ {Λᵞ}ᵥ = ηᵤᵥ.
  ```

  Next, we will need the fact that the section is null, which means
  that

  ```math
  ηᵦᵧ σᵝ σᵞ = 0.
  ```

  Finally, we can differentiate that expression to find that

  ```math
  ηᵦᵧ σᵝ \frac{∂σᵞ}{∂xᴬ} = 0.
  ```

  We can now use the expression for the derivative of ``σ'`` to expand
  the expression for ``g'_{AB}``, and use these results to simplify,
  then find

  ```math
  g'_{AB} = κ² ηᵦᵧ \frac{∂{σ}ᵝ}{∂xᴬ} \frac{∂{σ}ᵞ}{∂xᴮ} = κ² g_{AB}.
  ```

  That is, ``{dΩ'}² = κ² dΩ²``.  This is the key result: Lorentz
  transformations transform the unit sphere metric on the celestial
  sphere by a conformal factor ``κ²``.
* Combine that with the fact that Bondi gauge requires the angular
  part of the metric to be *asymptotically* ``r² dΩ²``, and if the
  transformation is an isometry we must have the *asymptotic* relation
  ``r² dΩ² \sim {r'}² {dΩ'}²``, and we find that ``r \sim κ r'``,
  which is Sachs's *definition* of ``κ``.

!!! info "To do"

    Here's a more mechanical and unenlightening derivation, though it
    may be more familiar, so it could be useful to have both.
    * Conformal factor under boost:
      - aberration formula ``\cos θ' = (\cos θ - β) / (1 - β \cos θ)``
      - differentiate to find ``\sin θ'\, dθ' = κ² \sin θ\, dθ``
      - use ``\sin² θ' = 1-\cos² θ'`` to find ``\sin θ' = κ \sin θ``
      - also have ``dθ' = κ dθ``
      - combine to show that a boost along ``z`` transforms the unit
        sphere metric as ``{dΩ'}² = κ² dΩ²``.
      - rotations preserve the unit sphere metric, so arbitrary
        Lorentz transformations transform the unit sphere metric as
        ``{dΩ'}² = κ² dΩ²``, with ``κ`` as we defined it.

!!! info "To do"

    Show that the conformal factor of a product of Lorentz
    transformations is the product of the conformal factors:
    ``κ(Λ₂ Λ₁) = κ(Λ₂) κ(Λ₁)``.  This is a consequence of the
    group structure, but it is not *entirely* trivial to show.
    Specifically, it's technically a *crossed* homomorphism,
    because the second factor is evaluated at the transformed
    point; if we ignore the evaluation point, then it's just a
    homomorphism.

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
t'(p) = t(p) - εᵅ α(p),
```

where ``p`` refers to a point on ``ℐ``, and ``α`` is a real-valued
function that is independent of ``t``.  The ``εᵅ`` is just a
bookkeeping device to keep track of the differing conventions for the
sign in the formula above.

## BMS

We now have two groups: the (proper orthochronous) Lorentz group ``ℒ``
and the supertranslations ``𝒮``.  The supertranslation ``α ∈ 𝒮``
acts solely on the time coordinate as

```math
t'(t, 𝐤) = t - εᵅ α(𝐤).
```

The Lorentz transformation ``Λ ∈ ℒ`` acts on the spatial coordinates
as well as the time coordinate, the latter as

```math
t'(t, 𝐤) = κ(Λ, 𝐤) t.
```

We can combine these two transformations, applying ``α`` first and
then ``Λ``, to get the general transformation of the time coordinate:

```math
t'(t, 𝐤) = κ(Λ, 𝐤) [t - εᵅ α(𝐤)].
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
&= κ(Λ₂, 𝐤') [t' - εᵅ α₂(𝐤')] \\
&= κ(Λ₂, Λ₁𝐤) \left\{κ(Λ₁, 𝐤) [t - εᵅ α₁(𝐤)] - εᵅ α₂(Λ₁𝐤)\right\} \\
&= κ(Λ₂Λ₁, 𝐤) \left\{t - εᵅ [α₁(𝐤) + α₂(Λ₁𝐤) / κ(Λ₁, 𝐤)]\right\} \\
\end{aligned}
```

That is, the result of a sequence of two transformations is another
transformation of the same form:

```math
t''(t, 𝐤) = κ(Λ, 𝐤) [t - εᵅ α(𝐤)],
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

It's easy to see that multiplication by ``1/κ = γ(1-v⃗⋅n̂)`` preserves
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
    supertranslation ``α' = γ\,δt\,(1 - v⃗⋅n̂) = δt/κ``.

    Our group multiplication law is consistent, because we take
    ``Λ₂=Λ⁻¹``, ``Λ₁=Λ``, ``α₁=0``, and ``α₂=δt`` to find that
    ```math
    (Λ₂, α₂) (Λ₁, α₁) = Λ⁻¹ ∘ δt ∘ Λ = (1, α ∘ Λ / κ) = (1, δt / κ).
    ```
    The last equality is because ``α`` is constant, so evaluating
    at ``Λ𝐤`` is the same as evaluating at ``𝐤``, and that value is
    ``δt``.  That factor of  ``1/κ`` is given by our composition
    law, and is exactly the same as the result of the Poincaré
    transformation.

## [Representations of the supertranslation](@id bms_representations)

Strictly speaking, there are *two* BMS groups: ``\text{BMS}⁺`` acting
on ``ℐ⁺`` and ``\text{BMS}⁻`` acting on ``ℐ⁻``.  They are isomorphic —
canonically so, under the antipodal matching of [Future and past null
infinity](@ref scri_pm_conventions) [Strominger_2014,
Strominger_2017](@cite) — so nothing about the group structure
distinguishes them; what differs is purely the *representation* of
elements.  The supertranslation ``α`` is stored as mode weights of a
function of the labels ``n̂`` on the celestial sphere, and the two
ends of null infinity are labeled antipodally: at ``ℐ⁺`` a generator
is labeled by its outgoing propagation direction, while at ``ℐ⁻`` it
is labeled by the direction in which an observer *sees* it — opposite
to the propagation.  The same abstract supertranslation therefore has
two mode representations, related by composition with the antipodal
map ``A: n̂ ↦ -n̂``:

```math
α ↦ α ∘ A,
\qquad
(α ∘ A)_{ℓ,m} = (-1)^ℓ\, α_{ℓ,m},
```

since ``Y_{ℓ,m}(-n̂) = (-1)^ℓ\, Y_{ℓ,m}(n̂)``.  The Lorentz part needs
no such choice — a rotor is representation-neutral — but its
*realization* as a map of labels does depend on the labeling: on the
antipodal labels, ``Λ`` acts through the conjugated map ``A ∘ Λ ∘ A``,
which flips the boost part (``v⃗ ↦ -v⃗``) while leaving rotations
alone — exactly why the conformal factor and the direction map involve
the ``ε^ℐ`` sign.

The [`BMS`](@ref) type therefore retains that choice as data — the
`Eᴵ` type parameter — declaring which labeling the stored modes refer
to, alongside the analogous supertranslation-sign convention `Eᵅ` (the
sign with which ``α`` enters the time law ``t' = κ(t - εᵅα)``).  Both
default to `+1`: future null infinity, and the plus sign.  The
conversion constructor `BMS(g; εᵅ, εᴵ)` (see [`BMS`](@ref)) converts
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
±n̂)``.  Note that ``κ`` itself needs no such decoration: as defined
[above](@ref "Lorentz transformations ``ℒ`` and the conformal factor
``κ``"), it is a single function of ``Λ`` and a null vector — indeed,
being scale-invariant in ``𝐤``, a function of the null *ray* — and
the two labelings merely feed it antipodal rays, ``κ(Λ, (1, -n̂))`` on
``ℐ⁻`` whereas it is ``κ(Λ, (1, +n̂))`` on ``ℐ⁺``.  The same is true
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
antipodally: if ``g`` maps ``(t, n̂) ↦ (t', n̂')`` on ``ℐ⁺``, then
``P(g)`` maps ``(t, -n̂) ↦ (t', -n̂')`` on ``ℐ⁻``.

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

Note that the ``ε^ℐ`` relevant to `DataComponents` is conceptually
distinct from the ``ε^ℐ`` relevant to the BMS element.  It may be
quite reasonable, for example, to always consider the effect a BMS
element has on future null infinity, but still be concerned with the
transformations of data on past null infinity.  This is the reason why
both `BMS` and `DataComponents` each have their own `εᴵ` parameter.

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
