# The BMS Group

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

## Lorentz transformations ``ℒ`` and the conformal factor ``K``

* Lorentz transformations preserve ``𝒩⁺`` but not the section with
  ``ℓ⁰=1``.
* Define ``K`` in terms of time component ``{Λ⁰}_μ ℓ^μ = γ(1-v⃗⋅n̂) =
  1/K``.
* Lorentz transformation ``Λ`` induces a transformation of the section
  as ``σ' = K (Λ ∘ σ)`` — where we have to rescale by ``K`` to get
  back to the preferred section with ``ℓ'^{0'}=1``.
* Differentiating ``{σ'}ᵝ = K {Λᵝ}ᵧσᵞ`` by ``xᴬ`` gives us two terms,
  the first differentiating ``K`` and the second differentiating
  ``σ``:

  ```math
  \frac{∂{σ'}ᵝ}{∂xᴬ} = \frac{∂K}{∂xᴬ} {Λᵝ}ᵧσᵞ + K {Λᵝ}ᵧ \frac{∂σᵞ}{∂xᴬ}.
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
  g'_{AB} = K² ηᵦᵧ \frac{∂{σ}ᵝ}{∂xᴬ} \frac{∂{σ}ᵞ}{∂xᴮ} = K² g_{AB}.
  ```

  That is, ``{dΩ'}² = K² dΩ²``.  This is the key result: Lorentz
  transformations transform the unit sphere metric on the celestial
  sphere by a conformal factor ``K²``.
* Combine that with the fact that Bondi gauge requires the angular
  part of the metric to be *asymptotically* ``r² dΩ²``, and if the
  transformation is an isometry we must have the *asymptotic* relation
  ``r² dΩ² \sim {r'}² {dΩ'}²``, and we find that ``r \sim K r'``,
  which is Sachs's *definition* of ``K``.

!!! info "To do"

    Here's a more mechanical and unenlightening derivation, though it
    may be more familiar, so it could be useful to have both.
    * Conformal factor under boost:
      - aberration formula ``\cos θ' = (\cos θ - β) / (1 - β \cos θ)``
      - differentiate to find ``\sin θ'\, dθ' = K² \sin θ\, dθ``
      - use ``\sin² θ' = 1-\cos² θ'`` to find ``\sin θ' = K \sin θ``
      - also have ``dθ' = K dθ``
      - combine to show that a boost along ``z`` transforms the unit
        sphere metric as ``{dΩ'}² = K² dΩ²``.
      - rotations preserve the unit sphere metric, so arbitrary
        Lorentz transformations transform the unit sphere metric as
        ``{dΩ'}² = K² dΩ²``, with ``K`` as we defined it.

!!! info "To do"

    Show that the conformal factor of a product of Lorentz
    transformations is the product of the conformal factors:
    ``K(Λ₂ Λ₁) = K(Λ₂) K(Λ₁)``.  This is a consequence of the
    group structure, but it is not *entirely* trivial to show.

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
acts solely on the time coordinate as ``t' = t - εᵅ α``.  The Lorentz
transformation ``Λ ∈ ℒ`` acts on the spatial coordinates as well as
the time coordinate, the latter as ``t' = Kt``.  We can combine these
two transformations, applying ``α`` first and then ``Λ``, to get the
general transformation of the time coordinate:

```math
t' = K(t - εᵅ α),
```

where we implicitly evaluate all fields at the same point on ``ℐ``.
To be a little more explicit, we can treat these quantities as
functions on null rays: ``α(λ𝐤) = α(𝐤)`` for any scalar ``λ > 0``
and null vector ``𝐤``.  Then the transformation law is

```math
\begin{aligned}
t'(𝐤) &= K(𝐤)(t(𝐤) - εᵅ α(𝐤)), \\
𝐤' &= Λ 𝐤.
\end{aligned}
```

Recall that ``K(𝐤) = {Λ⁰}ᵦkᵝ/k⁰``.  Now, we can investigate how
repeated transformations compose.  Take ``α₁, α₂ ∈ 𝒮`` and ``Λ₁, Λ₂ ∈
ℒ``.  We have

```math
\begin{aligned}
t'(𝐤) &= K₂(𝐤') \left\{K₁(𝐤)[t(𝐤) - εᵅ α₁(𝐤)] - εᵅ α₂(𝐤')\right\} \\
      &= K₂(𝐤') K₁(𝐤) t(𝐤) - K₂(𝐤') εᵅ α₁(𝐤) - K₂(𝐤') εᵅ α₂(𝐤') \\
      &= [K₂(𝐤') K₁(𝐤)] \left\{t(𝐤) - εᵅ \left[α₁(𝐤) + α₂(𝐤')/K₁(𝐤)\right] \right\}.
\end{aligned}
```

That is, the result of a sequence of two transformations is another
transformation of the same form, with

```math
\begin{aligned}
Λ &= Λ₂ Λ₁, \\
α(𝐤) &= α₁(𝐤) + α₂(Λ₁ 𝐤)/K₁(𝐤).
\end{aligned}
```

Note that ``α₂`` is to be evaluated at the same geometric point on
``ℐ`` as ``α₁`` and ``K₁``, but is presumably expressed in the
*transformed* coordinates.

We can formalize this a little more directly by forming the (outer)
[semidirect
product](https://en.wikipedia.org/wiki/Semidirect_product#Outer_semidirect_product)
[Aluffi_2009; page 230](@cite) of these two groups:

```math
\text{BMS} = ℒ ⋉_φ 𝒮,
```

where the homomorphism ``φ: ℒ → \text{Aut}(𝒮)`` is given by

```math
φ(Λ)(α) = α/K.
```

The identity element of this group is the pair ``(1, 0)``, and a
general element is

```math
(Λ, α) = (Λ, 0) ∘ (1, α).
```

The group operation is defined as

```math
(Λ₂, α₂) (Λ₁, α₁) = (Λ₂ Λ₁, α₁ + α₂/K₁).
```

It's easy to see that multiplication by ``1/K = γ(1-v⃗⋅n̂)`` preserves
the defining properties of the supertranslations, so ``φ(Λ)`` is
indeed an automorphism.  And repeated multiplication by ``1/Kᵢ`` is
consistent with the composition of Lorentz transformations, so ``φ``
is a homomorphism.

The groups ``ℒ`` and ``𝒮`` are isomorphic to subgroups of
``\text{BMS}`` in the obvious way: we can identify any Lorentz
transformation with the element ``(Λ, 0)``, and any supertranslation
with the element ``(1, α)``.  We will abuse notation slightly and
refer to the elements ``Λ`` and ``α`` themselves as elements of the
BMS group.

Note that ``𝒮`` is a [normal
subgroup](https://en.wikipedia.org/wiki/Normal_subgroup) of
``\text{BMS}``, but ``ℒ`` is not.  That is, for ``Λ ∈ ℒ`` and ``α ∈
𝒮``, the element ``Λ α Λ⁻¹`` is still an element of ``𝒮``, but ``α Λ
α⁻¹`` is not an element of ``ℒ``.  This fact will be useful later.

The supertranslations form an *abelian* subgroup of ``\text{BMS}``
because

```math
(1, α₂) (1, α₁)
= (1, α₁ + α₂)
= (1, α₂ + α₁)
= (1, α₁) (1, α₂).
```

We have essentially *constructed* the BMS group, but
[Sachs_1962a](@citet) actually *derived* it from the asymptotic metric
conditions, and analyzed the group structure after the fact
[Sachs_1962b](@cite).

[Sachs_1962a](@citet) was the first to describe the BMS group.  (He
referred to it as the "Generalized Bondi-Metzner group", or GBM group;
later authors renamed it the Bondi-Metzner-Sachs group to honor his
contribution.)  His slightly later paper [Sachs_1962b](@cite) was more
specifically about the BMS group itself, and proved some important
properties, including:

  1. The supertranslations form an abelian normal subgroup ``N`` of
     the generalized Bondi-Metzner group; the factor group is
     isomorphic to the orthochronous homogeneous Lorentz group.
  2. The translations form a normal four-dimensional subgroup of the
     proper Bondi Metzner group.
  3. If ``N'`` is a four dimensional normal subgroup of the proper GBM
     group then ``N'`` is contained in the supertranslation group
     ``N``.
  4. The only normal four dimensional subgroup of the GBM group is the
     translation group.

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
