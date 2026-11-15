# Aberration of Gravitational Waves

The approach on this page is *not* the primary method used in this
package; it is retained as an independent second line of reasoning.
The primary method — extracting the ``K`` factor of the Iwasawa
``KAN`` decomposition, which requires no angles, no branches, and no
transcendental functions — is described on [The Lorentz Group](@ref)
page, and is what [`Scri.aberration`](@ref) implements.  The
derivation below works instead with explicit angles on the sphere,
which is more laborious and numerically more delicate, but entirely
independent; it survives as the oracle implementation
(`AberrationOracle.aberration` in `test/aberration.jl`) against which
the primary implementation is validated.

## Setup

Consider two inertial observers, A and B, both of whom pass through
the spacetime origin on their respective worldlines.  B moves at
constant velocity ``v⃗`` relative to A, with speed ``β = |v⃗| < 1``.
Associated to A is a set of mode weights — a spin-weighted function on
the sphere, expressed as a function of A's proper time — describing
either a field that A **emits** (to future null infinity ``ℐ^+``) or a
field that A **receives** (from past null infinity ``ℐ^-``).  The goal
of this section is to find B's description of the same physical field.

We use the sign ``ℐ ∈ \{+1,-1\}`` of [Future and past null
infinity](@ref scri_pm_conventions) to track the two cases:

```math
ℐ = \begin{cases}
  +1 & \text{emitted field ($\mathscr{I}^+$)} \\
  -1 & \text{received field ($\mathscr{I}^-$)}
\end{cases}
```

This single variable turns out to control the sign of every
frame-dependent quantity in the transformation.  In the code it enters
[`Scri.aberration`](@ref) through the boost in its second argument,
`Boost(ℐ * v⃗) * R`, and it is the `emitted` keyword of the test oracle
(`emitted = true` ⟺ ``ℐ = +1``).

## Waves and null wavevectors

A monochromatic plane wave in flat Minkowski spacetime has the form

```math
Φ(x) = A \exp(i\, k_μ x^μ),
```

where the phase ``φ(x) = k_μ x^μ`` is a Lorentz scalar.  Its gradient
``k_μ = ∂_μ φ`` is the **covariant wavevector**.  Inserting this into
the wave equation ``□Φ = 0`` gives

```math
k^μ k_μ = 0,
```

the **null condition**, with metric signature ``{-}{+}{+}{+}``.

### Frame decomposition and the ℐ convention

Given observer A with four-velocity ``t^μ = (1,0,0,0)`` (at rest in
A's frame), we define

```math
ω = -k_μ t^μ = k^t,
```

the frequency as measured by A.  The null condition then forces the
spatial wavevector ``\vec{k} = (k^x, k^y, k^z)`` to satisfy
``|\vec{k}| = ω``.  We write

```math
k^μ = ω(1,\, ℐ\hat{k}),
```

where ``\hat{k}`` is a unit 3-vector and ``ℐ = \pm 1``, as defined
above:

- **``ℐ = +1``** (outgoing, ``ℐ^+``): ``\hat{k}`` points **away from**
  A — it is the propagation direction of the wave.  (Think of
  gravitational waves emitted by A; the chosen portion of the field
  travels outward in direction ``\hat{k}``.)
- **``ℐ = -1``** (incoming, ``ℐ^-``): ``\hat{k}`` points from A
  **toward the source**, and the wave propagates in the opposite
  direction, ``-\hat{k}``.  (Think of a plane wave whose source is far
  away in the direction ``\hat{k}``.)

This is identical to the null section ``σ_ℐ: \hat{k} \mapsto (1,
ℐ\hat{k})`` used on the [BMS group page](@ref bms_group), which for
``ℐ = +1`` is the section introduced on the [null-cones page](@ref
"The null cone and the celestial sphere") in its discussion of the
celestial sphere.  The same rotor language connects
naturally to the [Spacetime Algebra](@extref Quaternionic :doc:`spacetime_algebra`) developed in Quaternionic.jl, where null vectors arise from
combinations of boost and rotation generators in the even subalgebra
of ``\mathrm{Cl}(3,1)``.

## Lorentz transformation of the wavevector

B's four-velocity in A's frame is ``(t_B)^μ = (γ, γ\vec{v})``, where
``γ = 1/\sqrt{1-β^2}``.  B measures the frequency

```math
ω' = -k_μ (t_B)^μ.
```

Expanding with ``k^μ = ω(1, ℐ\hat{k})`` and writing ``\cos Θ = \hat{k}
\cdot \hat{v}`` (where ``\hat{v} = \vec{v}/β`` is the unit boost
direction):

```math
ω' = γω(1 - ℐβ\cos Θ).
```

The angle ``Θ`` is *not* an independent definition; it is fixed by the
inner product of the spatial wavevector with the boost direction,
``\cos Θ = ℐ(k^i \hat{v}_i)/ω``, which comes directly from the
contraction ``k_μ (t_B)^μ``.  When ``ℐ = +1``, this is ``\cos Θ =
\hat{k}\cdot\hat{v}`` with ``\hat{k}`` the propagation direction; when
``ℐ = -1`` it is ``\cos Θ = \hat{k}\cdot\hat{v}`` with ``\hat{k}`` the
direction toward the source.

### The conformal factor

The ratio

```math
κ = \frac{ω}{ω'} = \frac{1}{γ(1 - ℐ\vec{v}\cdot\hat{k})}
```

is exactly the conformal factor introduced on the [null-cones
page](@ref "Lorentz transformations ``ℒ`` and the conformal factor
``κ``"): for ``ℐ = +1`` it reduces to ``κ =
1/[γ(1-\vec{v}\cdot\hat{k})],`` while for ``ℐ = -1`` it becomes ``κ =
1/[γ(1+\vec{v}\cdot\hat{k})]``.

## The aberration formula

In addition to changing the frequency, B's boost changes the
**direction** from which the wave appears to come.  For a boost along
``\hat{z}``, the new direction ``Θ'`` (measured from ``\hat{z}``)
satisfies

```math
\cos Θ' = \frac{\cos Θ - ℐβ}{1 - ℐβ\cos Θ}.
```

The ``ℐ = +1`` case is the standard relativistic aberration formula
[MisnerThorneWheeler_1973](@cite), and the sign of ``ℐ`` extends it to
incoming radiation.

This is derived by applying the Lorentz boost directly to the spatial
wavevector components: with ``k^μ = ω(1, ℐ\sinΘ, 0, ℐ\cos Θ)``
(choosing the boost along ``z`` and the wavevector in the
``x``-``z`` plane, for simplicity), the
boosted spatial ``z``-component is

```math
k'^z = γ(k^z - β k^t) = γω(ℐ\cos Θ - β),
```

which is ``ℐω'\cos Θ'``.  Dividing by ``ω' = γω(1 - ℐβ\cos Θ)`` and
multiplying by ``ℐ`` gives the formula above.

### Half-angle form

The half-angle substitution ``\cos Θ = (1-t^2)/(1+t^2)``, ``t =
\tan(Θ/2)``, allows us to write the aberration formula as

```math
\tan\!\frac{Θ'}{2} = e^{+ℐφ}\,\tan\!\frac{Θ}{2},
\qquad φ = \operatorname{atanh}(β).
```

This compact form, involving the **rapidity** ``φ``, is the one used
in the rotor implementation.

### Geometric sign

The direction of the shift depends on ``ℐ``:

- **``ℐ = +1``** (emitted): ``Θ' > Θ`` for ``β > 0``.  In B's frame
  the wave appears at a *larger* angle from the boost axis than in A's
  frame.  Equivalently, in A's frame the radiation is beamed *toward*
  the boost axis.  This is relativistic beaming familiar from
  astrophysical jets [MisnerThorneWheeler_1973](@cite).

- **``ℐ = -1``** (received): ``Θ' < Θ`` for ``β > 0``.  Moving
  observers see incoming sources shifted *toward* their direction of
  motion — the classical stellar-aberration effect.  This is the
  convention used in Penrose-Rindler Vol. 1 around Eq. (1.3.5)
  [PenroseRindler_1984](@cite), which works on the past celestial
  sphere.

## General boost: the rotor formulation

For a boost in an arbitrary direction ``\vec{v}``, the aberration is a
rotation of ``\hat{k}`` in the plane spanned by ``\hat{k}`` and
``\vec{v}``.  Let ``Θ'`` denote the boosted-frame angle between
``\hat{k}'`` and the boost axis (this is the angle *B* observes), and
let ``Θ`` denote the corresponding rest-frame angle (the angle *A*
uses).  They are related by the half-angle formula

```math
\tan\!\frac{Θ}{2} = e^{-ℐφ}\,\tan\!\frac{Θ'}{2}.
```

The rotation that maps ``\hat{k}'`` to the correct rest-frame
direction is generated by the **aberration rotor**

```math
B' = \exp\!\left(\frac{\hat{k}' \times \vec{v}}{|\hat{k}' \times \vec{v}|}\,\frac{Θ'-Θ}{2}\right).
```

For ``ℐ = +1``, ``Θ' > Θ`` so the exponent is positive; for ``ℐ =
-1``, ``Θ' < Θ`` and the exponent is negative, reversing the rotation.

Crucially, ``B'`` does not merely map the *direction* ``\hat{k}' \to
\hat{k}``:  it also rotates the **tangent frame** at each point on the
sphere, producing the spin-weight phase factor that enters the mode
transformation of spin-weighted functions.  This is why we work with
full rotors rather than unit 3-vectors.

The oracle (`AberrationOracle.aberration` in `test/aberration.jl`)
implements this as the rest-frame rotor ``\operatorname{ab}(R_{\mathrm{pix}},
\vec{v}) = B' R_{\mathrm{pix}}``, where
``R_{\mathrm{pix}}`` is the pixel rotor encoding both the direction
and the tangent frame orientation at a grid point.  The derivation of
``B'`` is given in Appendix C of [Boyle_2015](@cite), Eqs. (C6)–(C8).
The primary implementation, [`Scri.aberration`](@ref), obtains the
same rotor as the ``K`` factor of the ``KAN`` decomposition; the
testitem `"aberration: KAN K factor matches Boyle (2015) oracle"`
verifies that the two agree — as full rotors, so that the
tangent-frame part is pinned down along with the direction — across
random rotors, velocities, and both signs of ``ℐ``.

!!! warning "The oracle's Taylor branch is inaccurate for large β near the poles"

    To remain well-conditioned where ``\hat{k}' × \vec{v} → 0``, this
    implementation switches to a Taylor expansion when ``β\sin Θ' <
    ϵ^{1/3}``.  That expansion is a series in ``β`` alone, so while it
    is excellent for small ``β``, it is *wrong* when ``β`` is large
    and only ``\sin Θ'`` is small — e.g., at ``Θ' = 10^{-7}``, ``β =
    0.99``, it returns a rotor component of ``9.25 × 10^{-10}`` where
    the correct value is ``e^{-ℐφ} Θ'/2 ≈ 3.54 × 10^{-9}``.  The
    ``K``-factor implementation is purely algebraic, has no branches,
    and gets this regime right (verified against `BigFloat`
    evaluations, which keep the oracle in its exact branch).  This
    defect — discovered only when the two independent implementations
    were compared — is one more reason this page's approach was
    demoted to oracle duty.

## Geometric observations and tests

Each of the following observations is a direct consequence of the
formula above.  They are precise enough to serve as tests, and each
maps to a `@testitem` in `test/aberration.jl` or `src/aberration.jl`.

### ``β = 0``: no boost, no aberration

At ``β = 0``, ``φ = 0`` and ``Θ = Θ'`` for both values of ``ℐ``.  The
aberration rotor is the identity, so ``\operatorname{ab}(R_{\mathrm{pix}}, \vec{v}) =
R_{\mathrm{pix}}``
regardless of the boost direction.  This holds for all input rotors
and all choices of ``\hat{v}``.

*Test*: `"aberration: identity at β=0"` (in `src/aberration.jl`; both
signs of ``ℐ``).

### Pole invariance: no tangent rotation along the boost axis

When ``\hat{k}' \parallel \pm\hat{v}``, the cross product ``\hat{k}'
\times \vec{v} = 0`` and the rotation axis in ``B'`` vanishes.
Consequently ``B' = 1`` and the rotor is unchanged.

This covers both the north pole (``\hat{k}' = +\hat{v}``) and the
south pole (``\hat{k}' = -\hat{v}``), and holds for both ``ℐ = +1``
and ``ℐ = -1`` since the cross product is independent of the sign
convention.

*Test*: `"aberration: pole invariance — no tangent rotation along
boost axis"` (in `src/aberration.jl`; both signs of ``ℐ``).

### Equatorial formula: ``\cos Θ = ℐβ``

Take the boost along ``\hat{z}`` and consider an equatorial direction
``\hat{k}'`` (``Θ' = π/2``).  Setting ``\cos Θ' = 0`` in the aberration
formula gives

```math
\cos Θ = ℐβ.
```

For ``ℐ = +1``: the rest-frame direction is in the *northern*
hemisphere, ``\cos Θ = +β > 0``.  For ``ℐ = -1``: the rest-frame
direction is in the *southern* hemisphere, ``\cos Θ = -β < 0``.

*Test*: `"aberration: geometric sign — equatorial pixel: cosΘ = ℐβ"`.

### Azimuthal symmetry

For a boost along ``\hat{z}``, any rotation ``R_z`` about ``\hat{z}``
commutes with the aberration correction:

```math
\operatorname{ab}(R_z R_{\mathrm{pix}},\, \vec{v}) = R_z\, \operatorname{ab}(R_{\mathrm{pix}},\, \vec{v}).
```

This follows because ``R_z`` acts on ``\hat{k}'`` by rotating it about
the boost axis, which leaves the angle ``Θ'`` — and therefore the
magnitude of the aberration — unchanged, while rotating the axis
``\hat{k}' \times \vec{v}`` by the same ``R_z``.

*Test*: `"aberration: azimuthal symmetry — z-rotation commutes with z-boost"`.

### Round-trip: ``ℐ = +1`` and ``ℐ = -1`` are inverses

Applying the ``ℐ = +1`` aberration and then the ``ℐ = -1`` aberration
with the same ``\vec{v}`` recovers the identity:

```math
\operatorname{ab}(\operatorname{ab}(R, \vec{v};\, ℐ=+1),\, \vec{v};\, ℐ=-1) = R.
```

This is immediate from the half-angle formula: ``e^{-ℐφ}`` with ``ℐ =
+1`` followed by ``e^{-ℐφ}`` with ``ℐ = -1`` gives ``e^{-φ} e^{+φ} =
1``.  (For the ``K``-factor implementation the same identity follows
from uniqueness of the Iwasawa decomposition: the second boost cancels
the first, and the leftover ``AN`` factor is absorbed.)

*Test*: `"aberration: round-trip with inverse gives identity"` (both
orders of the signs, plus a cross-check against the independent
`boosted_rotor` helper).

### Wigner rotation: composition of non-collinear boosts

Two successive boosts in different directions do not compose to a pure
boost; their composition in the Spin group combines a boost with a
Wigner rotation.  Specifically, change frame first by ``\vec{v}_2`` and
then by ``\vec{v}_1``, read in the frame the first change produces.
The frame rotor of the two changes together is
``B(\vec{v}_2)\,B(\vec{v}_1)``, and if it is written as
``B(v⃗_{\mathrm{eff}})\,R_W``, as Quaternionic.jl's `vR` decomposition
does, then

```math
\operatorname{ab}\!\left(\operatorname{ab}(R, \vec{v}_1),\, \vec{v}_2\right)
= \operatorname{ab}\!\left(R_W R,\, \vec{v}_{\mathrm{eff}}\right).
```

The Wigner rotation ``R_W`` is a purely spatial rotation that has no
classical analogue; it appears because the Lorentz group is not simply
the direct product of rotations and boosts.

*Test*: `"aberration: two-boost composition matches vR decomposition (Wigner
rotation)"`.

## Working in the unprimed frame

At this point, we are in a good position to understand a subtlety in
the implementation of the transformation.  The field is expressed in
either frame on a series of time slices, and on each time slice the
field components are decomposed in spin-weighted spherical harmonics.
During the transformation, these have to be evaluated at a series of
"pixels" — rotors from ``\mathrm{Spin}(3)`` [Boyle_2016](@cite), which
correspond to the flagpoles and flagplanes of the [null-cones
page](@ref "Flagpoles and flagplanes").  Specifically, there is a
natural choice of pixels in the *transformed* frame which will allow
us to compute the mode weights.  We will denote those pixels as
``𝐑'_p``.  But the field at each such pixel has to be computed from
the original data.  So we need to evaluate it at some corresponding
``𝐑_p``, and the flagpole and flagplane should be the same for both:

```math
\begin{aligned}
𝐑'ₚ \boldsymbol{ℓ}' 𝐑̃'ₚ
&∝ 𝐑ₚ \boldsymbol{ℓ} 𝐑̃ₚ,
\\
𝐑'ₚ \left(𝐱' ∧ \boldsymbol{ℓ}'\right) 𝐑̃'ₚ
&∝ 𝐑ₚ \left(𝐱 ∧ \boldsymbol{ℓ}\right) 𝐑̃ₚ.
\end{aligned}
```

The primed reference vectors are the new legs, which the frame rotor
``Λ = B(v⃗)\,R`` builds from the old ones (see [Frames and
Components](@ref frames_and_components)):

```math
\boldsymbol{ℓ}' = Λ \boldsymbol{ℓ} Λ̃,
\qquad
𝐱' = Λ 𝐱 Λ̃.
```

The pixel rotor ``𝐑'_p`` is built on the primed basis in the same
way.  But the `R′ₚ` actually generated by the code is just an array of
four numbers, and throughout the code the frame of those four numbers
is never explicitly specified.  Nothing in the code that generates the
array refers to either basis, but the pixels are chosen in the new
frame, so the numbers are the components of ``𝐑'_p`` in the primed
basis,

```math
𝐑'ₚ = R'ₚʷ + R'ₚˣ 𝐢' + R'ₚʸ 𝐣' + R'ₚᶻ 𝐤',
```

while the rest of the calculation treats every array as components in
the *unprimed* basis.  The rotor that the array represents there is

```math
r = R'ₚʷ + R'ₚˣ 𝐢 + R'ₚʸ 𝐣 + R'ₚᶻ 𝐤 = Λ̃\, 𝐑'ₚ\, Λ,
```

because the basis bivectors transform like the legs, ``𝐢' =
Λ\,𝐢\,Λ̃``, and so on.  So ``𝐑'_p\,Λ = Λ\,r``, and the flagpole
condition becomes

```math
𝐑'ₚ \boldsymbol{ℓ}' 𝐑̃'ₚ
= \left(Λ\, r\right) \boldsymbol{ℓ} \left(Λ\, r\right)^{\sim}
∝ 𝐑ₚ \boldsymbol{ℓ} 𝐑̃ₚ,
```

and similarly for the flagplane.  The [Lorentz page](@ref "Reversing
the order") argued that all we actually need is the ``K`` factor of
the rotor ``Λ\,r``, so

```math
𝐑ₚ = K\left(Λ\, r\right).
```

The rotor ``r`` is what the Lorentz page writes as ``𝐑'_p``, the
pixel's point frame read in the old components, so ``Λ\,r`` is that
page's point frame ``𝐐``.  Here ``K`` is the function that extracts
the ``K`` factor from a ``KAN`` decomposition.  (Note that the
equation above is in "math" face, while the one below is in "code"
face, to distinguish what happens analytically from what happens
practically.)  In code, ``r`` is just the array `R′ₚ`, so this is

```julia
Rₚ = K(Λℐ * R′ₚ)
```

which is exactly what [`aberration`](@ref Scri.aberration)`(R′ₚ, Λℐ)`
computes, where `Λℐ` is ``Λ`` at ``ℐ⁺``.  At ``ℐ⁻`` the code passes
the parity conjugate ``Λ_ℐ`` = `ℂconj(Λ) = Boost(-v⃗) * R` instead of
``Λ``, as the docstring of `aberration` explains.
Again, the ``A`` and ``N`` factors that the extraction discards are
exactly the ones accounted for elsewhere: the conformal factor ``κ``,
and the component mixing, which re-enters as a null rotation about the
generator ``𝐧``, as discussed on [the tetrad page](@ref "BMS Action on
the Tetrad").
