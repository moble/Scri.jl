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

We introduce a sign variable ``ε ∈ \{+1,-1\}`` to track the two cases:

```math
ε = \begin{cases}
  +1 & \text{emitted field ($\mathscr{I}^+$)} \\
  -1 & \text{received field ($\mathscr{I}^-$)}
\end{cases}
```

This single variable turns out to control the sign of every
frame-dependent quantity in the transformation.  It is the same sign
as the ``ℐ`` of [Future and past null infinity](@ref
scri_pm_conventions); in the code it is the `ℐ` argument of
[`Scri.aberration`](@ref), and the `emitted` keyword of the test
oracle (`emitted = true` ⟺ ``ε = +1``).

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

### Frame decomposition and the ε convention

Given observer A with four-velocity ``t^μ = (1,0,0,0)`` (at rest in
A's frame), we define

```math
ω = -k_μ t^μ = k^t,
```

the frequency as measured by A.  The null condition then forces the
spatial wavevector ``\vec{k} = (k^x, k^y, k^z)`` to satisfy
``|\vec{k}| = ω``.  We write

```math
k^μ = ω(1,\, ε\hat{n}),
```

where ``\hat{n}`` is a unit 3-vector and ``ε = \pm 1``, as defined
above:

- **``ε = +1``** (outgoing, ``ℐ^+``): ``\hat{n}`` points **away from**
  A — it is the propagation direction of the wave.  (Think of
  gravitational waves emitted by A; the chosen portion of the field
  travels outward in direction ``\hat{n}``.)
- **``ε = -1``** (incoming, ``ℐ^-``): ``\hat{n}`` points **toward** A
  — still the propagation direction of the wave.  (Think of a plane
  wave whose source is far away in the direction ``-\hat{n}``.)

This is identical to the null section ``σ_ε: \hat{n} \mapsto (1,
ε\hat{n})`` introduced in the [BMS group page](@ref "Lorentz
transformations ``ℒ`` and the conformal factor ``κ``") in its
discussion of the celestial sphere.  The same rotor language connects
naturally to the [Spacetime Algebra](@extref Quaternionic :doc:`spacetime_algebra`) developed in Quaternionic.jl, where null vectors arise from
combinations of boost and rotation generators in the even subalgebra
of ``\mathrm{Cl}(3,1)``.

## Lorentz transformation of the wavevector

B's four-velocity in A's frame is ``(t_B)^μ = (γ, γ\vec{v})``, where
``γ = 1/\sqrt{1-β^2}``.  B measures the frequency

```math
ω' = -k_μ (t_B)^μ.
```

Expanding with ``k^μ = ω(1, ε\hat{n})`` and writing ``\cos Θ = \hat{n}
\cdot \hat{v}`` (where ``\hat{v} = \vec{v}/β`` is the unit boost
direction):

```math
ω' = γω(1 - εβ\cos Θ).
```

The angle ``Θ`` is *not* an independent definition: it is the inner
product of the spatial wavevector with the boost direction, ``\cos Θ =
(k^i \hat{v}_i)/ω``, read off directly from the contraction ``k_μ
(t_B)^μ``.  When ``ε = +1``, this is ``\cos Θ = \hat{n}\cdot\hat{v}``
with ``\hat{n}`` the propagation direction; when ``ε = -1`` it is
``\cos Θ = \hat{n}\cdot\hat{v}`` with ``\hat{n}`` the direction toward
the source.

### The conformal factor

The ratio

```math
κ = \frac{ω}{ω'} = \frac{1}{γ(1 - ε\vec{v}\cdot\hat{n})}
```

is exactly the conformal factor introduced in the BMS page: for ``ε =
+1`` it reduces to ``κ = 1/[γ(1-\vec{v}\cdot\hat{n})]`` while for ``ε
= -1`` it becomes ``κ = 1/[γ(1+\vec{v}\cdot\hat{n})]``.

## The aberration formula

In addition to changing the frequency, B's boost changes the
**direction** from which the wave appears to come.  For a boost along
``\hat{z}``, the new direction ``Θ'`` (measured from ``\hat{z}``)
satisfies

```math
\cos Θ' = \frac{\cos Θ - εβ}{1 - εβ\cos Θ}.
```

[⚠️ *Verify*: MTW §22.5 (or a nearby exercise) for the standard form
of this formula; the ``ε = +1`` case is the standard result.]  [⚠️
*Verify*: Schutz §2.7–2.8 derives ``ω' = γω(1-β\cos Θ)`` for ``ε =
+1``; check whether the aberration formula itself also appears there.]

This is derived by applying the Lorentz boost directly to the spatial
wavevector components: with ``k^μ = ω(1, ε\sinΘ, 0, ε\cos Θ)``
(choosing the boost in the ``x``-``z`` plane for simplicity), the
boosted spatial ``z``-component is

```math
k'^z = γ(k^z - β k^t) = γω(ε\cos Θ - β),
```

and dividing by ``ω' = γω(1 - εβ\cos Θ)`` gives the formula above.

### Half-angle form

The half-angle substitution ``\cos Θ = (1-t^2)/(1+t^2)``, ``t =
\tan(Θ/2)``, turns the aberration formula into

```math
\tan\!\frac{Θ'}{2} = e^{-εφ}\,\tan\!\frac{Θ}{2},
\qquad φ = \operatorname{atanh}(β).
```

This compact form, involving the **rapidity** ``φ``, is the one used
in the rotor implementation.

### Geometric sign

The direction of the shift depends on ``ε``:

- **``ε = +1``** (emitted): ``Θ' > Θ`` for ``β > 0``.  In B's frame
  the wave appears at a *larger* angle from the boost axis than in A's
  frame.  Equivalently, in A's frame the radiation is beamed *toward*
  the boost axis.  This is relativistic beaming familiar from
  astrophysical jets [MisnerThorneWheeler_1973](@cite).

- **``ε = -1``** (received): ``Θ' < Θ`` for ``β > 0``.  Moving
  observers see incoming sources shifted *toward* their direction of
  motion — the classical stellar-aberration effect.  This is the
  convention used in Penrose-Rindler Vol. 1 around Eq. (1.3.5)
  [PenroseRindler_1984](@cite), which works on the past celestial
  sphere.  [⚠️ *Verify*: confirm Eq. (1.3.5) is for incoming null
  vectors, giving opposite sign to the ``ε = +1`` case.]

## General boost: the rotor formulation

For a boost in an arbitrary direction ``\vec{v}``, the aberration is a
rotation of ``\hat{n}`` in the plane spanned by ``\hat{n}`` and
``\vec{v}``.  Let ``Θ'`` denote the boosted-frame angle between
``\hat{n}'`` and the boost axis (this is the angle *B* observes), and
let ``Θ`` denote the corresponding rest-frame angle (the angle *A*
uses).  They are related by the half-angle formula

```math
\tan\!\frac{Θ}{2} = e^{-εφ}\,\tan\!\frac{Θ'}{2}.
```

The rotation that maps ``\hat{n}'`` to the correct rest-frame
direction is generated by the **aberration rotor**

```math
B' = \exp\!\left(\frac{\hat{n}' \times \vec{v}}{|\hat{n}' \times \vec{v}|}\,\frac{Θ'-Θ}{2}\right).
```

For ``ε = +1``, ``Θ' > Θ`` so the exponent is positive; for ``ε =
-1``, ``Θ' < Θ`` and the exponent is negative, reversing the rotation.

Crucially, ``B'`` does not merely map the *direction* ``\hat{n}' \to
\hat{n}``:  it also rotates the **tangent frame** at each point on the
sphere, producing the spin-weight phase factor that enters the mode
transformation of spin-weighted functions.  This is why we work with
full rotors rather than unit 3-vectors.

The oracle (`AberrationOracle.aberration` in `test/aberration.jl`)
implements this for the product ``R' = B' R_{\mathrm{pix}}``, where
``R_{\mathrm{pix}}`` is the pixel rotor encoding both the direction
and the tangent frame orientation at a grid point.  The derivation of
``B'`` is given in Appendix C of [Boyle_2015](@cite), Eqs. (C6)–(C8).
The primary implementation, [`Scri.aberration`](@ref), obtains the
same rotor as the ``K`` factor of the ``KAN`` decomposition; the
testitem `"aberration: KAN K factor matches Boyle (2015) oracle"`
verifies that the two agree — as full rotors, so that the
tangent-frame part is pinned down along with the direction — across
random rotors, velocities, and both signs of ``ε``.

!!! warning "The oracle's Taylor branch is inaccurate for large β near the poles"

    To remain well-conditioned where ``\hat{n}' × \vec{v} → 0``, this
    implementation switches to a Taylor expansion when ``β\sin Θ' <
    ϵ^{1/3}``.  That expansion is a series in ``β`` alone, so while it
    is excellent for small ``β``, it is *wrong* when ``β`` is large
    and only ``\sin Θ'`` is small — e.g., at ``Θ' = 10^{-7}``, ``β =
    0.99``, it returns a rotor component of ``9.25 × 10^{-10}`` where
    the correct value is ``e^{-εφ} Θ'/2 ≈ 3.54 × 10^{-9}``.  The
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

At ``β = 0``, ``φ = 0`` and ``Θ = Θ'`` for both values of ``ε``.  The
aberration rotor is the identity, so ``R' = R_{\mathrm{pix}}``
regardless of the boost direction.  This holds for all input rotors
and all choices of ``\hat{v}``.

*Test*: `"aberration: identity at β=0"` (in `src/aberration.jl`; both
signs of ``ε``).

### Pole invariance: no tangent rotation along the boost axis

When ``\hat{n}' \parallel \pm\hat{v}``, the cross product ``\hat{n}'
\times \vec{v} = 0`` and the rotation axis in ``B'`` vanishes.
Consequently ``B' = 1`` and the rotor is unchanged.

This covers both the north pole (``\hat{n}' = +\hat{v}``) and the
south pole (``\hat{n}' = -\hat{v}``), and holds for both ``ε = +1``
and ``ε = -1`` since the cross product is independent of the sign
convention.

*Test*: `"aberration: pole invariance — no tangent rotation along
boost axis"` (in `src/aberration.jl`; both signs of ``ε``).

### Equatorial formula: ``\cos Θ = εβ``

Take the boost along ``\hat{z}`` and consider an equatorial direction
``\hat{n}'`` (``Θ' = π/2``).  Setting ``\cos Θ' = 0`` in the aberration
formula gives

```math
\cos Θ = εβ.
```

For ``ε = +1``: the rest-frame direction is in the *northern*
hemisphere, ``\cos Θ = +β > 0``.  For ``ε = -1``: the rest-frame
direction is in the *southern* hemisphere, ``\cos Θ = -β < 0``.

*Test*: `"aberration: geometric sign — equatorial pixel: cosΘ = ℐβ"`.

### Azimuthal symmetry

For a boost along ``\hat{z}``, any rotation ``R_z`` about ``\hat{z}``
commutes with the aberration correction:

```math
R'(R_z R_{\mathrm{pix}},\, \vec{v}) = R_z\, R'(R_{\mathrm{pix}},\, \vec{v}).
```

This follows because ``R_z`` acts on ``\hat{n}'`` by rotating it about
the boost axis, which leaves the angle ``Θ'`` — and therefore the
magnitude of the aberration — unchanged, while rotating the axis
``\hat{n}' \times \vec{v}`` by the same ``R_z``.

*Test*: `"aberration: azimuthal symmetry — z-rotation commutes with z-boost"`.

### Round-trip: ``ε = +1`` and ``ε = -1`` are inverses

Applying the ``ε = +1`` aberration and then the ``ε = -1`` aberration
with the same ``\vec{v}`` recovers the identity:

```math
R'(R'(R, \vec{v};\, ε=+1),\, \vec{v};\, ε=-1) = R.
```

This is immediate from the half-angle formula: ``e^{-εφ}`` with ``ε =
+1`` followed by ``e^{-εφ}`` with ``ε = -1`` gives ``e^{-φ} e^{+φ} =
1``.  (For the ``K``-factor implementation the same identity follows
from uniqueness of the Iwasawa decomposition: the second boost cancels
the first, and the leftover ``AN`` factor is absorbed.)

*Test*: `"aberration: round-trip with inverse gives identity"` (both
orders of the signs, plus a cross-check against the independent
`boosted_rotor` helper).

### Wigner rotation: composition of non-collinear boosts

Two successive boosts in different directions do not compose to a pure
boost; their composition in the Spin group is a boost followed by a
Wigner rotation.  Specifically, if the composite Lorentz
transformation decomposes as ``(v⃗_{\mathrm{eff}}, R_W)`` (in the
sense of Quaternionic.jl's `vR` decomposition), then

```math
R'\!\left(R'(R, \vec{v}_1),\, \vec{v}_2\right)
= R'\!\left(R_W R,\, \vec{v}_{\mathrm{eff}}\right).
```

The Wigner rotation ``R_W`` is a purely spatial rotation that has no
classical analogue; it appears because the Lorentz group is not simply
the direct product of rotations and boosts.

*Test*: `"aberration: two-boost composition matches vR decomposition (Wigner
rotation)"`.
