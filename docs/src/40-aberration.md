# Aberration of gravitational waves

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
ε = \begin{cases}+1 & \text{emitted field (ℐ^+)} \\ -1 & \text{received field (ℐ^-)}\end{cases}
```

This single variable turns out to control the sign of every
frame-dependent quantity in the transformation.  In the code, it is
exposed as the `emitted` keyword argument of [`Scri.aberration`](@ref).

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

where ``\hat{n}`` is a unit 3-vector and ``ε = \pm 1``:

- **``ε = +1``** (outgoing, ``ℐ^+``): ``\hat{n}`` points **away** from A — it is the
  propagation direction of the wave.  (Think of gravitational waves emitted by A; each
  frequency component travels outward in direction ``\hat{n}``.)
- **``ε = -1``** (incoming, ``ℐ^-``): ``\hat{n}`` points **toward** A — it is the
  direction from which the wave arrives.  (Think of a plane wave whose source is far
  away in the direction ``\hat{n}``.)

This is identical to the null section ``σ_ε: \hat{n} \mapsto (1,
ε\hat{n})`` introduced in the [BMS group page](@ref "Lorentz
transformations ``ℒ`` and the conformal factor ``K``") in its
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

Expanding with ``k^μ = ω(1, ε\hat{n})`` and writing ``\cosΘ = \hat{n}
\cdot \hat{v}`` (where ``\hat{v} = \vec{v}/β`` is the unit boost
direction):

```math
ω' = γω(1 - εβ\cosΘ).
```

The angle ``Θ`` is *not* an independent definition: it is the inner
product of the spatial wavevector with the boost direction, ``\cosΘ =
(k^i \hat{v}_i)/ω``, read off directly from the contraction ``k_μ
(t_B)^μ``.  When ``ε = +1``, this is ``\cosΘ = \hat{n}\cdot\hat{v}``
with ``\hat{n}`` the propagation direction; when ``ε = -1`` it is
``\cosΘ = \hat{n}\cdot\hat{v}`` with ``\hat{n}`` the direction toward
the source.

### The conformal factor

The ratio

```math
K = \frac{ω}{ω'} = \frac{1}{γ(1 - ε\vec{v}\cdot\hat{n})}
```

is exactly the conformal factor introduced in the BMS page: for ``ε =
+1`` it reduces to ``K = 1/[γ(1-\vec{v}\cdot\hat{n})]`` while for ``ε
= -1`` it becomes ``K = 1/[γ(1+\vec{v}\cdot\hat{n})]``.

## The aberration formula

In addition to changing the frequency, B's boost changes the
**direction** from which the wave appears to come.  For a boost along
``\hat{z}``, the new direction ``Θ'`` (measured from ``\hat{z}``)
satisfies

```math
\cosΘ' = \frac{\cosΘ - εβ}{1 - εβ\cosΘ}.
```

[⚠️ *Verify*: MTW §22.5 (or a nearby exercise) for the standard form
of this formula; the ``ε = +1`` case is the standard result.]  [⚠️
*Verify*: Schutz §2.7–2.8 derives ``ω' = γω(1-β\cosΘ)`` for ``ε =
+1``; check whether the aberration formula itself also appears there.]

This is derived by applying the Lorentz boost directly to the spatial
wavevector components: with ``k^μ = ω(1, ε\sinΘ, 0, ε\cosΘ)``
(choosing the boost in the ``x``-``z`` plane for simplicity), the
boosted spatial ``z``-component is

```math
k'^z = γ(k^z - β k^t) = γω(ε\cosΘ - β),
```

and dividing by ``ω' = γω(1 - εβ\cosΘ)`` gives the formula above.

### Half-angle form

The half-angle substitution ``\cosΘ = (1-t^2)/(1+t^2)``, ``t =
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
``\vec{v}``.  Let ``Θ̑`` denote the boosted-frame angle between
``\hat{n}'`` and the boost axis (this is the angle *B* observes), and
let ``Θ`` denote the corresponding rest-frame angle (the angle *A*
uses).  They are related by the half-angle formula

```math
\tan\!\frac{Θ}{2} = e^{-εφ}\,\tan\!\frac{Θ̑}{2}.
```

The rotation that maps ``\hat{n}'`` to the correct rest-frame
direction is generated by the **aberration rotor**

```math
B' = \exp\!\left(\frac{\hat{n}' \times \vec{v}}{|\hat{n}' \times \vec{v}|}\,\frac{Θ̑-Θ}{2}\right).
```

For ``ε = +1``, ``Θ̑ > Θ`` so the exponent is positive; for ``ε =
-1``, ``Θ̑ < Θ`` and the exponent is negative, reversing the rotation.

Crucially, ``B'`` does not merely map the *direction* ``\hat{n}' \to
\hat{n}``:  it also rotates the **tangent frame** at each point on the
sphere, producing the spin-weight phase factor that enters the mode
transformation of spin-weighted functions.  This is why we work with
full rotors rather than unit 3-vectors.

The function [`Scri.aberration`](@ref) implements this for the product ``R' = B'
R_{\mathrm{pix}}``, where ``R_{\mathrm{pix}}`` is the pixel rotor
encoding both the direction and the tangent frame orientation at a
grid point.  The derivation of ``B'`` is given in Appendix C of
[Boyle_2015](@cite), Eqs. (C6)–(C8).

## Geometric observations and tests

Each of the following observations is a direct consequence of the
formula above.  They are precise enough to serve as tests, and each
maps to a `@testitem` in `test/test_aberration.jl`.

### ``β = 0``: no boost, no aberration

At ``β = 0``, ``φ = 0`` and ``Θ = Θ̑`` for both values of ``ε``.  The
aberration rotor is the identity, so ``R' = R_{\mathrm{pix}}``
regardless of the boost direction.  This holds for all input rotors
and all choices of ``\hat{v}``.

*Tests*: `"aberration: identity at β=0"` and `"aberration: emitted=false — identity at
β=0"`.

### Pole invariance: no tangent rotation along the boost axis

When ``\hat{n}' \parallel \pm\hat{v}``, the cross product ``\hat{n}'
\times \vec{v} = 0`` and the rotation axis in ``B'`` vanishes.
Consequently ``B' = 1`` and the rotor is unchanged.

This covers both the north pole (``\hat{n}' = +\hat{v}``) and the
south pole (``\hat{n}' = -\hat{v}``), and holds for both ``ε = +1``
and ``ε = -1`` since the cross product is independent of the sign
convention.

*Tests*: `"aberration: pole invariance"` and `"aberration: emitted=false — pole
invariance"`.

### Equatorial formula: ``\cosΘ = εβ``

Take the boost along ``\hat{z}`` and consider an equatorial direction
``\hat{n}'`` (``Θ̑ = π/2``).  Setting ``\cosΘ' = 0`` in the aberration
formula gives

```math
\cosΘ = εβ.
```

For ``ε = +1``: the rest-frame direction is in the *northern*
hemisphere, ``\cosΘ = +β > 0``.  For ``ε = -1``: the rest-frame
direction is in the *southern* hemisphere, ``\cosΘ = -β < 0``.

*Tests*: `"aberration: geometric sign — equatorial pixel maps closer to boost
axis (future ℐ⁺)"` and `"aberration: emitted=false — equatorial pixel maps
farther from boost axis (past ℐ⁻)"`.

### Azimuthal symmetry

For a boost along ``\hat{z}``, any rotation ``R_z`` about ``\hat{z}``
commutes with the aberration correction:

```math
R'(R_z R_{\mathrm{pix}},\, \vec{v}) = R_z\, R'(R_{\mathrm{pix}},\, \vec{v}).
```

This follows because ``R_z`` acts on ``\hat{n}'`` by rotating it about
the boost axis, which leaves the angle ``Θ̑`` — and therefore the
magnitude of the aberration — unchanged, while rotating the axis
``\hat{n}' \times \vec{v}`` by the same ``R_z``.

*Test*: `"aberration: azimuthal symmetry — z-rotation commutes with z-boost"`.

### Round-trip: ``ε = +1`` and ``ε = -1`` are inverses

Applying the ``ε = +1`` aberration and then the ``ε = -1`` aberration
with the same ``\vec{v}`` recovers the identity:

```math
R'(R'(R, \vec{v};\, \text{emitted}=\mathtt{true}),\, \vec{v};\, \text{emitted}=\mathtt{false}) = R.
```

This is immediate from the half-angle formula: ``e^{-εφ}`` with ``ε =
+1`` followed by ``e^{-εφ}`` with ``ε = -1`` gives ``e^{-φ} e^{+φ} =
1``.

*Tests*: `"aberration: round-trip with inverse gives identity"` (using the old
`boosted_rotor` helper) and `"aberration: emitted=false round-trip gives
identity"` (using both keyword values).

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
