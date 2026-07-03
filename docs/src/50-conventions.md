# Conventions

As with all of GR, the literature on spacetime asymptotics is,
unfortunately, full of subtly varying conventions for objects that are
fundamentally the same.  The metric signature, the tetrad
normalization, curvature quantities, and even the ``ð`` operator all
differ between references.  These choices are not just internal; they
actually affect the transformation laws of the fields.  For example,
the Newman–Penrose Weyl components ``Ψ_n`` are defined as contractions
of the Weyl tensor with various elements of the tetrad, ``(l, m,
\bar{m}, n)``.  But these are null vectors, so the normalizations are
not fixed and naturally (even accounting for different letters used to
represent the elements) different authors choose different
normalizations for the tetrad elements.  For example, two reasonable
choices found in the literature require

```math
\begin{aligned}
l \qquad &↔ \qquad l\sqrt{2}, \\
\hphantom{\sqrt{2}}n \qquad &↔ \qquad n/\sqrt{2},
\end{aligned}
```

which changes the transformation law for ``ψ₃`` in a nontrivial way:

```math
ψ₃' = \frac{e^{-iλ}}{κ³} \left[ψ₃ + \frac{ðα}{2κ} ψ₄\right]
\qquad ↔ \qquad
ψ₃' = \frac{e^{-iλ}}{κ³} \left[ψ₃ + \frac{ðα}{2κ} \frac{ψ₄}{\sqrt{2}}\right].
```

Assuming wrong conventions about the input data can lead to incorrect
results, so it is important to know which conventions are being used.
This package provides a `Conventions` type to specify the conventions
of the input data, so that the correct transformation laws are
applied.  The overhead of using non-default conventions is very small.

---

An important source of these conventions is [BoyleEtAl_2019](@citet),
which describes conventions used for SXS waveforms in Appendix C.  BMS
conventions are described in [MitmanEtAl_2024](@citet).
[Iozzo_2021](@citet) creates a framework for comparing conventions
across the literature.

The conventions for the geometric algebra are described in the
documentation of `Quaternionic.jl`, [here for the
fundamentals](https://moble.github.io/Quaternionic.jl/stable/geometric_algebra/)
and [here specifically for the spacetime
algebra](https://moble.github.io/Quaternionic.jl/stable/spacetime_algebra/).

!!! warn

    Note LISA Rosetta Stone when it comes out.

Nominal values for solar and planetary quantities are given in
[PršaEtAl_2016](@citet).  Other constants can be found in the
[`lisaconstants` package](https://pypi.org/project/lisaconstants/),
which otherwise integrates with [`astropy`](https://www.astropy.org/).

## Metric, curvature, and perturbations

The signature is ``{-}{+}{+}{+}``, and we use units where ``G=c=1``.

The Christoffel symbols and the Riemann, Ricci, and Weyl tensors
follow the Misner-Thorne-Wheeler conventions
[MisnerThorneWheeler_1973](@cite) — Eqs. (14.36), (8.44), (8.47), and
(13.50) of that reference, respectively:

```math
\begin{gathered}
{\Gamma^a}_{bc} = \frac{1}{2} g^{ad} \bigl( \partial_b g_{cd} + \partial_c g_{bd} - \partial_d g_{bc} \bigr), \\
{R^a}_{bcd} = \partial_c {\Gamma^a}_{bd} - \partial_d {\Gamma^a}_{bc} + {\Gamma^a}_{ce} {\Gamma^e}_{bd} - {\Gamma^a}_{de} {\Gamma^e}_{bc}, \\
R_{ab} = {R^c}_{acb}, \\
C_{abcd} = R_{abcd} - \frac{1}{2} (g_{ac} R_{bd} - g_{ad} R_{bc} + g_{bd} R_{ac} - g_{bc} R_{ad}) + \frac{1}{6} R (g_{ac} g_{bd} - g_{ad} g_{bc}).
\end{gathered}
```

The Newman-Penrose Weyl components are defined as

```math
\begin{aligned}
\Psi_0 &= C_{abcd} \ell^a m^b \ell^c m^d, \\
\Psi_1 &= C_{abcd} \ell^a n^b \ell^c m^d, \\
\Psi_2 &= C_{abcd} \ell^a m^b \bar{m}^c n^d, \\
\Psi_3 &= C_{abcd} \ell^a n^b \bar{m}^c n^d, \\
\Psi_4 &= C_{abcd} n^a \bar{m}^b n^c \bar{m}^d.
\end{aligned}
```

The (Maxwell/Faraday) field-strength tensor ``F_{ab}`` is similarly
decomposed into components as

```math
\begin{aligned}
φ_0 &= F_{ab}\, ℓ^a m^b, \\
φ_1 &= \tfrac{1}{2} F_{ab}\, (ℓ^a n^b + \bar{m}^a m^b), \\
φ_2 &= F_{ab}\, \bar{m}^a n^b.
\end{aligned}
```

These are consistent with the Weyl components above: each ``φ_n``
carries spin weight ``1-n`` (just as ``Ψ_n`` carries ``2-n``), and
stepping down the tower swaps an ``ℓ``-type dyad slot (``o``) for an
``n``-type one (``ι``).  The overall sign is the ``c_φ`` of the
[convention table](@ref "Convention parameters") below, and —
paralleling the Weyl conversion ``∝ (c_l c_m)^{2-n}`` — the
inter-convention factor is

```math
φ_n^{[X]} = c_φ\, (c_l c_m)^{1-n}\, φ_n^{[\mathrm{SpEC}]},
```

with the dyad scaling ``c_l c_m`` raised to the spin weight
``1-n`` and no Riemann-sign factor (the field strength does not see the
curvature convention).

The metric perturbation is defined as

```math
h_{ab} = g_{ab} - \eta_{ab},
```

where ``\eta_{ab}`` is the Minkowski metric.  The strain components
are defined as

```math
\begin{aligned}
h_+ &= \frac{1}{2} (h_{\hat{\theta}\hat{\theta}} - h_{\hat{\phi}\hat{\phi}}), \\
h_\times &= h_{\hat{\theta}\hat{\phi}}, \\
h &= h_+ - i h_\times.
\end{aligned}
```

where the hats indicate orthonormal components in the spherical basis.
Near *future* null infinity, we have the asymptotic relation

```math
\Psi_4 \sim -\ddot{h},
```

where the dots indicate time derivatives.

## The eth operator

Throughout this package, ``ð`` is the **Newman–Penrose** eth (the
spin-raising operator), *not* the Geroch–Held–Penrose (GHP) one.
Restricted to the unit round sphere (with the boost weight dropping
out), the two differ by a factor of ``\sqrt{2}``:

```math
ð_{\mathrm{NP}} = \sqrt{2}\, ð_{\mathrm{GHP}}.
```

While spin-weighted spherical functions [*cannot actually be
defined*](@cite Boyle_2016) on the sphere ``𝕊²`` itself, we can often
just about get away with writing them as functions on *coordinates
over the sphere*.  This is the standard approach in the literature,
and as such the Newman–Penrose eth is defined as acting on a quantity
``{}_s f`` of spin weight ``s`` via

```math
ð\, {}_s f = -(\sin θ)^{s}\left(∂_θ + \frac{i}{\sin θ}∂_ϕ\right)
\left[(\sin θ)^{-s}\, {}_s f\right],
```

raising the spin weight by one; equivalently, on the spin-weighted
spherical harmonics,

```math
ð\, {}_s Y_{ℓ,m} = \sqrt{(ℓ-s)(ℓ+s+1)}\; {}_{s+1} Y_{ℓ,m}.
```

For a spin-0 function this is just ``ð f = -\left(∂_θ + \frac{i}{\sin
θ}∂_ϕ\right) f``, which ties ``ð`` directly to the angular dyad ``m``
of [the standard tetrad](@ref "The standard tetrad").  Since ``m̃ =
\frac{1}{\sqrt{2}}\left(∂_θ + \frac{i}{\sin θ}∂_ϕ\right)``,

```math
ð f = -\sqrt{2}\, m̃(f),
\qquad
ð̄ f = -\sqrt{2}\, m̄̃(f)
\qquad (s = 0).
```

Equivalently the GHP eth is simply ``ð_{\mathrm{GHP}} f = -m̃(f)``.
This relation and the tetrad normalizations are what fix the factors
of ``\sqrt{2}`` — and ultimately the ``1/2`` in the [Weyl mixing
parameter](@ref "BMS action on fields") — whenever ``ð`` of a
coordinate is re-expressed through the tetrad.

## Convention parameters

Every quantity above is defined in one specific convention — the
**SpEC** convention together with the **Newman–Penrose** ``ð`` — which
is what `Scri.jl` computes with natively.  Other codes and papers use
different signs and scales; the differences are captured by a small set
of ``c``-parameters (the `Scri.Conventions` struct), following the
conventions appendix and [Iozzo_2021](@citet).  The defaults reproduce
the package's native convention.

| parameter | meaning               | default |
|:----------|:----------------------|:-------:|
| ``c_s``   | metric signature      | ``1``   |
| ``c_l``   | ``l``-leg scale       | ``1``   |
| ``c_m``   | ``m``-leg scale       | ``1``   |
| ``c_R``   | Riemann sign          | ``1``   |
| ``c_Ψ``   | Weyl sign             | ``1``   |
| ``c_σ``   | shear sign            | ``1``   |
| ``c_h``   | strain scale          | ``1``   |
| ``c_φ``   | Faraday sign          | ``1``   |
| ``c_ð``   | eth coefficient       | ``1``   |
| ``c_α``   | time-law sign         | ``1``   |

Note that ``c_m`` is the spin-phase *factor* ``e^{iΘ}`` of the ``m``
leg, not the angle ``Θ``: a phase of ``Θ = π`` is set with ``c_m =
-1``.  This keeps the common conventions at exactly ``±1``, which
`Scri.jl` stores as the type-level singletons `One`/`MinusOne` so that
a convention factor of ``±1`` compiles away (to at most a negation)
wherever a convention is used in a formula; ``c_m`` may still be any
real or complex number for an exotic convention.  ``c_α`` is the sign
in the time-transformation law ``t′ = κ(t − c_α α)`` (see [The BMS
group](@ref bms_group)); it is a convention like the others.

Two defaults deserve emphasis.  ``c_s=+1`` is the ``−+++`` signature
fixed [above](@ref "Metric, curvature, and perturbations"); note that
the SpEC convention's ``c_s=+1`` already *implies* ``−+++``.  And
``c_ð=1`` is the Newman–Penrose ``ð`` of [The eth operator](@ref), in
which ``ð`` carries **no** ``\sqrt2`` — ``ð f = -(∂_θ + i\cscθ\,∂_φ)f``
on a spin-0 ``f`` (a Geroch–Held–Penrose ``ð`` would be ``c_ð =
1/\sqrt2``).

Inter-convention conversion of the Weyl components and the strain (from
the appendix) is

```math
ψ_n^{[X]} = c_s c_Ψ c_R\,(c_l c_m)^{2-n}\, ψ_n^{[\mathrm{SpEC}]},
\qquad
h^{[X]} = c_s c_h^{-1} c_m^{-2}\, h^{[\mathrm{SpEC}]},
```

A *transformation* that stays within one convention is more economical
than these inter-convention factors suggest: as worked out in ["BMS
action on fields"](@ref convention_dependence_tetrad_future), the
overall signs ``c_s, c_Ψ, c_R, c_φ`` cancel out of the (homogeneous)
peeling towers, leaving only the **dyad scaling ``c_l c_m``**,
which rescales the mixing parameter ``b``; and the eth coefficient
``c_ð`` does not enter a transform at all — it changes only how ``ð``
and the (geometric) ``b`` and shear shift are *written*, never their
values, which the code computes with the Newman–Penrose ``ð``
natively.  Only the *inhomogeneous* shear/strain shifts carry a full
conversion factor (``F_σ``, ``F_h``).

## The Null Cone and Transformations

A general BMS element is the composition of three simpler
transformations, which we choose to apply in this order:
supertranslation ``α`` first, then spatial rotation ``R``, then boost
``b``.  Each of these three parts is specified in the coordinate
system of the (inertial) initial observer ``A``.  That is, we will
describe these transformations in turn as if they were the only one
being applied.

The supertranslation ``α`` shifts the (retarded or advanced) time
coordinate by a real, angle-dependent function of the direction in
``A``'s frame.  We parameterize whichever part of ``ℐ`` we are
interested in with the original time coordinate ``t`` and the null ray
``𝐧`` from the observer to that point.  The effect of the
supertranslation is

```math
t'(t, 𝐧) = t - εᵅ α(𝐧),
```

where ``εᵅ = +1`` for a supertranslation of future null infinity and
``εᵅ = -1`` for a supertranslation of past null infinity.  The null
ray is unchanged by the supertranslation, so the angular coordinates
are unchanged.

We represent the spatial rotation ``R`` as a rotor ``R ∈
\mathrm{Spin}(3)`` that rotates the angular coordinate system.  The
spatial basis vectors of observer ``B`` are related to those of ``A``
by ``𝐞' = R 𝐞 R̃`` (where ``𝐞`` stands in for any basis vector), so
a null ray at ``𝐧`` in ``A``'s frame appears at ``R 𝐧 R̃`` in
``B``'s frame.  The retarded time is unchanged.

Finally, the boost is parameterized by a velocity ``v⃗`` expressed in
``A``'s frame.  Observer ``A``'s worldline is parallel to the time
axis ``𝐭``.  Observer ``B``'s worldline is parallel to

```math
𝐭' = γ(𝐭 + v⃗),
```

which lies inside the future light cone and has positive components
along both ``𝐭`` and ``v⃗``.  The boost rescales the retarded time by
the conformal factor ``κ``:

```math
t' = κ(𝐧)\, t,
\qquad
κ(𝐧) = \frac{1}{γ(1 - ε^ℐ\, v⃗ ⋅ 𝐧)},
```

where ``𝐧`` here is the null ray *before* any boost or rotation has
been applied, and ``ε^ℐ = ±1`` selects future or past null infinity (see
[Future and past null infinity](@ref scri_pm_conventions) below).  The full BMS action of
the element ``(Λ = b∘R, α)`` on ``A``'s coordinates ``(t, 𝐧)`` is then

```math
t' = κ(Λ, 𝐧) \bigl[t - εᵅ α(𝐧)\bigr],
\qquad
𝐧' = \frac{Λ𝐧}{(Λ𝐧)^0},
```

where ``κ(Λ, 𝐧) = n⁰/(Λ𝐧)⁰`` is the conformal factor for the
combined Lorentz transformation ``Λ`` evaluated on the null ray, which
is given in $A$'s coordinate frame as ``𝐧 = (1, ε^ℐ n̂)`` up to
normalization.

## [Future and past null infinity](@id scri_pm_conventions)

Every transformation above lives on a chosen piece of ``ℐ``, and most
expressions carry a sign ``ε^ℐ`` that records *which* piece:

```math
ε^ℐ = \begin{cases} +1 & \text{future null infinity } ℐ⁺ \text{ (outgoing radiation)} \\ -1 & \text{past null infinity } ℐ⁻ \text{ (incoming radiation)}. \end{cases}
```

The subtlety is how the celestial sphere ``𝕊²`` is *labeled* at each
end.  At ``ℐ⁺`` we label a point by the direction in which outgoing
radiation **propagates**, so the section is the future-pointing null
ray ``𝐧 = (1, n̂)``.  At ``ℐ⁻`` we follow
[PenroseRindler_1984](@citet), who work on the observer's **past light
cone** — the sphere of directions from which light **arrives** (the
astronomer's sky).  To see a source you look *opposite* to the light's
direction of travel, so the ``ℐ⁻`` section is ``𝐧 = (1, -n̂)``, i.e.
the labeling is **antipodal** to the ``ℐ⁺`` one.  Equivalently, a
single free null geodesic of Minkowski space joins a point ``n̂`` of
``ℐ⁺`` to the antipodal point ``-n̂`` of ``ℐ⁻``; identifying the two
spheres this way is the *antipodal matching* used to relate ``ℐ⁺`` and
``ℐ⁻`` (e.g., in the analysis of gravitational scattering and soft
theorems [Strominger_2014, Strominger_2017](@cite)).

Both choices are captured by writing the section as ``𝐧 = (1, ε^ℐ
n̂)``, which is the single place the convention enters.  Propagating
it through the conformal factor gives

```math
κ(𝐧) = \frac{1}{γ(1 - ε^ℐ\, v⃗ ⋅ n̂)}.
```

The sign flip is exactly the past-cone Doppler factor: an observer
boosted with velocity ``v⃗`` sees radiation arriving from ``n̂``
shifted by ``γ(1 + v⃗ ⋅ n̂)`` (blueshift toward the direction of
motion), which is ``κ^{-1}`` at ``ε^ℐ = -1``.  This is the same sign
difference Penrose and Rindler note at their Eq. (1.3.5), and it is
the convention implemented throughout this package — in
[`aberration`](@ref Scri.aberration) (its `emitted` keyword), in the
conformal factor used by [`transform!`](@ref), in the ``ε^ℐ`` type
parameter of [`DataComponents`](@ref Scri.DataComponents) (data live
on a fixed null infinity), and in the ``ε^ℐ`` argument of the
[`BMS`](@ref) operations.  Because ``κ`` itself is a spin-0
(orientation-blind) quantity, the choice also governs the *handedness*
of the spin-weighted ``ð`` operator that appears in the
component-mixing law, which is why the Weyl/Faraday peeling tower runs
in the opposite direction at ``ℐ⁻`` (see [Computing ``ðt'/κ``](@ref
computing_eth_tprime_over_kappa) and [`DataComponents`](@ref
Scri.DataComponents)).
