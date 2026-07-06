# Conventions

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
| ``c_R``   | Riemann sign          | ``1``   |
| ``c_l``   | ``l``-leg scale       | ``1``   |
| ``c_m``   | ``m``-leg scale       | ``1``   |
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
unit-modulus complex number for an exotic convention.  ``c_α`` is the
sign in the time-transformation law ``t′ = κ(t − c_α α)`` (see [The BMS
group](@ref bms_group)); it is a convention like the others.

Every ``c``-parameter is defined in **export form**: it takes a
quantity *from* the SpEC convention *to* the convention ``X`` being
described, as in ``q^{[X]} = c_q q^{[\mathrm{SpEC}]}`` — and per tetrad
leg, ``ℓ^{[X]} = c_l ℓ^{[\mathrm{SpEC}]}`` and ``m^{[X]} = c_m
m^{[\mathrm{SpEC}]}``, with ``n^{[X]} = n^{[\mathrm{SpEC}]}/c_l``
forced by the ``ℓ·n`` normalization.  Because ``ℓ`` is a real null
vector, ``c_l`` must be real; and because ``m·m̄`` is fixed, ``c_m``
must have unit modulus.  The `Conventions` constructor validates both.
(Unit modulus is also what makes ``c_m^{-1} = c̄_m``, so factors like
the ``c_m^{-2}`` below are unambiguous.)

Two defaults deserve emphasis.  ``c_s=+1`` is the ``−+++`` signature
fixed above; note that the SpEC convention's ``c_s=+1`` already
*implies* ``−+++``.  And ``c_ð=1`` is the Newman–Penrose ``ð`` of [The
eth operator](@ref the-operator-eth), in which ``ð`` carries **no**
``\sqrt2`` — ``ð f = -(∂_θ + i\cscθ\,∂_φ)f`` on a spin-0 ``f`` (a
Geroch–Held–Penrose ``ð`` would be ``c_ð = 1/\sqrt2``).

Inter-convention conversion of the Weyl and Faraday components (from
the appendix) is

```math
ψ_n^{[X]} = c_s c_Ψ c_R\,(c_l c_m)^{2-n}\, ψ_n^{[\mathrm{SpEC}]},
\qquad
φ_n^{[X]} = c_φ\,(c_l c_m)^{1-n}\, φ_n^{[\mathrm{SpEC}]},
```

while the radiative quantities convert as

```math
σ^{[X]} = c_σ c_l^{\,ℐ} c_m^{2}\, σ^{[\mathrm{SpEC}]},
\qquad
h^{[X]} = c_s c_h^{-1} c_m^{-2}\, h^{[\mathrm{SpEC}]},
\qquad
N^{[X]} = c_s c_h^{-1} c_m^{-2}\, N^{[\mathrm{SpEC}]}.
```

Note that there is no ``c_s`` in the shear factor (the two raised ``m``
indices in ``σ = -c_σ mᵃmᵇ∇_a l_b`` contribute ``c_s^2 = 1``, and the
defining relations fix the one-form ``l_b`` directly), and the ``c_l``
inverts on ``ℐ⁻``, where the radiative shear is built from the ``n``
leg; see ["BMS action on fields"](@ref convention_dependence_fields)
for the derivations.  The News inherits the strain factor exactly,
because ``N = ∂_u h`` and the coordinates — including ``u`` — are
shared by all conventions; only the tetrad and field *definitions*
differ.  In code, these factors are `Scri.weyl_factor`,
`Scri.faraday_factor`, `Scri.shear_factor`, `Scri.strain_factor`, and
`Scri.news_factor` (dispatched per component by
`Scri.conversion_factor`), and `Scri.represent!` converts a data array
between two conventions by multiplying each component slice by the
appropriate ratio of factors.

A *transformation* that stays within one convention is more economical
than these inter-convention factors suggest: as worked out in ["BMS
action on fields"](@ref convention_dependence_fields), the
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
t'(t, 𝐧) = t - c_α α(𝐧),
```

where ``c_α = +1`` for a supertranslation of future null infinity and
``c_α = -1`` for a supertranslation of past null infinity.  The null
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
κ(𝐧) = \frac{1}{γ(1 - ℐ\, v⃗ ⋅ 𝐧)},
```

where ``𝐧`` here is the null ray *before* any boost or rotation has
been applied, and ``ℐ = ±1`` selects future or past null infinity (see
[Future and past null infinity](@ref scri_pm_conventions) below).  The full BMS action of
the element ``(Λ = b∘R, α)`` on ``A``'s coordinates ``(t, 𝐧)`` is then

```math
t' = κ(Λ, 𝐧) \bigl[t - c_α α(𝐧)\bigr],
\qquad
𝐧' = \frac{Λ𝐧}{(Λ𝐧)^0},
```

where ``κ(Λ, 𝐧) = n⁰/(Λ𝐧)⁰`` is the conformal factor for the
combined Lorentz transformation ``Λ`` evaluated on the null ray, which
is given in $A$'s coordinate frame as ``𝐧 = (1, ℐ n̂)`` up to
normalization.

## [Future and past null infinity](@id scri_pm_conventions)

Every transformation above lives on a chosen piece of ``ℐ``, and most
expressions carry a sign ``ℐ`` that records *which* piece:

```math
ℐ = \begin{cases} +1 & \text{future null infinity } ℐ⁺ \text{ (outgoing radiation)} \\ -1 & \text{past null infinity } ℐ⁻ \text{ (incoming radiation)}. \end{cases}
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

Both choices are captured by writing the section as ``𝐧 = (1, ℐ
n̂)``, which is the single place the convention enters.  Propagating
it through the conformal factor gives

```math
κ(𝐧) = \frac{1}{γ(1 - ℐ\, v⃗ ⋅ n̂)}.
```

The sign flip is exactly the past-cone Doppler factor: an observer
boosted with velocity ``v⃗`` sees radiation arriving from ``n̂``
shifted by ``γ(1 + v⃗ ⋅ n̂)`` (blueshift toward the direction of
motion), which is ``κ^{-1}`` at ``ℐ = -1``.  This is the same sign
difference Penrose and Rindler note at their Eq. (1.3.5), and it is
the convention implemented throughout this package — in
[`aberration`](@ref Scri.aberration) (its `ℐ` argument), in the
conformal factor used by [`transform!`](@ref), in the ``ℐ`` type
parameter of [`DataComponents`](@ref Scri.DataComponents) (data live
on a fixed null infinity), and in the ``ℐ`` representation
parameter of [`BMS`](@ref) elements (see [Representations of the
supertranslation](@ref bms_representations)).  Because ``κ`` itself is a spin-0
(orientation-blind) quantity, the choice also governs the *handedness*
of the spin-weighted ``ð`` operator that appears in the
component-mixing law, which is why the Weyl/Faraday peeling tower runs
in the opposite direction at ``ℐ⁻`` (see [Computing ``ðt'/κ``](@ref
computing_eth_tprime_over_kappa) and [`DataComponents`](@ref
Scri.DataComponents)).
