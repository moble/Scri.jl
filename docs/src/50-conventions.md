# Conventions

An important source of these conventions is [BoyleEtAl_2019](@citet),
which describes conventions used for SXS waveforms in Appendix C.  BMS
conventions are described in [MitmanEtAl_2024](@citet).  The
conventions for the geometric algebra are described in the
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
the conformal factor ``K``:

```math
t' = K(𝐧)\, t,
\qquad
K(𝐧) = \frac{1}{γ(1 - ε^ℐ\, v⃗ ⋅ 𝐧)},
```

where ``𝐧`` here is the null ray *before* any boost or rotation has
been applied, and ``ε^ℐ = ±1`` selects future or past null infinity (see
[Future and past null infinity](@ref scri_pm_conventions) below).  The full BMS action of
the element ``(Λ = b∘R, α)`` on ``A``'s coordinates ``(t, 𝐧)`` is then

```math
t' = K(Λ, 𝐧) \bigl[t - εᵅ α(𝐧)\bigr],
\qquad
𝐧' = \frac{Λ𝐧}{(Λ𝐧)^0},
```

where ``K(Λ, 𝐧) = n⁰/(Λ𝐧)⁰`` is the conformal factor for the
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
K(𝐧) = \frac{1}{γ(1 - ε^ℐ\, v⃗ ⋅ n̂)}.
```

The sign flip is exactly the past-cone Doppler factor: an observer
boosted with velocity ``v⃗`` sees radiation arriving from ``n̂``
shifted by ``γ(1 + v⃗ ⋅ n̂)`` (blueshift toward the direction of
motion), which is ``K^{-1}`` at ``ε^ℐ = -1``.  This is the same sign
difference Penrose and Rindler note at their Eq. (1.3.5), and it is
the convention implemented throughout this package — in
[`aberration`](@ref Scri.aberration) (its `emitted` keyword), in the
conformal factor used by [`transform!`](@ref), in the ``ε^ℐ`` type
parameter of [`DataComponents`](@ref Scri.DataComponents) (data live
on a fixed null infinity), and in the ``ε^ℐ`` argument of the
[`BMS`](@ref) operations.  Because ``K`` itself is a spin-0
(orientation-blind) quantity, the choice also governs the *handedness*
of the spin-weighted ``ð`` operator that appears in the
component-mixing law, which is why the Weyl/Faraday peeling tower runs
in the opposite direction at ``ℐ⁻`` (see [Computing ``ðt'/K``](@ref
computing_eth_tprime_over_K) and [`DataComponents`](@ref
Scri.DataComponents)).
