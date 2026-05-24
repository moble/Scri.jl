# Conventions

An important source of these conventions is [BoyleEtAl_2019](@citet),
which describes conventions used for SXS waveforms in Appendix C.  BMS
conventions are described in [MitmanEtAl_2024](@citet).  The
conventions for the geometric algebra

!!! warn

    Note LISA Rosetta Stone when it comes out.

Nominal values for solar and planetary quantities are given in
[PršaEtAl_2016](@citet).  Other constants can be found in the
[`lisaconstants` package](https://pypi.org/project/lisaconstants/),
which otherwise integrates with [`astropy`](https://www.astropy.org/).

The signature is ``{-}{+}{+}{+}``, and we use units where ``G=c=1``.

The Christoffel symbols and the Riemann, Ricci, and Weyl tensors
follow the Misner-Thorne-Wheeler conventions
[MisnerThorneWheeler_1973](@cite) — Eqs. (14.36), (8.44), (8.47), and
(13.50) of that reference, respectively.  The Newman-Penrose Weyl
components are defined as

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
We have the asymptotic relation

```math
\Psi_4 \sim -\ddot{h},
```

where the dots indicate time derivatives.
