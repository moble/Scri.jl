# [The Operator ``ð``](@id the-operator-eth)

Throughout this package, the operator[^1] ``ð`` is the spin-raising
**Newman–Penrose** operator, *not* the Geroch–Held–Penrose (GHP) one.
Restricted to the unit round sphere (with the boost weight dropping
out), the two differ by a factor of ``\sqrt{2}``:

```math
ð_{\mathrm{NP}} = \sqrt{2}\, ð_{\mathrm{GHP}}.
```

While spin-weighted spherical functions [*cannot actually be
defined*](@cite Boyle_2016) on the sphere ``𝕊²`` itself, we can often
just about get away with writing them as functions on *coordinates
over the sphere*.  (That approach becomes meaningless once we perform
any transformation.)  This is the standard approach in the literature,
and as such the Newman–Penrose ð is defined as acting on a quantity
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

Equivalently the GHP ð is simply ``ð_{\mathrm{GHP}} f = -m̃(f)``.
This relation and the tetrad normalizations are what fix the factors
of ``\sqrt{2}`` — and ultimately the ``1/2`` in the [Weyl mixing
parameter](@ref "BMS action on fields") — whenever ``ð`` of a
coordinate is re-expressed through the tetrad.

[^1]: This character is the lowercase "eth", which looks like a
    partial derivative with a diagonal slash (not a horizontal cross)
    on its ascender.  It represents the *voiced* dental fricative — so
    the "th" sounds like the one in "this" or "that", as opposed to
    the one in "thin" or "thick".  You should feel your vocal cords
    vibrate when you pronounce it.
