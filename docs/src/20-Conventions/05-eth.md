# [The Operator ``ð``](@id the-operator-eth)

Throughout this package, the operator[^1] ``ð`` is the spin-raising
**Newman–Penrose** operator, *not* the Geroch–Held–Penrose (GHP) one.
Restricted to the unit round sphere (with the boost weight dropping
out), the two differ by a factor of ``\sqrt{2}``:

```math
ð_{\mathrm{NP}} = \sqrt{2}\, ð_{\mathrm{GHP}}.
```

We choose the Newman–Penrose version as our default.  The distinction
should not be visible to users, but anyone interested in derivations
elsewhere in this documentation or in the implementation should be
aware of it.

Spin-weighted spherical functions [*cannot actually be defined*](@cite
Boyle_2016) over the sphere ``𝕊²`` itself; they are more properly
defined over the sphere ``𝕊³`` — or better yet ``\mathrm{Spin}(3)
\mathrel{\cong_{\text{Grp}}} \mathrm{SU}(2)``, the group of unit
quaternions.  Then we can easily define ``\eth`` in those terms.
Start with the right-derivative operator with respect to the generator
``𝔤`` acting on a function ``f`` and evaluated at ``Q``:

```math
R_𝔤 f(Q) = -i \left.\frac{d}{dϵ}\right|_{ϵ=0} f\left(Q e^{-ϵ𝔤/2}\right).
```

The ``ð`` operator is defined as

```math
ð = c_ð \left(R_{x} + i R_{y}\right),
```

where ``x`` and ``y`` generate rotations about their corresponding
axes.  This is the form we will use when [deriving ``ðt'/2κ.``](@ref
computing_eth_tprime_over_2kappa)

Of course, this may look unfamiliar.  The more common approach in the
literature takes advantage of the fact that *sometimes* we can *just
about* get away with writing spin-weighted spherical functions as
functions on *coordinates over the sphere*.  (That approach becomes
meaningless once we perform any transformation.)  In this way, Newman
and Penrose originally defined ð as acting on a quantity ``{}_s f`` of
spin weight ``s`` via

```math
ð\, {}_s f = -c_ð\, (\sin θ)^{s}\left(∂_θ + \frac{i}{\sin θ}∂_ϕ\right)
\left[(\sin θ)^{-s}\, {}_s f\right],
```

raising the spin weight by one.  This simpler definition is actually
implied by the more general and precise definition above, if we
parameterize the quaternions by spherical coordinates, which
conventionally represent an initial rotation about the ``y`` axis by
``θ`` followed by a rotation about the ``z`` axis by ``ϕ``.  Then, we
can write ``Q = e^{ϕ 𝐤/2} e^{θ𝐣/2}``, and the right-derivative
operator simplifies to precisely the Newman–Penrose form.

For a spin-0 function this is just ``ð f = -c_ð\left(∂_θ +
\frac{i}{\sin θ}∂_ϕ\right) f``, which ties ``ð`` to the angular dyad
``m`` of [the standard tetrad](@ref "Tetrad").  Since ``m̃ =
\frac{c_m}{\sqrt{2}}\left(∂_θ + \frac{i}{\sin θ}∂_ϕ\right)``,

```math
ð f = -c_ð\frac{\sqrt{2}}{c_m}\, m̃(f) \qquad (s = 0).
```

This relation and the tetrad normalizations are what fix the factors
of ``\sqrt{2}`` — and ultimately the ``1/2`` in the [Weyl mixing
parameter](@ref "BMS Action on Fields") — whenever ``ð`` acting on a
coordinate function is re-expressed through the tetrad.

It is helpful to note that, acting on the spin-weighted spherical
harmonics,

```math
ð\, {}_s Y_{ℓ,m} = c_ð\, \sqrt{(ℓ-s)(ℓ+s+1)}\; {}_{s+1} Y_{ℓ,m}.
```

The spin-weighted spherical harmonics can be more precisely defined on
``\mathrm{Spin}(3)`` [Boyle_2016](@cite), but again we obtain the
standard spherical-coordinate form is we restrict to those coordinates
as above.

[^1]: This character ``ð`` is the lowercase "eth", which looks like a
    partial derivative with a diagonal slash (not a horizontal cross)
    on its ascender.  It represents the *voiced* dental fricative — so
    the "th" sounds like the one in "this" or "that".  You should feel
    your vocal cords vibrate when you pronounce it.  This is as
    opposed to the "thorn" character ``þ``, which is used to represent
    a different derivative operator [GHP_1973](@cite) and represents
    the *unvoiced* dental fricative, as in "thin" or "thick", which
    should feel almost the same, but your vocal cords should not
    vibrate.
