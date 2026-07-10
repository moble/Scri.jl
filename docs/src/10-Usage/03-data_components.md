```@meta
CurrentModule = Scri
```

# Data Components

The waveform array passed to [`transform!`](@ref) is three-dimensional,
with its third axis running over the different physical fields — strain,
Weyl components, Faraday components, and so on.  The
[`DataComponents`](@ref Scri.DataComponents) descriptor declares *which*
fields those are, and in *what order*, so that the transform knows the
spin weight, conformal weight, and mixing law of each slice.  It also
records which end of null infinity the data live on and which
[conventions](@ref conventions-type) they are expressed in.

## The components

Each component is named by a symbol; the valid names are the
Newman–Penrose Weyl components ``ψ₀, …, ψ₄``, the Maxwell/Faraday
components ``φ₀, φ₁, φ₂``, the shears ``σ`` (on ``ℐ⁺``) and ``λ`` (on
``ℐ⁻``), the complex strain ``h``, and the news ``News``.  Their spin
weights are those of [the definitions page](@ref "Tensor Components"):

| component              | spin weight ``s``     |
|:-----------------------|:---------------------:|
| ``ψ₀, ψ₁, ψ₂, ψ₃, ψ₄`` | ``+2, +1, 0, -1, -2`` |
| ``φ₀, φ₁, φ₂``         | ``+1, 0, -1``         |
| ``σ``                  | ``+2``                |
| ``λ``, ``h``, ``News`` | ``-2``                |

You may spell the names as strings in almost any reasonable form — the
constructor parses them into the canonical symbols:

```julia
DataComponents("Psi_3", "psi4", "sigma")   # → DataComponents(:ψ₃, :ψ₄, :σ)
```

## Future versus past null infinity

The `ℐ` sign selects the piece of null infinity the data describe: `+1`
for ``ℐ⁺`` (outgoing radiation) and `-1` for ``ℐ⁻`` (incoming), with
`+1` the default.  The two ends carry *different* radiative shears — the
``ℓ``-congruence ``σ`` on ``ℐ⁺`` and the ``n``-congruence ``λ`` on
``ℐ⁻`` (see [Tensor Components](@ref)) — so the constructor rejects
``σ`` on ``ℐ⁻`` and ``λ`` on ``ℐ⁺``:

```julia
DataComponents(:ψ₄, :σ)              # ℐ⁺ (default): shear σ
DataComponents(:ψ₀, :ψ₁, :λ; ℐ=-1)   # ℐ⁻: shear λ, and note ψ₀ leads the tower
```

The `ℐ` of a `DataComponents` is conceptually distinct from the `ℐ` of
a [`BMS`](@ref) element: one says which data you hold, the other says
which labeling a supertranslation is expressed in (see [Representations
of the supertranslation](@ref bms_representations)).

## The Weyl and Faraday hierarchy

Because the [peeling tower](@ref "Mixing: the peeling tower") mixes a
component into all the components below it, transforming any Weyl
component requires every component above it in the tower to be present.
On ``ℐ⁺`` the tower is fixed by ``ψ₄`` and fills downward, so ``ψᵢ``
requires all ``ψⱼ`` with ``j > i``; the same holds for the Faraday
components (``φ₀`` requires ``φ₁`` and ``φ₂``).  On ``ℐ⁻`` the tower is
fixed by ``ψ₀`` instead and runs the other way, so ``ψᵢ`` requires all
``ψⱼ`` with ``j < i``.  The constructor enforces these dependencies:

```julia
DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄)  # OK: full tower on ℐ⁺
DataComponents(:ψ₃, :ψ₄)                 # OK: ψ₃ needs only ψ₄ above it on ℐ⁺
DataComponents(:ψ₀, :ψ₁; ℐ=-1)           # OK: ψ₀ needs only ψ₁ above it on ℐ⁻
DataComponents(:ψ₃)                      # ERROR: ψ₃ requires ψ₄
```

The strain, news, and shears carry no such requirement — their laws
are [self-contained](@ref "Strain, shear, and news").

## Conventions travel with the data

The `conventions` field records the [convention](@ref conventions-type)
the component data are written in; it defaults to the package-native SXS
conventions.  The conventions ride along with the descriptor:
[`transform!`](@ref) applies the transformation laws *native to that
convention* (so no round-trip conversion is needed), while
[`represent!`](@ref Scri.represent!) re-expresses a data array from one
convention in another.

```julia
DataComponents(:ψ₄, :h; conventions=Conventions(:MB))
```

Because the components, the `ℐ` sign, and the conventions are all
encoded at the type level, the per-component spin weights, array
indices, and convention factors are compile-time constants — the branch
for each component specializes away, and default conventions compile to
the identity.

```@docs; canonical=false
Scri.DataComponents
```
