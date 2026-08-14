"""
    impose_reality(αᵢₙ, ℓₘₐₓ, c_α)

Given a set of mode weights `αᵢₙ` of a spin-weight 0 field, return a new set of mode weights
that satisfy the reality condition ``α_{ℓ,-m} = (-1)^m ᾱ_{ℓ,m}``.  Simultaneously, pad the
output array with zeros up to `ℓₘₐₓ`.  The input `αᵢₙ` is expected to be ordered by
increasing `ℓ`, starting from 0, then by increasing `m` within each `ℓ`.  The output array
is ordered in the same way, and has length `(ℓₘₐₓ + 1)^2`.

The final argument `c_α` multiplies the symmetrized modes, folding the [supertranslation
sign convention](@ref conventions-overview) into the output so that downstream code can
uniformly interpret the result through ``t′ = κ(t - α)``.  Pass `1` when no sign
adjustment is wanted.
"""
function impose_reality(αᵢₙ, ℓₘₐₓ, c_α)
    Nᵢₙ = length(αᵢₙ)
    Lᵢₙ = isqrt(Nᵢₙ)
    if Lᵢₙ^2 != Nᵢₙ
        throw(ArgumentError("Input `αᵢₙ` has $Nᵢₙ elements, which is not a perfect square"))
    end
    if ℓₘₐₓ < Lᵢₙ - 1
        throw(ArgumentError("Input `ℓₘₐₓ` is too small to accommodate input `αᵢₙ`"))
    end
    α = zeros(eltype(αᵢₙ), (ℓₘₐₓ + 1)^2)
    for ℓ ∈ 0:(Lᵢₙ - 1)
        # The m=0 modes are purely real, so we just take the real part.
        i₀ = ℓ * (ℓ + 1) + 1
        α[i₀] = c_α * real(αᵢₙ[i₀])
        for m ∈ 1:ℓ
            i₊ = ℓ * (ℓ + 1) + m + 1
            i₋ = ℓ * (ℓ + 1) - m + 1
            α[i₊] = c_α * (αᵢₙ[i₊] + (-1)^m * conj(αᵢₙ[i₋])) / 2
            α[i₋] = (-1)^m * conj(α[i₊])
        end
    end
    return α
end

"""
    compute_t′(t, αₚ, Rₚ, v⃗, ℐ=1)

Compute the new time samples `t′` corresponding to the input time samples `t` after a BMS
transformation with supertranslation `αₚ` and boost velocity `v⃗`.  The `Rₚ` describe the
locations of the pixels.  The null-infinity sign `ℐ = ±1` enters the conformal factor as
`1/κ = γ(1 - ℐ v⃗⋅k̂)` (`+1` for ``ℐ⁺``, `-1` for ``ℐ⁻``).

The objective is to create a new time grid that has the same number of samples as `t` and
has roughly the same spacing, while accounting for the fact that some parts of the cylinder
(the subset of ℐ⁺) on which we have data will need to be dropped because the new slices at
the ends will be "tilted", and thus won't have a complete sphere of data on which to compute
modes.

The result is a simple rescaling:

    scale = (t′ₘₐₓ - t′ₘᵢₙ) / (tₘₐₓ - tₘᵢₙ)
    t′ = @. t′ₘᵢₙ + scale * (t - tₘᵢₙ)

For convenience, this function also returns `tᵪ`, the crossover time at which `t′=t` (the
fixed point of the t ↦ t′ map), which we can derive from the above formula as

    tᵪ = (t′ₘᵢₙ - scale * tₘᵢₙ) / (1 - scale)
"""
function compute_t′(t, αₚ, Rₚ, v⃗, ℐ=One())
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    t′ₘᵢₙ, t′ₘₐₓ = compute_t′_bounds(t, αₚ, Rₚ, v⃗, ℐ)
    if t′ₘₐₓ ≤ t′ₘᵢₙ
        error(
            "\n\tThere are no complete slices in the t′ coordinate system " *
            "for t ∈ [$(tₘᵢₙ), ..., $(tₘₐₓ)], β = $(absvec(v⃗)) and α as given." *
            "\n\tYou may wish to decrease β or move the origin (zero) of " *
            "the time coordinate closer to the average value of t.",
        )
    end
    return rescale_t′(t, t′ₘᵢₙ, t′ₘₐₓ)
end

"""
    rescale_t′(t, t′ₘᵢₙ, t′ₘₐₓ)

Map the input time samples `t` affinely onto the span `[t′ₘᵢₙ, t′ₘₐₓ]`, returning the new
grid `t′` together with the crossover time `tᵪ` at which `t′ = t` (the fixed point of the
affine map).  This is the common grid-construction step of the [`compute_t′`](@ref)
methods.
"""
function rescale_t′(t, t′ₘᵢₙ, t′ₘₐₓ)
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    scale = (t′ₘₐₓ - t′ₘᵢₙ) / (tₘₐₓ - tₘᵢₙ)
    t′ = @. t′ₘᵢₙ + scale * (t - tₘᵢₙ)
    t′[end] = t′ₘₐₓ  # ensure exact endpoint to avoid extrapolation
    tᵪ = (t′ₘᵢₙ - scale * tₘᵢₙ) / (1 - scale)
    return t′, tᵪ
end

"""
    compute_t′(t, β, δt)

Construct a single output time grid `t′` that is guaranteed to consist of complete slices
for *every* BMS transformation whose boost speed is at most `β` and whose supertranslation
satisfies `|α(k̂)| ≤ δt` in every direction `k̂` — regardless of the boost direction, the
frame rotation, and the choice of `ℐ = ±1`.

This is the natural choice for the `t′` keyword of [`transform!`](@ref) in an optimization
loop over BMS parameters, where the grid must remain valid at every parameter value the
optimizer tries, and where differentiation with
[ForwardDiff](https://juliadiff.org/ForwardDiff.jl/) requires the grid to be held fixed.
Estimate an upper bound `β` on the boost speed and an upper bound `δt` on the magnitude of
the supertranslation (including any time translation) that the optimizer may reach, and use
the resulting grid for every evaluation.

The grid is built exactly as in the exact-transformation method — an affine rescaling of `t`
onto the worst-case span of [`compute_t′_bounds`](@ref) — and the return value is likewise
`(t′, tᵪ)`, with `tᵪ` the fixed point of the affine map.  Because the span is a worst case,
it is generally smaller than the exact span for any particular transformation.
"""
function compute_t′(t, β::Real, δt::Real)
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    t′ₘᵢₙ, t′ₘₐₓ = compute_t′_bounds(t, β, δt)
    if t′ₘₐₓ ≤ t′ₘᵢₙ
        error(
            "\n\tThere are no complete slices in the t′ coordinate system " *
            "for t ∈ [$(tₘᵢₙ), ..., $(tₘₐₓ)] in the worst case allowed by " *
            "β = $β and δt = $δt." *
            "\n\tYou may wish to decrease β or δt, or move the origin (zero) of " *
            "the time coordinate closer to the average value of t.",
        )
    end
    return rescale_t′(t, t′ₘᵢₙ, t′ₘₐₓ)
end

@doc raw"""
    v⃗dotk̂(Rₚᵢ, vˣ, vʸ, vᶻ)

Compute ``v⃗ ⋅ k̂`` for the direction ``k̂ = Rₚᵢ 𝐤 R̃ₚᵢ`` singled out by the pixel rotor
`Rₚᵢ` — the third column of the rotation matrix of `Rₚᵢ`, contracted with the boost
velocity components.  This is the angular dependence entering the inverse conformal factor
``1/κ = γ(1 - ℐ v⃗⋅k̂)``.  Generic over the component types, so `Rₚᵢ` and the velocity may
carry dual numbers.
"""
@inline function v⃗dotk̂(Rₚᵢ, vˣ, vʸ, vᶻ)
    (Rʷ, Rˣ, Rʸ, Rᶻ) = components(Rₚᵢ)
    return (
        2vˣ * (Rʷ * Rʸ + Rˣ * Rᶻ) +
        2vʸ * (Rʸ * Rᶻ - Rʷ * Rˣ) +
        vᶻ * (Rʷ^2 + Rᶻ^2 - Rˣ^2 - Rʸ^2)
    )
end

"""
    compute_t′_bounds(t, αₚ, Rₚ, v⃗, ℐ=1)

Return the extreme values `(t′ₘᵢₙ, t′ₘₐₓ)` that a transformed time grid may take while every
pixel of the sphere still maps into the span of the input time samples `t` — that is, such
that `t′ κ⁻¹(k̂) + α(k̂) ∈ [t[begin], t[end]]` for every pixel direction.  The inputs are as
for [`compute_t′`](@ref), which uses these bounds to build its default grid; the same bounds
validate a user-supplied grid in [`transform!`](@ref).  A returned pair with `t′ₘₐₓ ≤ t′ₘᵢₙ`
means no complete slice exists.
"""
function compute_t′_bounds(t, αₚ, Rₚ, v⃗, ℐ=One())
    β = absvec(v⃗)
    γ = 1 / √(1 - β^2)
    vˣ, vʸ, vᶻ = vec(v⃗)
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    T = promote_type(eltype(t), eltype(αₚ), typeof(γ))
    t′ₘᵢₙ, t′ₘₐₓ = typemin(T), typemax(T)
    @inbounds @simd for p ∈ eachindex(αₚ, Rₚ)
        κ⁻¹ = γ * (1 - ℐ * v⃗dotk̂(Rₚ[p], vˣ, vʸ, vᶻ))
        t′ₘᵢₙ = max(t′ₘᵢₙ, (tₘᵢₙ - αₚ[p]) / κ⁻¹)
        t′ₘₐₓ = min(t′ₘₐₓ, (tₘₐₓ - αₚ[p]) / κ⁻¹)
    end
    return t′ₘᵢₙ, t′ₘₐₓ
end

"""
    compute_t′_bounds(t, β, δt)

Return the extreme values `(t′ₘᵢₙ, t′ₘₐₓ)` of a transformed time grid that remain valid for
*every* BMS transformation whose boost speed is at most `β` and whose supertranslation
satisfies `|α(k̂)| ≤ δt` in every direction `k̂`.  This is the worst-case counterpart of the
exact-transformation method above, for use when the transformation is not yet known — most
commonly in an optimization loop over BMS parameters, where a single fixed `t′` grid must
stay within the valid span at every parameter value the optimizer tries.

The bounds follow from the per-direction constraint `t′/κ(k̂) + α(k̂) ∈ [t[begin], t[end]]`:
the conformal factor `1/κ = γ(1 - ℐ v⃗⋅k̂)` ranges over `[γ(1-β), γ(1+β)]` as the direction
varies, and the supertranslation over `[-δt, δt]`, so taking the adverse extreme of each
gives a span contained in the exact span of any admissible transformation.  The boost
direction, the frame rotation, and the choice of `ℐ = ±1` only move these extremes around
the sphere, so none of them affects the result.  A returned pair with `t′ₘₐₓ ≤ t′ₘᵢₙ` means
that some admissible transformation may leave no complete slice.
"""
function compute_t′_bounds(t, β::Real, δt::Real)
    0 ≤ β < 1 || throw(ArgumentError("Boost-speed bound must satisfy 0 ≤ β < 1; got β=$β"))
    δt ≥ 0 || throw(ArgumentError("Supertranslation bound must satisfy δt ≥ 0; got δt=$δt"))
    γ = 1 / √(1 - β^2)
    κ⁻¹₋, κ⁻¹₊ = γ * (1 - β), γ * (1 + β)
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    t′ₘᵢₙ = max((tₘᵢₙ + δt) / κ⁻¹₋, (tₘᵢₙ + δt) / κ⁻¹₊)
    t′ₘₐₓ = min((tₘₐₓ - δt) / κ⁻¹₋, (tₘₐₓ - δt) / κ⁻¹₊)
    return t′ₘᵢₙ, t′ₘₐₓ
end

"""
    validate_t′(t′, t, αₚ, Rₚ, v⃗, ℐ=1)

Check that a user-supplied output time grid `t′` is usable by [`transform!`](@ref) or
[`transform_objective`](@ref): it must be nonempty, strictly increasing, and lie within
the pixel-wise bounds of [`compute_t′_bounds`](@ref) (up to a small roundoff allowance),
so that no pixel ever requires extrapolation beyond the input time span.  Throws an
`ArgumentError` otherwise; returns `t′` for convenience.

Note that this function does *not* require `length(t′) == length(t)`; that restriction
applies only to [`transform!`](@ref) — which writes its result back into `data` in place —
and is enforced at its call site.  [`transform_objective`](@ref) and
[`pixel_waveform`](@ref) accept output grids of any length.
"""
function validate_t′(t′, t, αₚ, Rₚ, v⃗, ℐ=One())
    isempty(t′) && throw(ArgumentError("Supplied `t′` must have at least one sample"))
    issorted(t′; lt=(≤)) ||
        throw(ArgumentError("Supplied `t′` must be strictly increasing"))
    t′ₘᵢₙ, t′ₘₐₓ = compute_t′_bounds(t, αₚ, Rₚ, v⃗, ℐ)
    if t′ₘₐₓ ≤ t′ₘᵢₙ
        throw(
            ArgumentError(
                "There are no complete slices in the t′ coordinate system for the given " *
                "inputs; no `t′` grid can be valid.  You may wish to decrease β or move " *
                "the origin (zero) of the time coordinate closer to the average value of t.",
            ),
        )
    end
    # Allow a tiny roundoff excursion past the exact bounds; the spline evaluation handles
    # a distance-O(ϵ) extrapolation gracefully (error O(τ³)).
    tol = 8 * eps(typeof(t′ₘₐₓ)) * max(abs(t′ₘᵢₙ), abs(t′ₘₐₓ), t′ₘₐₓ - t′ₘᵢₙ)
    if t′[begin] < t′ₘᵢₙ - tol || t′[end] > t′ₘₐₓ + tol
        throw(
            ArgumentError(
                "Supplied `t′` ∈ [$(t′[begin]), ..., $(t′[end])] exceeds the valid range " *
                "[$t′ₘᵢₙ, ..., $t′ₘₐₓ] for which every pixel maps into the input time span",
            ),
        )
    end
    return t′
end

@doc raw"""
    compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, ℐ=1)

Compute the null-rotation parameter ``ðt'/2κ`` (where ``t`` can represent either ``u`` or
``v``, depending on `ℐ`) on the boosted grid `Rₚ`, returned as a `2 × length(Rₚ)` matrix
split into a part independent of the time `t` (row 1) and the coefficient of `t` (row 2):

    ðt′╱2κ(t)[i] = ðt′╱2κ[1, i] + t * ðt′╱2κ[2, i].

The time-coefficient is just ``ðκ/2κ``, which has the closed form below in terms of
``λ = R̃ᵢ v⃗ Rᵢ`` (the rest-frame components of the boost velocity at pixel `i`); the
time-independent part folds in the ``c_α``-corrected supertranslation `αₚ` and its
eth-derivative `ðαₚ`:

    ðt′╱2κ[2, i] = (λˣ + im * λʸ) / 2(λᶻ - ℐ),
    ðt′╱2κ[1, i] = -(ðt′╱2κ[2, i] * αₚ[i] + ðαₚ[i] / 2).

See the documentation page ["Computing ``ðt'/κ``"](@ref computing_eth_tprime_over_2kappa) for
the derivation.  Note in particular that the boost × supertranslation cross term ``ðt'╱2κ[2,
i]·αₚ[i]`` enters with a **minus** sign.
"""
function compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, ℐ=One())
    vˣ, vʸ, vᶻ = vec(v⃗)
    T = promote_type(basetype(eltype(Rₚ)), eltype(αₚ), real(eltype(ðαₚ)), typeof(vˣ))
    ðt′╱2κ = Matrix{Complex{T}}(undef, 2, length(Rₚ))
    @inbounds @simd ivdep for i ∈ eachindex(Rₚ, αₚ, ðαₚ)
        Rₚᵢʷ, Rₚᵢˣ, Rₚᵢʸ, Rₚᵢᶻ = components(Rₚ[i])
        λˣ = (
            (Rₚᵢʷ^2 + Rₚᵢˣ^2 - Rₚᵢʸ^2 - Rₚᵢᶻ^2)*vˣ +
            (-Rₚᵢʷ*Rₚᵢʸ + Rₚᵢˣ*Rₚᵢᶻ)*2vᶻ +
            (Rₚᵢˣ*Rₚᵢʸ + Rₚᵢʷ*Rₚᵢᶻ)*2vʸ
        )
        λʸ = (
            (Rₚᵢʷ^2 - Rₚᵢˣ^2 + Rₚᵢʸ^2 - Rₚᵢᶻ^2)*vʸ +
            (Rₚᵢˣ*Rₚᵢʸ - Rₚᵢʷ*Rₚᵢᶻ)*2vˣ +
            (Rₚᵢʸ*Rₚᵢᶻ + Rₚᵢʷ*Rₚᵢˣ)*2vᶻ
        )
        λᶻ = (
            (Rₚᵢʷ^2 + Rₚᵢᶻ^2 - Rₚᵢˣ^2 - Rₚᵢʸ^2)*vᶻ +
            (-Rₚᵢʷ*Rₚᵢˣ + Rₚᵢʸ*Rₚᵢᶻ)*2vʸ +
            (Rₚᵢˣ*Rₚᵢᶻ + Rₚᵢʷ*Rₚᵢʸ)*2vˣ
        )
        ðκ╱2κ = (λˣ + im * λʸ) / 2(λᶻ - ℐ)
        ðt′╱2κ[2, i] = ðκ╱2κ
        ðt′╱2κ[1, i] = -ðκ╱2κ * αₚ[i] - ðαₚ[i] / 2
    end
    return ðt′╱2κ
end

"""
    diagnostics(data, data_components)

Compute the per-``ℓ`` power monitors ``E(ℓ) = Σₘ |f_{ℓ,m}|²`` for each component of the
input data, as functions of time.

`data` is a complex mode-weight array with dimensions `(Nᵐ, Nᵗ, Nᵈ)`, exactly as for
[`transform!`](@ref), and `data_components` is the corresponding [`DataComponents`](@ref)
descriptor.  Returns a `Dict{Symbol,Matrix}` mapping each component symbol (e.g., `:h` or
`:ψ₄`) to a real `Nᵗ × (ℓₘₐₓ+1)` matrix whose `[j, ℓ+1]` entry is ``E(ℓ)`` at time index
`j`.

For a component of spin weight `s`, the `ℓ < |s|` rows of `data` hold the null-space
diagnostic `ξ` of the [augmented SSHT](@ref "Augmented Direct SSHT") after a transform, so
the corresponding `E(ℓ)` values measure power that no band-limited field can represent;
together with `E(ℓₘₐₓ)`, these indicate whether a transformation was adequately resolved.
See [Choosing ``ℓ_\\mathrm{max}``](@ref) for how to interpret them.

Loading `Plots` activates an extension providing `diagnostics(t, data, data_components)`,
which returns a corresponding `Dict` of plots of these powers against time.
"""
function diagnostics(data, data_components::DataComponents{C,ℐ}) where {C,ℐ}
    Nᵐ, Nᵗ, Nᵈ = size(data)
    L = isqrt(Nᵐ)
    if L^2 != Nᵐ
        throw(ArgumentError("Input `data` has $Nᵐ modes, which is not a perfect square"))
    end
    ℓₘₐₓ = L - 1
    diag = Dict{Symbol,Matrix{real(eltype(data))}}()
    for (d, comp) ∈ enumerate(C)
        power = Matrix{real(eltype(data))}(undef, Nᵗ, ℓₘₐₓ+1)
        for ℓ ∈ 0:ℓₘₐₓ
            mode_indices = (ℓ ^ 2 + 1):((ℓ + 1) ^ 2)
            power[:, ℓ + 1] = sum(abs2, (@view data[mode_indices, :, d]); dims=1)[1, :]
        end
        diag[comp] = power
    end
    return diag
end

@doc raw"""
    rotate_modes(α, Q, ℓₘₐₓ)

Rigidly rotate the spin-weight-0 mode weights `α` by the rotor `Q`: if the input modes
represent the function ``f``, the output modes represent ``f′(k̂) = f(Q⁻¹ k̂)`` — that is,
the function actively rotated by `Q`.  In terms of Wigner's ``𝔇`` matrices,

```math
f′_{ℓ,m′} = \sum_m \overline{𝔇^ℓ_{m′,m}(Q)}\, f_{ℓ,m}.
```

A rigid rotation preserves the band limit, so this is exact (up to roundoff) at the input
resolution; `ℓₘₐₓ` sets the output resolution, zero-padding above the input's band limit.
The mode ordering is as in [`BMS`](@ref).
"""
function rotate_modes(α::Vector{Complex{T}}, Q::Rotor, ℓₘₐₓ::Int) where {T}
    ℓᵅ = isqrt(length(α)) - 1
    D = D_matrices(Q, min(ℓᵅ, ℓₘₐₓ))
    out = zeros(Complex{T}, (ℓₘₐₓ + 1)^2)
    for ℓ ∈ 0:min(ℓᵅ, ℓₘₐₓ), m′ ∈ (-ℓ):ℓ
        s = zero(Complex{T})
        for m ∈ (-ℓ):ℓ
            s += conj(D[WignerDindex(ℓ, m′, m)]) * α[ℓ ^ 2 + ℓ + m + 1]
        end
        out[ℓ ^ 2 + ℓ + m′ + 1] = s
    end
    return out
end

"""
    primal_float(T)

The floating-point type underlying `T`, used to build *parameter-independent* objects —
the output pixel grid and its analysis factorizations in [`transform!`](@ref).  For
ordinary real types this is `T` itself.  The ForwardDiff package extension peels
dual-number types down to their value type, so that automatic differentiation does not
propagate (identically zero) perturbations into `qr`/`lu` factorizations of constant
matrices — which would poison every derivative with NaNs.
"""
primal_float(::Type{T}) where {T<:Real} = T
