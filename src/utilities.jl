"""
    impose_reality(αᵢₙ, ℓₘₐₓ, εᵅ)

Given a set of mode weights `αᵢₙ` of a spin-weight 0 field, return a new set of mode weights
that satisfy the reality condition ``α_{ℓ,-m} = (-1)^m ᾱ_{ℓ,m}``.  Simultaneously, pad the
output array with zeros up to `ℓₘₐₓ`.  The input `αᵢₙ` is expected to be ordered by
increasing `ℓ`, starting from 0, then by increasing `m` within each `ℓ`.  The output array
is ordered in the same way, and has length `(ℓₘₐₓ + 1)^2`.
"""
function impose_reality(αᵢₙ, ℓₘₐₓ, εᵅ)
    Nᵢₙ = length(αᵢₙ)
    Lᵢₙ = isqrt(Nᵢₙ)
    @assert Lᵢₙ^2 == Nᵢₙ "Input `αᵢₙ` has $Nᵢₙ elements, which is not a perfect square"
    if ℓₘₐₓ < Lᵢₙ - 1
        throw(ArgumentError("Input `ℓₘₐₓ` is too small to accommodate input `αᵢₙ`"))
    end
    α = zeros(eltype(αᵢₙ), (ℓₘₐₓ + 1)^2)
    for ℓ ∈ 0:(Lᵢₙ-1)
        # The m=0 modes are purely real, so we just take the real part.
        i₀ = ℓ * (ℓ + 1) + 1
        α[i₀] = εᵅ * real(αᵢₙ[i₀])
        for m ∈ 1:ℓ
            i₊ = ℓ * (ℓ + 1) + m + 1
            i₋ = ℓ * (ℓ + 1) - m + 1
            α[i₊] = εᵅ * (αᵢₙ[i₊] + (-1)^m * conj(αᵢₙ[i₋])) / 2
            α[i₋] = (-1)^m * conj(α[i₊])
        end
    end
    return α
end


"""
    compute_t′(t, αₚ, Rₚ, v⃗)

Compute the new time samples `t′` corresponding to the input time samples `t` after a BMS
transformation with supertranslation `αₚ` and boost velocity `v⃗`.  The `Rₚ` describe the
locations of the pixels.

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
function compute_t′(t, αₚ, Rₚ, v⃗)
    β = absvec(v⃗)
    γ = 1 / √(1 - β^2)
    vˣ, vʸ, vᶻ = vec(v⃗)
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    T = promote_type(eltype(t), eltype(αₚ), typeof(γ))
    t′ₘᵢₙ, t′ₘₐₓ = typemin(T), typemax(T)
    @inbounds @simd for p ∈ eachindex(αₚ, Rₚ)
        (Rʷ, Rˣ, Rʸ, Rᶻ) = components(Rₚ[p])
        v⃗dotn̂ = (
            2vˣ * (Rʷ * Rʸ + Rˣ * Rᶻ)
            + 2vʸ * (Rʸ * Rᶻ - Rʷ * Rˣ)
            + vᶻ * (Rʷ^2 + Rᶻ^2 - Rˣ^2 - Rʸ^2)
        )
        k⁻¹ = γ * (1 - v⃗dotn̂)
        t′ₘᵢₙ = max(t′ₘᵢₙ, (tₘᵢₙ - αₚ[p]) / k⁻¹)
        t′ₘₐₓ = min(t′ₘₐₓ, (tₘₐₓ - αₚ[p]) / k⁻¹)
    end
    if t′ₘₐₓ ≤ t′ₘᵢₙ
        error(
            "\n\tThere are no complete slices in the t′ coordinate system "
            * "for t ∈ [$(tₘᵢₙ), ..., $(tₘₐₓ)], β = $β and α as given."
            * "\n\tYou may wish to decrease β or move the origin (zero) of "
            * "the time coordinate closer to the average value of t."
        )
    end
    scale = (t′ₘₐₓ - t′ₘᵢₙ) / (tₘₐₓ - tₘᵢₙ)
    t′ = @. t′ₘᵢₙ + scale * (t - tₘᵢₙ)
    t′[end] = t′ₘₐₓ  # ensure exact endpoint to avoid extrapolation
    tᵪ = (t′ₘᵢₙ - scale * tₘᵢₙ) / (1 - scale)
    return t′, tᵪ
end
