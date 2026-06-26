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
    for ℓ ∈ 0:(Lᵢₙ - 1)
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
    compute_t′(t, αₚ, Rₚ, v⃗, εᴵ=1)

Compute the new time samples `t′` corresponding to the input time samples `t` after a BMS
transformation with supertranslation `αₚ` and boost velocity `v⃗`.  The `Rₚ` describe the
locations of the pixels.  The null-infinity sign `εᴵ = ±1` enters the conformal factor as
`1/κ = γ(1 - εᴵ v⃗⋅n̂)` (`+1` for ``ℐ⁺``, `-1` for ``ℐ⁻``).

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
function compute_t′(t, αₚ, Rₚ, v⃗, εᴵ=1)
    β = absvec(v⃗)
    γ = 1 / √(1 - β^2)
    vˣ, vʸ, vᶻ = vec(v⃗)
    tₘᵢₙ, tₘₐₓ = t[begin], t[end]
    T = promote_type(eltype(t), eltype(αₚ), typeof(γ))
    t′ₘᵢₙ, t′ₘₐₓ = typemin(T), typemax(T)
    @inbounds @simd for p ∈ eachindex(αₚ, Rₚ)
        (Rʷ, Rˣ, Rʸ, Rᶻ) = components(Rₚ[p])
        v⃗dotn̂ = (
            2vˣ * (Rʷ * Rʸ + Rˣ * Rᶻ) +
            2vʸ * (Rʸ * Rᶻ - Rʷ * Rˣ) +
            vᶻ * (Rʷ^2 + Rᶻ^2 - Rˣ^2 - Rʸ^2)
        )
        κ⁻¹ = γ * (1 - εᴵ * v⃗dotn̂)
        t′ₘᵢₙ = max(t′ₘᵢₙ, (tₘᵢₙ - αₚ[p]) / κ⁻¹)
        t′ₘₐₓ = min(t′ₘₐₓ, (tₘₐₓ - αₚ[p]) / κ⁻¹)
    end
    if t′ₘₐₓ ≤ t′ₘᵢₙ
        error(
            "\n\tThere are no complete slices in the t′ coordinate system " *
            "for t ∈ [$(tₘᵢₙ), ..., $(tₘₐₓ)], β = $β and α as given." *
            "\n\tYou may wish to decrease β or move the origin (zero) of " *
            "the time coordinate closer to the average value of t.",
        )
    end
    scale = (t′ₘₐₓ - t′ₘᵢₙ) / (tₘₐₓ - tₘᵢₙ)
    t′ = @. t′ₘᵢₙ + scale * (t - tₘᵢₙ)
    t′[end] = t′ₘₐₓ  # ensure exact endpoint to avoid extrapolation
    tᵪ = (t′ₘᵢₙ - scale * tₘᵢₙ) / (1 - scale)
    return t′, tᵪ
end

@doc raw"""
    compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, εᴵ=1)

Compute the null-rotation parameter ``ðt'/2κ`` (where ``t`` can represent either ``u`` or
``v``, depending on `εᴵ`) on the boosted grid `Rₚ`, returned as a `2 × length(Rₚ)` matrix
split into a part independent of the time `t` (row 1) and the coefficient of `t` (row 2):

    ðt′╱2κ(t)[i] = ðt′╱2κ[1, i] + t * ðt′╱2κ[2, i].

The time-coefficient is just ``ðκ/2κ``, which has the closed form below in terms of
``λ = R̃ᵢ v⃗ Rᵢ`` (the rest-frame components of the boost velocity at pixel `i`); the
time-independent part folds in the ``εᵅ``-corrected supertranslation `αₚ` and its
eth-derivative `ðαₚ`:

    ðt′╱2κ[2, i] = (λˣ + im * λʸ) / 2(λᶻ - εᴵ),
    ðt′╱2κ[1, i] = -(ðt′╱2κ[2, i] * αₚ[i] + ðαₚ[i] / 2).

See the documentation page ["Computing ``ðt'/κ``"](@ref computing_eth_tprime_over_kappa) for
the derivation.  Note in particular that the boost × supertranslation cross term ``ðt'╱2κ[2,
i]·αₚ[i]`` enters with a **minus** sign.
"""
function compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, εᴵ=1)
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
        ðt′╱2κ[2, i] = (λˣ + im * λʸ) / 2(λᶻ - εᴵ)
        ðt′╱2κ[1, i] = -(ðt′╱2κ[2, i] * αₚ[i] + ðαₚ[i] / 2)
    end
    return ðt′╱2κ
end

"""
    diagnostics(data, data_components)

Compute power monitors for the input data.
"""
function diagnostics(data, data_components::DataComponents{C,εᴵ}) where {C,εᴵ}
    Nᵐ, Nᵗ, Nᵈ = size(data)
    L = isqrt(Nᵐ)
    @assert L^2 == Nᵐ "Input `data` has $Nᵐ modes, which is not a perfect square"
    ℓₘₐₓ = L - 1
    diag = Dict{Symbol,Matrix{real(eltype(data))}}()
    for (d, comp) ∈ enumerate(C)
        power = Matrix{real(eltype(data))}(undef, Nᵗ, ℓₘₐₓ+1)
        for ℓ ∈ 0:ℓₘₐₓ
            mode_indices = (ℓ ^ 2 + 1):((ℓ + 1) ^ 2)
            power[:, ℓ + 1] = sum(abs2, (@view data[mode_indices, :, d]); dims=1)[1, :]
        end
        diag[comp] = power
        # for ℓ ∈ 0:ℓₘₐₓ
        #     mode_indices = ℓ^2+1:(ℓ+1)^2
        #     power = sum(abs2, (@view data[mode_indices, :, d]), dims=1)[1, :]
        #     diag[comp] = power
        # end
    end
    return diag
end
