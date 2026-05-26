const ValidDataComponents = (:Ψ₀, :Ψ₁, :Ψ₂, :Ψ₃, :Ψ₄, :σ, :h, :News, :φ₀, :φ₁, :φ₂)

"""
    DataComponents{C, εᴵ}

Encodes a fixed set of waveform data components at the type level.  `C` is an
`NTuple{N,Symbol}` whose elements are drawn from `ValidDataComponents`:

    $ValidDataComponents

The parameter `εᴵ` is the sign of the time direction for which these components are defined:
`+1` for ``ℐ⁺`` (outgoing) and `-1` for ``ℐ⁻`` (incoming).  The default value is `+1`, since
most applications will be for ``ℐ⁺``.

Examples:

    DataComponents(:Ψ₄)                          # gravitational waves only
    DataComponents(:Ψ₄, :Ψ₃, :Ψ₂)                # top three Weyl components
    DataComponents(:Ψ₀, :Ψ₁, :Ψ₂, :Ψ₃, :Ψ₄, :σ)  # full Weyl set with strain
    DataComponents(:φ₀, :φ₁, :φ₂)                # Faraday components
    DataComponents(:φ₀, :φ₁, :φ₂; εᴵ=-1)         # Faraday components on ℐ⁻
"""
struct DataComponents{C, Eᴵ}
    function DataComponents(cs::Symbol...; εᴵ=1)
        validate_data_components(cs, εᴵ)
        new{cs, εᴵ}()
    end
end


function validate_data_components(cs, εᴵ)
    @assert all(c -> c ∈ ValidDataComponents, cs) "" *
        "Invalid component in $cs; allowed: $ValidDataComponents"
    if εᴵ == 1
        :Ψ₀ ∈ cs && @assert :Ψ₁ ∈ cs "Ψ₀ requires Ψ₁ on ℐ⁺ (εᴵ=+1)"
        :Ψ₁ ∈ cs && @assert :Ψ₂ ∈ cs "Ψ₁ requires Ψ₂ on ℐ⁺ (εᴵ=+1)"
        :Ψ₂ ∈ cs && @assert :Ψ₃ ∈ cs "Ψ₂ requires Ψ₃ on ℐ⁺ (εᴵ=+1)"
        :Ψ₃ ∈ cs && @assert :Ψ₄ ∈ cs "Ψ₃ requires Ψ₄ on ℐ⁺ (εᴵ=+1)"
        :φ₀ ∈ cs && @assert :φ₁ ∈ cs "φ₀ requires φ₁ on ℐ⁺ (εᴵ=+1)"
        :φ₁ ∈ cs && @assert :φ₂ ∈ cs "φ₁ requires φ₂ on ℐ⁺ (εᴵ=+1)"
    elseif εᴵ == -1
        :Ψ₄ ∈ cs && @assert :Ψ₃ ∈ cs "Ψ₄ requires Ψ₃ on ℐ⁻ (εᴵ=-1)"
        :Ψ₃ ∈ cs && @assert :Ψ₂ ∈ cs "Ψ₃ requires Ψ₂ on ℐ⁻ (εᴵ=-1)"
        :Ψ₂ ∈ cs && @assert :Ψ₁ ∈ cs "Ψ₂ requires Ψ₁ on ℐ⁻ (εᴵ=-1)"
        :Ψ₁ ∈ cs && @assert :Ψ₀ ∈ cs "Ψ₁ requires Ψ₀ on ℐ⁻ (εᴵ=-1)"
        :φ₂ ∈ cs && @assert :φ₁ ∈ cs "φ₂ requires φ₁ on ℐ⁻ (εᴵ=-1)"
        :φ₁ ∈ cs && @assert :φ₀ ∈ cs "φ₁ requires φ₀ on ℐ⁻ (εᴵ=-1)"
    else
        throw(ArgumentError("Invalid εᴵ = $εᴵ; must be ±1"))
    end
end


"""
    component_index(dc::DataComponents{C}, ::Val{S})

Return the index of component `S` within `dc`, or `nothing` if absent.  The return type is
inferred as a compile-time constant when `dc` has a concrete type.
"""
@inline component_index(::DataComponents{C}, ::Val{S}) where {C, S} =
    findfirst(==(S), C)


"""
    has_component(dc::DataComponents{C}, ::Val{S})

Return `true` if component `S` is present in `dc`.
"""
@inline has_component(::DataComponents{C}, ::Val{S}) where {C, S} = S ∈ C


"""
    ncomponents(dc::DataComponents)

Return the number of components in `dc`.
"""
ncomponents(::DataComponents{C}) where C = length(C)


"""
    spin_weight(::Val{S})

Return the spin weight of field component `S`.
"""
spin_weight(::Val{:Ψ₀}) =  2
spin_weight(::Val{:Ψ₁}) =  1
spin_weight(::Val{:Ψ₂}) =  0
spin_weight(::Val{:Ψ₃}) = -1
spin_weight(::Val{:Ψ₄}) = -2
spin_weight(::Val{:σ})  =  2
spin_weight(::Val{:h})  = -2
spin_weight(::Val{:News})  = -2
spin_weight(::Val{:φ₀}) =  1
spin_weight(::Val{:φ₁}) =  0
spin_weight(::Val{:φ₂}) = -1
Base.@constprop :aggressive spin_weight(s)::Int = spin_weight(Val(s))

"""
    conformal_weight(::Val{S})

Return the conformal weight (spin weight + boost weight) of field component `S`, which is
the power of the conformal factor ``k`` in the BMS transformation law.

Note that we are assuming that these fields represent the asymptotic values of the physical
fields at null infinity, so they have already been rescaled by the appropriate power of the
conformal factor to be finite and nonzero at null infinity.  The transformation law for the
physical fields (finite-radius Weyl and Faraday spinors) would have a different conformal
weight, but the asymptotic fields are the ones we are transforming.

More specifically, the *asymptotic* Weyl spinor ``ψ`` and Faraday spinor ``φ`` are related
to the finite-radius Weyl spinor ``Ψ`` by ``ψ = ωΨ`` and the finite-radius Faraday spinor
``Φ`` by ``φ = ωΦ``.  The factor ``ω`` is the conformal factor that goes to zero (but has
nonzero derivative) at null infinity, which transforms as ``ω′=kω``, so we pick up a factor
of ``k⁻¹`` in the transformation laws for ``ψ`` and ``φ`` compared to ``Ψ`` and ``Φ``.
Since ``Ψ`` and ``Φ`` are the physical quantities, they do not change under coordinate
transformations.

Meanwhile the basis spinors each transform with a factor of ``1/√k``.  The Weyl components
``Ψₙ`` are defined by contracting the Weyl spinor with four basis spinors, so they pick up a
factor of ``k⁻²`` from the basis spinors and an additional factor of ``k⁻¹`` from the
conformal factor, for a total of ``k⁻³``.  Similarly, the Faraday components ``Φₙ`` are
defined by contracting the Faraday spinor with two basis spinors, so they pick up a factor
of ``k⁻¹`` from the conformal factor and an additional factor of ``k⁻¹`` from the basis
spinors, for a total of ``k⁻²``.
"""
conformal_weight(::Val{:Ψ₀}) = -3
conformal_weight(::Val{:Ψ₁}) = -3
conformal_weight(::Val{:Ψ₂}) = -3
conformal_weight(::Val{:Ψ₃}) = -3
conformal_weight(::Val{:Ψ₄}) = -3
conformal_weight(::Val{:σ})  = -1
conformal_weight(::Val{:h})  = -1
conformal_weight(::Val{:News})  = -2
conformal_weight(::Val{:φ₀}) = -2
conformal_weight(::Val{:φ₁}) = -2
conformal_weight(::Val{:φ₂}) = -2


"""
    mix_components!(dataᵢⱼ, k⁻¹, ðu′╱k, ð²α, dc)

Apply the BMS component-mixing transformation to `dataᵢⱼ`.  `k⁻¹` is the inverse conformal
factor for this pixel, `ðu′╱k` is the eth-derivative of the retarded time in the new frame
divided by ``k``, and `ð²α` is the sign-adjusted second anti-eth-derivative of the
supertranslation (used for the strain/shear component).

Note that Julia specializes on the concrete type of `dc`.  This means that the indexes into
`dataᵢⱼ` for the various components are known at compile time, and the branches for which
components are present or absent will be resolved at compile time.  The result is that —
even though the function body looks unwieldy and slow — it compiles down to minimal code
with no branches and only the necessary components, making it very fast in practice.  The
`@inline` annotation encourages this specialization and inlining, especially when just a few
components are being processed.
"""
@inline function mix_components!(
    dataᵢⱼ::AbstractVector{Complex{T}},
    k⁻¹,
    ðu′╱k,
    ð²α,
    dc::DataComponents{C, Eᴵ}
) where {T, C, Eᴵ}
    k⁻² = k⁻¹ * k⁻¹
    k⁻³ = k⁻² * k⁻¹
    ð̄²α = conj(ð²α)

    iΨ₄ = component_index(dc, Val(:Ψ₄))
    iΨ₃ = component_index(dc, Val(:Ψ₃))
    iΨ₂ = component_index(dc, Val(:Ψ₂))
    iΨ₁ = component_index(dc, Val(:Ψ₁))
    iΨ₀ = component_index(dc, Val(:Ψ₀))
    iσ  = component_index(dc, Val(:σ))
    ih  = component_index(dc, Val(:h))
    iNews = component_index(dc, Val(:News))
    iφ₂ = component_index(dc, Val(:φ₂))
    iφ₁ = component_index(dc, Val(:φ₁))
    iφ₀ = component_index(dc, Val(:φ₀))

    @inbounds begin
        ψ₀ = isnothing(iΨ₀) ? 0 : dataᵢⱼ[iΨ₀]
        ψ₁ = isnothing(iΨ₁) ? 0 : dataᵢⱼ[iΨ₁]
        ψ₂ = isnothing(iΨ₂) ? 0 : dataᵢⱼ[iΨ₂]
        ψ₃ = isnothing(iΨ₃) ? 0 : dataᵢⱼ[iΨ₃]
        ψ₄ = isnothing(iΨ₄) ? 0 : dataᵢⱼ[iΨ₄]
        σ  = isnothing(iσ)  ? 0 : dataᵢⱼ[iσ]
        h  = isnothing(ih)  ? 0 : dataᵢⱼ[ih]
        News = isnothing(iNews) ? 0 : dataᵢⱼ[iNews]
        φ₀ = isnothing(iφ₀) ? 0 : dataᵢⱼ[iφ₀]
        φ₁ = isnothing(iφ₁) ? 0 : dataᵢⱼ[iφ₁]
        φ₂ = isnothing(iφ₂) ? 0 : dataᵢⱼ[iφ₂]

        if Eᴵ == +1
            if !isnothing(iΨ₀)
                dataᵢⱼ[iΨ₀] = k⁻³ * (
                    ψ₀ - ðu′╱k * (4ψ₁ - ðu′╱k * (6ψ₂ - ðu′╱k * (4ψ₃ - ðu′╱k * ψ₄)))
                )
            end
            if !isnothing(iΨ₁)
                dataᵢⱼ[iΨ₁] = k⁻³ * (
                    ψ₁ - ðu′╱k * (3ψ₂ - ðu′╱k * (3ψ₃ - ðu′╱k * ψ₄))
                )
            end
            if !isnothing(iΨ₂)
                dataᵢⱼ[iΨ₂] = k⁻³ * (
                    ψ₂ - ðu′╱k * (2ψ₃ - ðu′╱k * ψ₄)
                )
            end
            if !isnothing(iΨ₃)
                dataᵢⱼ[iΨ₃] = k⁻³ * (
                    ψ₃ - ðu′╱k * ψ₄
                )
            end
            if !isnothing(iΨ₄)
                dataᵢⱼ[iΨ₄] = k⁻³ * (
                    ψ₄
                )
            end
            if !isnothing(iσ)
                dataᵢⱼ[iσ] = k⁻¹ * (σ + ð²α)
            end
            if !isnothing(ih)
                dataᵢⱼ[ih] = k⁻¹ * (h + ð̄²α)
            end
            if !isnothing(iNews)
                dataᵢⱼ[iNews] = k⁻² * News
            end
            if !isnothing(iφ₀)
                dataᵢⱼ[iφ₀] = k⁻² * (
                    φ₀ - ðu′╱k * (2φ₁ - ðu′╱k * φ₂)
                )
            end
            if !isnothing(iφ₁)
                dataᵢⱼ[iφ₁] = k⁻² * (
                    φ₁ - ðu′╱k * φ₂
                )
            end
            if !isnothing(iφ₂)
                dataᵢⱼ[iφ₂] = k⁻² * (
                    φ₂
                )
            end
        else  # Eᴵ == -1
            if !isnothing(iΨ₄)
                dataᵢⱼ[iΨ₄] = k⁻³ * (
                    ψ₄ - ðu′╱k * (4ψ₃ - ðu′╱k * (6ψ₂ - ðu′╱k * (4ψ₁ - ðu′╱k * ψ₀)))
                )
            end
            if !isnothing(iΨ₃)
                dataᵢⱼ[iΨ₃] = k⁻³ * (
                    ψ₃ - ðu′╱k * (3ψ₂ - ðu′╱k * (3ψ₁ - ðu′╱k * ψ₀))
                )
            end
            if !isnothing(iΨ₂)
                dataᵢⱼ[iΨ₂] = k⁻³ * (
                    ψ₂ - ðu′╱k * (2ψ₁ - ðu′╱k * ψ₀)
                )
            end
            if !isnothing(iΨ₁)
                dataᵢⱼ[iΨ₁] = k⁻³ * (
                    ψ₁ - ðu′╱k * ψ₀
                )
            end
            if !isnothing(iΨ₀)
                dataᵢⱼ[iΨ₀] = k⁻³ * (
                    ψ₀
                )
            end
            if !isnothing(iσ)
                dataᵢⱼ[iσ] = k⁻¹ * (σ - ð²α)
            end
            if !isnothing(ih)
                dataᵢⱼ[ih] = k⁻¹ * (h - ð̄²α)
            end
            if !isnothing(iNews)
                dataᵢⱼ[iNews] = k⁻² * News
            end
            if !isnothing(iφ₂)
                dataᵢⱼ[iφ₂] = k⁻² * (
                    φ₂ - ðu′╱k * (2φ₁ - ðu′╱k * φ₀)
                )
            end
            if !isnothing(iφ₁)
                dataᵢⱼ[iφ₁] = k⁻² * (
                    φ₁ - ðu′╱k * φ₀
                )
            end
            if !isnothing(iφ₀)
                dataᵢⱼ[iφ₀] = k⁻² * (
                    φ₀
                )
            end
        end
    end
end
