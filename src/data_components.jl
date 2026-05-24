const ValidDataComponents = (:Ψ₀, :Ψ₁, :Ψ₂, :Ψ₃, :Ψ₄, :σ, :h, :News)

"""
    DataComponents{C}

Encodes a fixed set of waveform data components at the type level.  `C` is an
`NTuple{N,Symbol}` whose elements are drawn from `ValidDataComponents`:

    $ValidDataComponents

Examples:

    DataComponents(:Ψ₄)                          # gravitational waves only
    DataComponents(:Ψ₄, :Ψ₃, :Ψ₂)                # top three Weyl components
    DataComponents(:Ψ₀, :Ψ₁, :Ψ₂, :Ψ₃, :Ψ₄, :σ)  # full Weyl set with strain
"""
struct DataComponents{C}
    function DataComponents(cs::Symbol...)
        validate_data_components(cs)
        new{cs}()
    end
end


function validate_data_components(cs)
    @assert all(c -> c ∈ ValidDataComponents, cs) "" *
        "Invalid component in $cs; allowed: $ValidDataComponents"
    :Ψ₀ ∈ cs && @assert :Ψ₁ ∈ cs "Ψ₀ requires Ψ₁"
    :Ψ₁ ∈ cs && @assert :Ψ₂ ∈ cs "Ψ₁ requires Ψ₂"
    :Ψ₂ ∈ cs && @assert :Ψ₃ ∈ cs "Ψ₂ requires Ψ₃"
    :Ψ₃ ∈ cs && @assert :Ψ₄ ∈ cs "Ψ₃ requires Ψ₄"
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
Base.@constprop :aggressive spin_weight(s)::Int = spin_weight(Val(s))

"""
    conformal_weight(::Val{S})

Return the conformal weight (spin weight + boost weight) of field component `S`, which is
the power of the conformal factor ``k`` in the BMS transformation law.
"""
conformal_weight(::Val{:Ψ₀}) = -3
conformal_weight(::Val{:Ψ₁}) = -3
conformal_weight(::Val{:Ψ₂}) = -3
conformal_weight(::Val{:Ψ₃}) = -3
conformal_weight(::Val{:Ψ₄}) = -3
conformal_weight(::Val{:σ})  = -1
conformal_weight(::Val{:h})  = -1
conformal_weight(::Val{:News})  = -2


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
    dc::DataComponents{C}
) where {T, C}
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

    @inbounds begin
        ψ₀ = isnothing(iΨ₀) ? 0 : dataᵢⱼ[iΨ₀]
        ψ₁ = isnothing(iΨ₁) ? 0 : dataᵢⱼ[iΨ₁]
        ψ₂ = isnothing(iΨ₂) ? 0 : dataᵢⱼ[iΨ₂]
        ψ₃ = isnothing(iΨ₃) ? 0 : dataᵢⱼ[iΨ₃]
        ψ₄ = isnothing(iΨ₄) ? 0 : dataᵢⱼ[iΨ₄]
        σ  = isnothing(iσ)  ? 0 : dataᵢⱼ[iσ]
        h  = isnothing(ih)  ? 0 : dataᵢⱼ[ih]
        News = isnothing(iNews) ? 0 : dataᵢⱼ[iNews]

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
    end
end
