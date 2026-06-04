const ValidDataComponents = (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :h, :News, :φ₀, :φ₁, :φ₂)

# This utility function just translates any reasonable string representation of a data
# component's name into the canonical symbol.
function parse_data_component(s::AbstractString)
    Symbol(replace(s,
        '_' => "", '{' => "", '}' => "",
        '0' => '₀', '1' => '₁', '2' => '₂', '3' => '₃', '4' => '₄',
        "σ" => "σ", "Σ" => "σ", "sigma" => "σ", "Sigma" => "σ", "SIGMA" => "σ",
        "h" => "h", "H" => "h", "strain" => "h", "Strain" => "h", "STRAIN" => "h",
        "News" => "News", "news" => "News", "NEWS" => "News", "N" => "News", "n" => "News",
        "ψ" => "ψ", "Ψ" => "ψ", "psi" => "ψ", "Psi" => "ψ", "PSI" => "ψ",
        "φ" => "φ", "ϕ" => "φ", "Φ" => "φ", "phi" => "φ", "Phi" => "φ", "PHI" => "φ",
        "varphi" => "φ", "varPhi" => "φ", "varPHI" => "φ", "VARPHI" => "φ"
    ))
end

"""
    DataComponents{C, εᴵ}

Encodes a fixed set of waveform data components at the type level.  `C` is an
`NTuple{N,Symbol}` whose elements are drawn from `ValidDataComponents`:

    $ValidDataComponents

The parameter `εᴵ` is the sign of the time direction for which these components are defined:
`+1` for ``ℐ⁺`` (outgoing) and `-1` for ``ℐ⁻`` (incoming).  The default value is `+1`, since
most applications will be for ``ℐ⁺``.

The inputs may alternatively be strings; any reasonable spelling will be parsed into the
canonical symbol form.  For example, `DataComponents("Psi_3", "psi4", "sigma")` will be
parsed into `DataComponents(:ψ₃, :ψ₄, :σ)`.  The constructor will throw an error if the
input components are invalid.

Examples:

    DataComponents(:ψ₄)                          # gravitational waves only
    DataComponents(:ψ₄, :ψ₃, :ψ₂)                # top three Weyl components
    DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ)  # full Weyl set with strain
    DataComponents(:φ₀, :φ₁, :φ₂)                # Faraday components
    DataComponents(:φ₀, :φ₁, :φ₂; εᴵ=-1)         # Faraday components on ℐ⁻
"""
struct DataComponents{C, Eᴵ}
    function DataComponents(cs::Symbol...; εᴵ=1)
        validate_data_components(cs, εᴵ)
        new{cs, εᴵ}()
    end
    function DataComponents(cs::AbstractString...; εᴵ=1)
        DataComponents((parse_data_component.(cs))...; εᴵ)
    end
end

function validate_data_components(cs, εᴵ)
    @assert all(c -> c ∈ ValidDataComponents, cs) "" *
        "Invalid component in $cs; allowed: $ValidDataComponents"
    @assert length(Set(cs)) == length(cs) "Duplicate components in $cs"
    if εᴵ == 1
        :ψ₀ ∈ cs && @assert :ψ₁ ∈ cs "ψ₀ requires ψ₁ on ℐ⁺ (εᴵ=+1)"
        :ψ₁ ∈ cs && @assert :ψ₂ ∈ cs "ψ₁ requires ψ₂ on ℐ⁺ (εᴵ=+1)"
        :ψ₂ ∈ cs && @assert :ψ₃ ∈ cs "ψ₂ requires ψ₃ on ℐ⁺ (εᴵ=+1)"
        :ψ₃ ∈ cs && @assert :ψ₄ ∈ cs "ψ₃ requires ψ₄ on ℐ⁺ (εᴵ=+1)"
        :φ₀ ∈ cs && @assert :φ₁ ∈ cs "φ₀ requires φ₁ on ℐ⁺ (εᴵ=+1)"
        :φ₁ ∈ cs && @assert :φ₂ ∈ cs "φ₁ requires φ₂ on ℐ⁺ (εᴵ=+1)"
    elseif εᴵ == -1
        :ψ₄ ∈ cs && @assert :ψ₃ ∈ cs "ψ₄ requires ψ₃ on ℐ⁻ (εᴵ=-1)"
        :ψ₃ ∈ cs && @assert :ψ₂ ∈ cs "ψ₃ requires ψ₂ on ℐ⁻ (εᴵ=-1)"
        :ψ₂ ∈ cs && @assert :ψ₁ ∈ cs "ψ₂ requires ψ₁ on ℐ⁻ (εᴵ=-1)"
        :ψ₁ ∈ cs && @assert :ψ₀ ∈ cs "ψ₁ requires ψ₀ on ℐ⁻ (εᴵ=-1)"
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
spin_weight(::Val{:ψ₀}) =  2
spin_weight(::Val{:ψ₁}) =  1
spin_weight(::Val{:ψ₂}) =  0
spin_weight(::Val{:ψ₃}) = -1
spin_weight(::Val{:ψ₄}) = -2
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
``ψₙ`` are defined by contracting the Weyl spinor with *four* basis spinors, so they pick up
a factor of ``k⁻²`` from the basis spinors and an additional factor of ``k⁻¹`` from the
conformal factor, for a total of ``k⁻³``.  Similarly, the Faraday components ``φₙ`` are
defined by contracting the Faraday spinor with *two* basis spinors, so they pick up a factor
of ``k⁻¹`` from the conformal factor and an additional factor of ``k⁻¹`` from the basis
spinors, for a total of ``k⁻²``.
"""
conformal_weight(::Val{:ψ₀}) = -3
conformal_weight(::Val{:ψ₁}) = -3
conformal_weight(::Val{:ψ₂}) = -3
conformal_weight(::Val{:ψ₃}) = -3
conformal_weight(::Val{:ψ₄}) = -3
conformal_weight(::Val{:σ})  = -1
conformal_weight(::Val{:h})  = -1
conformal_weight(::Val{:News})  = -2
conformal_weight(::Val{:φ₀}) = -2
conformal_weight(::Val{:φ₁}) = -2
conformal_weight(::Val{:φ₂}) = -2


"""
    mix_components!(dataᵢⱼ, k⁻¹, ðt′╱k, ð²α, dc)

Apply the BMS component-mixing transformation to `dataᵢⱼ`.  `k⁻¹` is the inverse conformal
factor for this pixel, `ðt′╱k` is the eth-derivative of the retarded time in the new frame
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
    ðt′╱k,
    ð²α,
    dc::DataComponents{C, Eᴵ}
) where {T, C, Eᴵ}
    k⁻² = k⁻¹ * k⁻¹
    k⁻³ = k⁻² * k⁻¹
    ð̄²α = conj(ð²α)

    iψ₄ = component_index(dc, Val(:ψ₄))
    iψ₃ = component_index(dc, Val(:ψ₃))
    iψ₂ = component_index(dc, Val(:ψ₂))
    iψ₁ = component_index(dc, Val(:ψ₁))
    iψ₀ = component_index(dc, Val(:ψ₀))
    iσ  = component_index(dc, Val(:σ))
    ih  = component_index(dc, Val(:h))
    iNews = component_index(dc, Val(:News))
    iφ₂ = component_index(dc, Val(:φ₂))
    iφ₁ = component_index(dc, Val(:φ₁))
    iφ₀ = component_index(dc, Val(:φ₀))

    @inbounds begin
        ψ₀ = isnothing(iψ₀) ? 0 : dataᵢⱼ[iψ₀]
        ψ₁ = isnothing(iψ₁) ? 0 : dataᵢⱼ[iψ₁]
        ψ₂ = isnothing(iψ₂) ? 0 : dataᵢⱼ[iψ₂]
        ψ₃ = isnothing(iψ₃) ? 0 : dataᵢⱼ[iψ₃]
        ψ₄ = isnothing(iψ₄) ? 0 : dataᵢⱼ[iψ₄]
        σ  = isnothing(iσ)  ? 0 : dataᵢⱼ[iσ]
        h  = isnothing(ih)  ? 0 : dataᵢⱼ[ih]
        News = isnothing(iNews) ? 0 : dataᵢⱼ[iNews]
        φ₀ = isnothing(iφ₀) ? 0 : dataᵢⱼ[iφ₀]
        φ₁ = isnothing(iφ₁) ? 0 : dataᵢⱼ[iφ₁]
        φ₂ = isnothing(iφ₂) ? 0 : dataᵢⱼ[iφ₂]

        if Eᴵ == +1
            if !isnothing(iψ₀)
                dataᵢⱼ[iψ₀] = k⁻³ * (
                    ψ₀ - ðt′╱k * (4ψ₁ - ðt′╱k * (6ψ₂ - ðt′╱k * (4ψ₃ - ðt′╱k * ψ₄)))
                )
            end
            if !isnothing(iψ₁)
                dataᵢⱼ[iψ₁] = k⁻³ * (
                    ψ₁ - ðt′╱k * (3ψ₂ - ðt′╱k * (3ψ₃ - ðt′╱k * ψ₄))
                )
            end
            if !isnothing(iψ₂)
                dataᵢⱼ[iψ₂] = k⁻³ * (
                    ψ₂ - ðt′╱k * (2ψ₃ - ðt′╱k * ψ₄)
                )
            end
            if !isnothing(iψ₃)
                dataᵢⱼ[iψ₃] = k⁻³ * (
                    ψ₃ - ðt′╱k * ψ₄
                )
            end
            if !isnothing(iψ₄)
                dataᵢⱼ[iψ₄] = k⁻³ * (
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
                    φ₀ - ðt′╱k * (2φ₁ - ðt′╱k * φ₂)
                )
            end
            if !isnothing(iφ₁)
                dataᵢⱼ[iφ₁] = k⁻² * (
                    φ₁ - ðt′╱k * φ₂
                )
            end
            if !isnothing(iφ₂)
                dataᵢⱼ[iφ₂] = k⁻² * (
                    φ₂
                )
            end
        else  # Eᴵ == -1
            if !isnothing(iψ₄)
                dataᵢⱼ[iψ₄] = k⁻³ * (
                    ψ₄ - ðt′╱k * (4ψ₃ - ðt′╱k * (6ψ₂ - ðt′╱k * (4ψ₁ - ðt′╱k * ψ₀)))
                )
            end
            if !isnothing(iψ₃)
                dataᵢⱼ[iψ₃] = k⁻³ * (
                    ψ₃ - ðt′╱k * (3ψ₂ - ðt′╱k * (3ψ₁ - ðt′╱k * ψ₀))
                )
            end
            if !isnothing(iψ₂)
                dataᵢⱼ[iψ₂] = k⁻³ * (
                    ψ₂ - ðt′╱k * (2ψ₁ - ðt′╱k * ψ₀)
                )
            end
            if !isnothing(iψ₁)
                dataᵢⱼ[iψ₁] = k⁻³ * (
                    ψ₁ - ðt′╱k * ψ₀
                )
            end
            if !isnothing(iψ₀)
                dataᵢⱼ[iψ₀] = k⁻³ * (
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
                    φ₂ - ðt′╱k * (2φ₁ - ðt′╱k * φ₀)
                )
            end
            if !isnothing(iφ₁)
                dataᵢⱼ[iφ₁] = k⁻² * (
                    φ₁ - ðt′╱k * φ₀
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
