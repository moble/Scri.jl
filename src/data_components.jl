const ValidDataComponents = (:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ, :λ, :h, :News, :φ₀, :φ₁, :φ₂)

# This utility function just translates any reasonable string representation of a data
# component's name into the canonical symbol.
function parse_data_component(s::AbstractString)
    return Symbol(
        replace(
            s,
            '_' => "",
            '{' => "",
            '}' => "",
            '0' => '₀',
            '1' => '₁',
            '2' => '₂',
            '3' => '₃',
            '4' => '₄',
            "σ" => "σ",
            "Σ" => "σ",
            "sigma" => "σ",
            "Sigma" => "σ",
            "SIGMA" => "σ",
            "λ" => "λ",
            "Λ" => "λ",
            "lambda" => "λ",
            "Lambda" => "λ",
            "LAMBDA" => "λ",
            "h" => "h",
            "H" => "h",
            "strain" => "h",
            "Strain" => "h",
            "STRAIN" => "h",
            "News" => "News",
            "news" => "News",
            "NEWS" => "News",
            "N" => "News",
            "n" => "News",
            "ψ" => "ψ",
            "Ψ" => "ψ",
            "psi" => "ψ",
            "Psi" => "ψ",
            "PSI" => "ψ",
            "φ" => "φ",
            "ϕ" => "φ",
            "Φ" => "φ",
            "phi" => "φ",
            "Phi" => "φ",
            "PHI" => "φ",
            "varphi" => "φ",
            "varPhi" => "φ",
            "varPHI" => "φ",
            "VARPHI" => "φ",
        ),
    )
end

"""
    DataComponents{C, ℐ, V}

Encodes a fixed set of waveform data components at the type level.  `C` is an
`NTuple{N,Symbol}` whose elements are drawn from `ValidDataComponents`:

    $ValidDataComponents

The parameter `ℐ` is the sign of the time direction for which these components are defined:
`+1` for ``ℐ⁺`` (outgoing) and `-1` for ``ℐ⁻`` (incoming).  The default value is `+1`, since
most applications will be for ``ℐ⁺``.

The `conventions` field (of type `V<:Conventions`) records which convention the component
data are expressed in; it defaults to `Conventions()`, the package-native SXS conventions.
The conventions travel with the data descriptor — [`transform!`](@ref) applies the
transformation laws *native to that convention*, and [`represent!`](@ref) converts data
between conventions.

The inputs may alternatively be strings; any reasonable spelling will be parsed into the
canonical symbol form.  For example, `DataComponents("Psi_3", "psi4", "sigma")` will be
parsed into `DataComponents(:ψ₃, :ψ₄, :σ)`.  The constructor will throw an error if the
input components are invalid.

Examples:

    DataComponents(:ψ₄)                          # gravitational waves only
    DataComponents(:ψ₄, :ψ₃, :ψ₂)                # top three Weyl components
    DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :σ)  # full Weyl set with shear (ℐ⁺)
    DataComponents(:ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄, :λ; ℐ=-1)  # full Weyl set with shear (ℐ⁻)
    DataComponents(:φ₀, :φ₁, :φ₂)                # Faraday components
    DataComponents(:φ₀, :φ₁, :φ₂; ℐ=-1)         # Faraday components on ℐ⁻
    DataComponents(:ψ₄, :h; conventions=Conventions(:MB))  # data in the MB convention
"""
struct DataComponents{C,I,V<:Conventions}
    conventions::V
    function DataComponents(
        c₁::Symbol, cs::Symbol...; ℐ=1, conventions::Conventions=Conventions()
    )
        components = (c₁, cs...)
        validate_data_components(components, ℐ)
        return new{components,ℐ,typeof(conventions)}(conventions)
    end
    function DataComponents(
        c₁::AbstractString,  # Keep one argument explicit to avoid ambiguities with 0 args
        cs::AbstractString...;
        ℐ=1,
        conventions::Conventions=Conventions(),
    )
        return DataComponents(map(parse_data_component, (c₁, cs...))...; ℐ, conventions)
    end
end

function validate_data_components(cs, ℐ)
    @assert all(c -> c ∈ ValidDataComponents, cs) "" *
        "Invalid component in $cs; allowed: $ValidDataComponents"
    @assert length(Set(cs)) == length(cs) "Duplicate components in $cs"
    if ℐ == 1
        @assert :λ ∉ cs "λ is the ℐ⁻ radiative shear; use σ on ℐ⁺ (ℐ=+1)"
        :ψ₀ ∈ cs && @assert :ψ₁ ∈ cs "ψ₀ requires ψ₁ on ℐ⁺ (ℐ=+1)"
        :ψ₁ ∈ cs && @assert :ψ₂ ∈ cs "ψ₁ requires ψ₂ on ℐ⁺ (ℐ=+1)"
        :ψ₂ ∈ cs && @assert :ψ₃ ∈ cs "ψ₂ requires ψ₃ on ℐ⁺ (ℐ=+1)"
        :ψ₃ ∈ cs && @assert :ψ₄ ∈ cs "ψ₃ requires ψ₄ on ℐ⁺ (ℐ=+1)"
        :φ₀ ∈ cs && @assert :φ₁ ∈ cs "φ₀ requires φ₁ on ℐ⁺ (ℐ=+1)"
        :φ₁ ∈ cs && @assert :φ₂ ∈ cs "φ₁ requires φ₂ on ℐ⁺ (ℐ=+1)"
    elseif ℐ == -1
        @assert :σ ∉ cs "σ is the ℐ⁺ radiative shear; use λ on ℐ⁻ (ℐ=-1)"
        :ψ₄ ∈ cs && @assert :ψ₃ ∈ cs "ψ₄ requires ψ₃ on ℐ⁻ (ℐ=-1)"
        :ψ₃ ∈ cs && @assert :ψ₂ ∈ cs "ψ₃ requires ψ₂ on ℐ⁻ (ℐ=-1)"
        :ψ₂ ∈ cs && @assert :ψ₁ ∈ cs "ψ₂ requires ψ₁ on ℐ⁻ (ℐ=-1)"
        :ψ₁ ∈ cs && @assert :ψ₀ ∈ cs "ψ₁ requires ψ₀ on ℐ⁻ (ℐ=-1)"
        :φ₂ ∈ cs && @assert :φ₁ ∈ cs "φ₂ requires φ₁ on ℐ⁻ (ℐ=-1)"
        :φ₁ ∈ cs && @assert :φ₀ ∈ cs "φ₁ requires φ₀ on ℐ⁻ (ℐ=-1)"
    else
        throw(ArgumentError("Invalid ℐ = $ℐ; must be ±1"))
    end
end

"""
    component_index(dc::DataComponents{C}, ::Val{S})

Return the index of component `S` within `dc`, or `nothing` if absent.  The return type is
inferred as a compile-time constant when `dc` has a concrete type.
"""
@inline component_index(::DataComponents{C}, ::Val{S}) where {C,S} = findfirst(==(S), C)

"""
    has_component(dc::DataComponents{C}, ::Val{S})

Return `true` if component `S` is present in `dc`.
"""
@inline has_component(::DataComponents{C}, ::Val{S}) where {C,S} = S ∈ C

"""
    ncomponents(dc::DataComponents)

Return the number of components in `dc`.
"""
ncomponents(::DataComponents{C}) where {C} = length(C)

"""
    spin_weight(::Val{S})

Return the spin weight of field component `S`.
"""
spin_weight(::Val{:ψ₀}) = 2
spin_weight(::Val{:ψ₁}) = 1
spin_weight(::Val{:ψ₂}) = 0
spin_weight(::Val{:ψ₃}) = -1
spin_weight(::Val{:ψ₄}) = -2
spin_weight(::Val{:σ}) = 2
spin_weight(::Val{:λ}) = -2
spin_weight(::Val{:h}) = -2
spin_weight(::Val{:News}) = -2
spin_weight(::Val{:φ₀}) = 1
spin_weight(::Val{:φ₁}) = 0
spin_weight(::Val{:φ₂}) = -1
Base.@constprop :aggressive spin_weight(s)::Int = spin_weight(Val(s))

"""
    conformal_weight(::Val{S})

Return the conformal weight (spin weight + boost weight) of field component `S`, which is
the power of the conformal factor ``κ`` in the BMS transformation law.

Note that we are assuming that these fields represent the asymptotic values of the physical
fields at null infinity, so they have already been rescaled by the appropriate power of the
conformal factor to be finite and nonzero at null infinity.  The transformation law for the
physical fields (finite-radius Weyl and Faraday spinors) would have a different conformal
weight, but the asymptotic fields are the ones we are transforming.

More specifically, the *asymptotic* Weyl spinor ``ψ`` and Faraday spinor ``φ`` are related
to the finite-radius Weyl spinor ``Ψ`` by ``Ψ = ωψ`` and the finite-radius Faraday spinor
``Φ`` by ``Φ = ωφ``.  The factor ``ω`` is the conformal factor that goes to zero (but has
nonzero derivative) at null infinity, which transforms as ``ω′=κω``, so we pick up a factor
of ``κ⁻¹`` in the transformation laws for ``ψ`` and ``φ`` compared to ``Ψ`` and ``Φ``.
Since ``Ψ`` and ``Φ`` are the physical quantities, they do not change under coordinate
transformations.

Meanwhile the basis spinors each transform with a factor of ``1/√κ``.  The Weyl components
``ψₙ`` are defined by contracting the Weyl spinor with *four* basis spinors, so they pick up
a factor of ``κ⁻²`` from the basis spinors and an additional factor of ``κ⁻¹`` from the
conformal factor, for a total of ``κ⁻³``.  Similarly, the Faraday components ``φₙ`` are
defined by contracting the Faraday spinor with *two* basis spinors, so they pick up a factor
of ``κ⁻¹`` from the conformal factor and an additional factor of ``κ⁻¹`` from the basis
spinors, for a total of ``κ⁻²``.
"""
conformal_weight(::Val{:ψ₀}) = -3
conformal_weight(::Val{:ψ₁}) = -3
conformal_weight(::Val{:ψ₂}) = -3
conformal_weight(::Val{:ψ₃}) = -3
conformal_weight(::Val{:ψ₄}) = -3
conformal_weight(::Val{:σ}) = -1
conformal_weight(::Val{:λ}) = -1
conformal_weight(::Val{:h}) = -1
conformal_weight(::Val{:News}) = -2
conformal_weight(::Val{:φ₀}) = -2
conformal_weight(::Val{:φ₁}) = -2
conformal_weight(::Val{:φ₂}) = -2

"""
    mix_components!(dataᵢⱼ, κ⁻¹, ðt′╱2κ, ð²α, dc)

Apply the BMS component-mixing transformation to `dataᵢⱼ`.  `κ⁻¹` is the inverse conformal
factor for this pixel, `ðt′╱2κ` is the eth-derivative of the retarded time in the new frame
divided by ``2κ`` (the *geometric* null-rotation parameter ``b = ðu'/2κ``, computed with the
native Newman–Penrose ``ð``), and `ð²α` is the second eth-derivative of the
supertranslation.  The latter enters the radiative-shear/strain law with a factor of one
half (the ½ has the same dyad/√2 origin as the ½ in ``b``).  The radiative shear is a
different quantity on each null infinity: the ``l``-congruence shear ``σ`` (spin weight
``+2``) on ``ℐ⁺``, shifting by ``+F_σ\\, ð²α/2``; and the ``n``-congruence shear ``λ`` (NP's
``λ``, spin weight ``-2``) on ``ℐ⁻``, shifting by ``+F_λ\\, ð̄²α/2``.  The strain ``h`` (spin
weight ``-2``) is present on both, shifting by ``+F_h\\, ð̄²α/2`` on ``ℐ⁺`` and ``-F_h\\,
ð̄²α/2`` on ``ℐ⁻``.

The convention factors come from the `Conventions` carried by `dc`, and appear in the laws
exactly as in the ["Convention dependence" documentation](@ref
convention_dependence_fields): the towers run on ``c_l c_m\\, ðu′/2κ`` (on ``ℐ⁺``; the
parameter is ``ð̄v′/2κ / (c_l c_m)`` on ``ℐ⁻``, whose tower mixes downward), and the
shifts use ``F_σ`` from [`shear_factor`](@ref), ``F_λ`` from [`lambda_factor`](@ref), and
``F_h`` from [`strain_factor`](@ref).  At the default conventions every factor is a `One`
singleton and compiles away.

Note that Julia specializes on the concrete type of `dc`.  This means that the indexes into
`dataᵢⱼ` for the various components are known at compile time, and the branches for which
components are present or absent will be resolved at compile time.  The result is that —
even though the function body looks unwieldy and slow — it compiles down to minimal code
with no branches and only the necessary components, making it very fast in practice.  The
`@inline` annotation encourages this specialization and inlining, especially when just a few
components are being processed.
"""
@inline function mix_components!(
    dataᵢⱼ::AbstractVector{Complex{T}}, κ⁻¹, ðt′╱2κ, ð²α, dc::DataComponents{C,I}
) where {T,C,I}
    κ⁻² = κ⁻¹ * κ⁻¹
    κ⁻³ = κ⁻² * κ⁻¹
    ð̄²α = conj(ð²α)
    F_σ = shear_factor(dc.conventions)
    F_λ = lambda_factor(dc.conventions)
    F_h = strain_factor(dc.conventions)

    iψ₄ = component_index(dc, Val(:ψ₄))
    iψ₃ = component_index(dc, Val(:ψ₃))
    iψ₂ = component_index(dc, Val(:ψ₂))
    iψ₁ = component_index(dc, Val(:ψ₁))
    iψ₀ = component_index(dc, Val(:ψ₀))
    iσ = component_index(dc, Val(:σ))
    iλ = component_index(dc, Val(:λ))
    ih = component_index(dc, Val(:h))
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
        σ = isnothing(iσ) ? 0 : dataᵢⱼ[iσ]
        λ = isnothing(iλ) ? 0 : dataᵢⱼ[iλ]
        h = isnothing(ih) ? 0 : dataᵢⱼ[ih]
        News = isnothing(iNews) ? 0 : dataᵢⱼ[iNews]
        φ₀ = isnothing(iφ₀) ? 0 : dataᵢⱼ[iφ₀]
        φ₁ = isnothing(iφ₁) ? 0 : dataᵢⱼ[iφ₁]
        φ₂ = isnothing(iφ₂) ? 0 : dataᵢⱼ[iφ₂]

        if I == +1
            # The ℐ⁺ tower mixes upward, so the convention-X mixing parameter is the SXS
            # input rescaled by the neighboring-factor ratio Fₙ/Fₙ₊₁ = c_l c_m:
            # b^{[X]} = (c_l c_m)·(ðt′/2κ).
            ðu′╱2κ = dyad_factor(dc.conventions) * ðt′╱2κ
            if !isnothing(iψ₀)
                dataᵢⱼ[iψ₀] =
                    κ⁻³ *
                    (ψ₀ + ðu′╱2κ * (4ψ₁ + ðu′╱2κ * (6ψ₂ + ðu′╱2κ * (4ψ₃ + ðu′╱2κ * ψ₄))))
            end
            if !isnothing(iψ₁)
                dataᵢⱼ[iψ₁] = κ⁻³ * (ψ₁ + ðu′╱2κ * (3ψ₂ + ðu′╱2κ * (3ψ₃ + ðu′╱2κ * ψ₄)))
            end
            if !isnothing(iψ₂)
                dataᵢⱼ[iψ₂] = κ⁻³ * (ψ₂ + ðu′╱2κ * (2ψ₃ + ðu′╱2κ * ψ₄))
            end
            if !isnothing(iψ₃)
                dataᵢⱼ[iψ₃] = κ⁻³ * (ψ₃ + ðu′╱2κ * ψ₄)
            end
            if !isnothing(iψ₄)
                dataᵢⱼ[iψ₄] = κ⁻³ * (ψ₄)
            end
            if !isnothing(iσ)
                dataᵢⱼ[iσ] = κ⁻¹ * (σ + F_σ * ð²α / 2)
            end
            if !isnothing(ih)
                dataᵢⱼ[ih] = κ⁻¹ * (h + F_h * ð̄²α / 2)
            end
            if !isnothing(iNews)
                dataᵢⱼ[iNews] = κ⁻² * News
            end
            if !isnothing(iφ₀)
                dataᵢⱼ[iφ₀] = κ⁻² * (φ₀ + ðu′╱2κ * (2φ₁ + ðu′╱2κ * φ₂))
            end
            if !isnothing(iφ₁)
                dataᵢⱼ[iφ₁] = κ⁻² * (φ₁ + ðu′╱2κ * φ₂)
            end
            if !isnothing(iφ₂)
                dataᵢⱼ[iφ₂] = κ⁻² * (φ₂)
            end
        else  # I == -1
            # The ℐ⁻ generator is l̃, so the peeling tower is the l-fixed null rotation, whose
            # parameter is the conjugate ð̄v′╱2κ = conj(ðt′╱2κ) (spin weight -1).  Only then do
            # the two terms in each rung share a spin weight, as the tower runs from ψ₀ (s=+2)
            # down to ψ₄ (s=-2): e.g. ψ₁ (s=+1) = ψ₁ + ð̄v′╱2κ (s=-1) · ψ₀ (s=+2).  Because
            # this tower mixes downward, the dyad factor divides — Fₙ/Fₙ₋₁ = 1/(c_l c_m) —
            # and the conjugation acts only on the SXS input, never on the factor:
            # b̄^{[X]} = conj(ðt′/2κ)/(c_l c_m), which is NOT conj(b^{[X]}) unless c_l² = 1.
            ð̄v′╱2κ = conj(ðt′╱2κ) / dyad_factor(dc.conventions)
            if !isnothing(iψ₄)
                dataᵢⱼ[iψ₄] =
                    κ⁻³ *
                    (ψ₄ + ð̄v′╱2κ * (4ψ₃ + ð̄v′╱2κ * (6ψ₂ + ð̄v′╱2κ * (4ψ₁ + ð̄v′╱2κ * ψ₀))))
            end
            if !isnothing(iψ₃)
                dataᵢⱼ[iψ₃] = κ⁻³ * (ψ₃ + ð̄v′╱2κ * (3ψ₂ + ð̄v′╱2κ * (3ψ₁ + ð̄v′╱2κ * ψ₀)))
            end
            if !isnothing(iψ₂)
                dataᵢⱼ[iψ₂] = κ⁻³ * (ψ₂ + ð̄v′╱2κ * (2ψ₁ + ð̄v′╱2κ * ψ₀))
            end
            if !isnothing(iψ₁)
                dataᵢⱼ[iψ₁] = κ⁻³ * (ψ₁ + ð̄v′╱2κ * ψ₀)
            end
            if !isnothing(iψ₀)
                dataᵢⱼ[iψ₀] = κ⁻³ * (ψ₀)
            end
            if !isnothing(iλ)
                dataᵢⱼ[iλ] = κ⁻¹ * (λ + F_λ * ð̄²α / 2)
            end
            if !isnothing(ih)
                dataᵢⱼ[ih] = κ⁻¹ * (h - F_h * ð̄²α / 2)
            end
            if !isnothing(iNews)
                dataᵢⱼ[iNews] = κ⁻² * News
            end
            if !isnothing(iφ₂)
                dataᵢⱼ[iφ₂] = κ⁻² * (φ₂ + ð̄v′╱2κ * (2φ₁ + ð̄v′╱2κ * φ₀))
            end
            if !isnothing(iφ₁)
                dataᵢⱼ[iφ₁] = κ⁻² * (φ₁ + ð̄v′╱2κ * φ₀)
            end
            if !isnothing(iφ₀)
                dataᵢⱼ[iφ₀] = κ⁻² * (φ₀)
            end
        end
    end
end

"""
    represent!(data, dc::DataComponents, conventions::Conventions)

Re-express `data` — mode weights with dimensions `(Nᵐ, Nᵗ, Nᵈ)`, as for
[`transform!`](@ref), described by `dc` — in the given `conventions`, in place.  Returns
`(data, dc′)`, where `dc′` is a new `DataComponents` carrying the target conventions.

Each component slice is multiplied by the ratio of its [`conversion_factor`](@ref)s,
``F^{[\\text{to}]}/F^{[\\text{from}]}`` — a constant per component, so mode weights and pixel
values convert identically.  This changes only the *description* of the data, not the
physics: converting, transforming with the target-convention laws, and converting back is
identical to transforming with the source-convention laws.

Note that the supertranslation time-law sign `c_α` is *not* a property of the data, so it
does not enter here; but it is part of the target `Conventions`, and `transform!` will read
it from `dc′`.
"""
function represent!(
    data::AbstractArray{<:Complex,3}, dc::DataComponents{C,I}, conventions::Conventions
) where {C,I}
    @assert size(data, 3) == length(C) "Input `data` has $(size(data, 3)) components, " *
        "but `dc` has $(length(C))"
    for (k, S) ∈ enumerate(C)
        f =
            conversion_factor(conventions, Val(S)) /
            conversion_factor(dc.conventions, Val(S))
        if !(f isa One)  # One() would be a no-op; skip the memory traversal
            view(data, :, :, k) .*= f
        end
    end
    return data, DataComponents(C...; ℐ=I, conventions)
end
