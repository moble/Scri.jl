# The analytic gradient of the objective.  The expensive part of each pixel's transformation
# is differentiated backwards by hand; the cheap per-pixel geometry is differentiated with
# ForwardDiff's dual numbers, with a handful of partials.

@doc raw"""
    objective_and_gradient!(∇, θ, data, t, dc, plan; ℓₘₐₓ₀=nothing)

Return the objective of [`transform_objective`](@ref) at the BMS transformation described by
the parameter vector `θ`, and store its gradient with respect to `θ` in `∇`.  The cost is a
few times that of the objective alone, and it grows only slowly with the number of
parameters, whereas ForwardDiff's gradient costs at least one evaluation of the objective
for each parameter.

The parameters are
```math
θ = (r_x, r_y, r_z,\; η_x, η_y, η_z,\; a_1, a_2, …, a_{(ℓ_α+1)^2}),
```
read by [`bms_from_parameters`](@ref): the rotation vector ``\vec{r}`` gives the frame
rotation ``R = \exp(\vec{r}/2)``, the rapidity vector ``\vec{η}`` gives the boost velocity
``v⃗ = \tanh(|\vec{η}|)\, \vec{η}/|\vec{η}|``, and the ``a`` are the real degrees of
freedom of the supertranslation's modes up to ``ℓ_α``.  Every ``θ`` describes a valid
BMS transformation, the identity is ``θ = 0``, and nothing is singular there, so an
optimizer may start at zero.  Use [`parameters_from_bms`](@ref) to start elsewhere.  The
plan's output times `t′` must still be valid for every ``θ`` the optimizer tries, as for
`transform_objective`: a grid from `compute_t′(t, β_max, δt_max)` covers boosts up to speed
`β_max` and supertranslations up to `δt_max`, and beyond them this function throws an
`ArgumentError`.

The arguments `data`, `t`, `dc`, and `plan` are as for `transform_objective`, including
plans with `ℓ_fit`, which cost one more pass over the pixels.  The shape of `∇` must match
`θ`.
"""
function objective_and_gradient!(
    ∇::AbstractVector,
    θ::AbstractVector{<:Real},
    data::Array{Complex{T1},3},
    t::AbstractVector{T2},
    dc::DataComponents,
    plan::ObjectivePlan{TQ};
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {T1<:Real,T2<:Real,TQ<:AbstractFloat}
    P = length(θ)
    length(∇) == P || throw(ArgumentError("∇ has length $(length(∇)), but θ has length $P"))
    T = promote_type(float(eltype(θ)), T1, T2, TQ)
    Φ, ∂Φ╱∂θ = gradient_core(
        data, t, convert(Vector{T}, θ), dc, plan, ℓₘₐₓ₀, Val(parameter_chunk_width)
    )
    ∇ .= ∂Φ╱∂θ
    return Φ
end

# The number of parameters whose derivatives the dual numbers of the geometry hold at one
# time.  The geometry is computed once for each chunk of this many parameters, so the
# gradient is compiled once for any number of parameters, and its dual numbers stay narrow.
const parameter_chunk_width = 8

# Tags for the two kinds of dual number below: one holds derivatives with respect to the
# parameters θ, the other with respect to the eight real scalars of a single pixel.
struct ParameterTag end
struct PixelTag end

valueof(z::Complex) = complex(ForwardDiff.value(real(z)), ForwardDiff.value(imag(z)))
valueof(x::Real) = ForwardDiff.value(x)
reim_partials(z::Complex) = (ForwardDiff.partials(real(z)), ForwardDiff.partials(imag(z)))

function seeded(::Type{Tag}, x::T, p, ::Val{P}) where {Tag,T,P}
    return ForwardDiff.Dual{Tag}(
        x, ForwardDiff.Partials(ntuple(q -> q == p ? one(T) : zero(T), Val(P)))
    )
end

# The objective and its gradient with respect to the parameters θ.  The geometry of each
# pixel (its rotor in the original frame, the conformal factor, the supertranslation and
# its eth derivatives there, and the null-rotation parameter) is computed with dual numbers
# that hold the derivatives with respect to W of the parameters, once for each chunk of W
# parameters.  That is cheap: a few operations per pixel and parameter.  The expensive part,
# synthesizing each pixel's time series from the modes, interpolating, and mixing, runs on
# plain numbers, and its gradient is propagated backwards by hand: through the mixing by its
# matrix, through the spline by its transpose, and through the synthesis by one more product
# with the modes, as costly as the synthesis itself.  The synthesis and its transpose are
# matrix products over a block of pixels.  The dependence on the eight real scalars of a
# pixel is found by running the interpolation and mixing once more with dual numbers in
# those scalars.  The derivation is on the documentation page "The gradient of the
# objective".
function gradient_core(
    data::Array{Complex{T1},3},
    t::AbstractVector{T2},
    θ::Vector{T},
    dc::DataComponents{C,I},
    plan::ObjectivePlan{TQ},
    ℓₘₐₓ₀,
    ::Val{W},
) where {T1,T2,T,C,I,TQ,W}
    P = length(θ)
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    R′ₚ, w, t′, v, u, c = plan.R′ₚ, plan.w, plan.t′, plan.v, plan.u, plan.c
    Nᵖ = length(R′ₚ)
    Nᵗ′ = length(t′)

    # The transformation, with its derivatives with respect to one chunk of the parameters
    # at a time.  In the dual numbers of the chunk `ps`, slot `a` holds the derivative with
    # respect to parameter `ps[a]`.  The last chunk may have slots to spare, whose
    # derivatives are zero.
    chunks = [p:min(p + W - 1, P) for p ∈ 1:W:max(P, 1)]
    transformations = map(chunks) do ps
        return bms_from_parameters([
            seeded(ParameterTag, θ[p], p - first(ps) + 1, Val(W)) for p ∈ 1:P
        ])
    end
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(
        data, t, first(transformations)..., dc, ℓₘₐₓ₀
    )
    validate_plan(plan, dc, N₀)

    # The geometry of every pixel for each chunk, exactly as the objective computes it.
    cubic_spline_cache = CubicSplineCache(t)
    geometries = [
        pixel_geometry(t, t′, v⃗, R, αᵢₙ, R′ₚ, dc) for (v⃗, R, αᵢₙ) ∈ transformations
    ]
    Tₚ = eltype(first(geometries).Rₚ)

    unique_spins = unique(spins)
    spin_range = minimum(spins):maximum(spins)
    Dg = ForwardDiff.Dual{PixelTag,T,8}

    # The geometry without its derivatives, which is the same in every chunk.  `Rotor{T}`
    # keeps the rotors' components exactly, without renormalizing them, so harmonics
    # evaluated at these rotors are the values of the harmonics evaluated with derivatives.
    G⁰ = first(geometries)
    Rₚ⁰ = [Rotor{T}(valueof.(components(r))...) for r ∈ G⁰.Rₚ]
    κ⁻¹ₚ⁰ = valueof.(G⁰.κ⁻¹ₚ)
    αₚ⁰ = valueof.(G⁰.αₚ)
    ðt′╱2κₚ⁰ = valueof.(G⁰.ðt′╱2κₚ)
    ð²αₚ⁰ = valueof.(G⁰.ð²αₚ)
    scalars(j) = (κ⁻¹ₚ⁰[j], αₚ⁰[j], ðt′╱2κₚ⁰[1, j], ðt′╱2κₚ⁰[2, j], ð²αₚ⁰[j])

    # With ℓ_fit, the objective compares the projections of the transformed waveform and of
    # the target onto the modes up to ℓ_fit, which no single pixel determines, so a first
    # pass computes the transformed waveform's projection, μ̂, exactly as the objective does.
    # Its difference from the target's, Π, is what each pixel's ∂ν╱∂μ is made from in the
    # main pass.
    Π = if plan.ℓ_fit ≥ 0
        geometry⁰ = (; Rₚ=Rₚ⁰, κ⁻¹ₚ=κ⁻¹ₚ⁰, αₚ=αₚ⁰, ðt′╱2κₚ=ðt′╱2κₚ⁰, ð²αₚ=ð²αₚ⁰)
        μ̂ = transformed_projection(
            data, t, t′, cubic_spline_cache, geometry⁰, dc, spins, N₀, plan.Ȳ, T
        )
        [μ̂[k] .- plan.ĉ[k] for k ∈ 1:Nᵈ]
    else
        nothing
    end

    # The main pass takes the pixels a block at a time.  The time series of a block are
    # synthesized from the harmonics without derivatives.  With ℓ_fit, the conjugate of Π
    # is evaluated back at the block's pixels by one more matrix product per component,
    # into `wr̄ᵇ`: each pixel's weight times the conjugate of the projected residual there.
    blocks = synthesis_blocks(
        data,
        Rₚ⁰,
        spins;
        extra_bytes_per_pixel=isnothing(Π) ? 0 : Nᵗ′ * Nᵈ * sizeof(Complex{T}),
    )
    B = maximum(length, blocks)
    wr̄ᵇ = Array{Complex{T}}(undef, isnothing(Π) ? (0, 0, 0) : (B, Nᵗ′, Nᵈ))
    Π̄ = isnothing(Π) ? nothing : [conj.(Π[k]) for k ∈ 1:Nᵈ]
    Yᵇ = Array{Complex{T}}(undef, Nᵐ, B, length(spin_range))  # the harmonics, then ∂ν╱∂Y
    σᵇ = Array{Complex{T}}(undef, Nᵗ, B, Nᵈ)  # the synthesized series σ, then ∂ν╱∂σ
    ∂ν╱∂gᵇ = Matrix{T}(undef, 8, B)  # ∂ν╱∂g for each pixel of the block
    calcs⁰ = [sYlmCalculator(one(eltype(Rₚ⁰)), ℓₘₐₓ, spin_range) for _ ∈ 1:nthreads()]
    # One set of buffers for each parallel task.  The names are those of the gradient page,
    # for one pixel: `σⱼ`, `μⱼ`, and the derivatives of the pixel's share ν of the objective.
    # A name such as `∂ν╱∂σ` holds the derivative of ν with respect to the complex σ, with σ
    # and its conjugate treated as independent: it is (∂/∂(Re σ) − i ∂/∂(Im σ))/2, so that a
    # change δσ changes ν by 2 Re Σ ∂ν╱∂σ δσ.  Each is a covector, and each backward step
    # applies the transpose of its forward step's Jacobian, with no conjugation.
    buffers = [
        (;
            calc=sYlmCalculator(one(Tₚ), ℓₘₐₓ, spin_range),
            Yᵈ=Matrix{Complex{basetype(Tₚ)}}(undef, Nᵐ, length(spin_range)),
            σⱼ=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ),
            σ̈ⱼ=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ),  # the spline's second derivatives
            μⱼ=Matrix{Complex{Dg}}(undef, Nᵈ, Nᵗ′),  # μ, with its derivatives ∂μ/∂g
            ∂ν╱∂μ=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′),
            ∂ν╱∂ι=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′),
            ∂ν╱∂σ=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ),
            M=Matrix{Complex{T}}(undef, Nᵈ, Nᵈ),
            e=Vector{Complex{T}}(undef, Nᵈ),
            z=Vector{Complex{T}}(undef, Nᵈ),
            ∂ν╱∂σ̈=Vector{Complex{T}}(undef, Nᵗ),
            ζ=Vector{Complex{T}}(undef, Nᵗ),
            τ=Vector{T}(undef, Nᵗ′),
        ) for _ ∈ 1:nthreads()
    ]
    # Each pixel's share of the objective, ν[j], and the gradient of that share with respect
    # to the parameters, ∂ν╱∂θ[p, j], a covector.  Both are summed in pixel order at the end.
    ν = zeros(T, Nᵖ)
    ∂ν╱∂θ = zeros(T, P, Nᵖ)

    for block ∈ blocks
        nb = length(block)
        Yb, σb = view(Yᵇ, :, 1:nb, :), view(σᵇ, :, 1:nb, :)
        # The backward pass overwrites the block's σ with ∂ν╱∂σ, and then its harmonics with
        # ∂ν╱∂Y, so these are two more names for the same memory.
        ∂ν╱∂Yb, ∂ν╱∂σb = Yb, σb
        harmonics_block!(Yb, calcs⁰, Rₚ⁰, block)
        synthesize_block!(σb, data, Yb, spins, N₀)
        if !isnothing(Π)
            # wr̄ᵇ[q, i, k] = Σₙ Ȳ[k][n, block[q]] conj(Π[k][n, i])
            for k ∈ 1:Nᵈ
                mul!(view(wr̄ᵇ, 1:nb, :, k), transpose(view(plan.Ȳ[k], :, block)), Π̄[k])
            end
        end

        foreach_chunk(nb, buffers) do b, chunk
            (; σⱼ, σ̈ⱼ, μⱼ, ∂ν╱∂μ, ∂ν╱∂ι, ∂ν╱∂σ, M, e, z, ∂ν╱∂σ̈, ζ, τ) = b
            for q ∈ chunk
                j = block[q]
                permutedims!(σⱼ, view(σb, :, q, :), (2, 1))

                # Interpolate and mix with dual numbers in the pixel's eight real numbers g.
                κ⁻¹, α, b⁰, b¹, ð²α = scalars(j)
                g = (
                    seeded(PixelTag, T(κ⁻¹), 1, Val(8)),
                    seeded(PixelTag, T(α), 2, Val(8)),
                    complex(
                        seeded(PixelTag, T(real(b⁰)), 3, Val(8)),
                        seeded(PixelTag, T(imag(b⁰)), 4, Val(8)),
                    ),
                    complex(
                        seeded(PixelTag, T(real(b¹)), 5, Val(8)),
                        seeded(PixelTag, T(imag(b¹)), 6, Val(8)),
                    ),
                    complex(
                        seeded(PixelTag, T(real(ð²α)), 7, Val(8)),
                        seeded(PixelTag, T(imag(ð²α)), 8, Val(8)),
                    ),
                )
                transform_pixel!(μⱼ, σⱼ, σ̈ⱼ, cubic_spline_cache, t, t′, g..., dc)

                # Through ν.  For the full objective, ∂ν╱∂μ is the pixel's own weighted
                # residual, conjugated.  With ℓ_fit, the residual is replaced by its
                # projection onto the modes up to ℓ_fit, evaluated back at the pixel, as
                # computed for the whole block above.
                @inbounds for i ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                    if isnothing(Π)
                        r = valueof(μⱼ[k, i]) - c[j, i, k]
                        ν[j] += w[j] * v[i] * u[k] * abs2(r)
                        ∂ν╱∂μ[k, i] = w[j] * v[i] * u[k] * conj(r)
                    else
                        ∂ν╱∂μ[k, i] = v[i] * u[k] * wr̄ᵇ[q, i, k]
                    end
                end

                # The parameters reach the pixel through its eight numbers g, by
                # ∂ν╱∂g = 2 Re Σ ∂ν╱∂μ ∂μ/∂g, where the ∂μ/∂g are the partials of `μⱼ`.
                ∂ν╱∂g = view(∂ν╱∂gᵇ, :, q)
                fill!(∂ν╱∂g, 0)
                @inbounds for i ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                    ∂Reμ╱∂g, ∂Imμ╱∂g = reim_partials(μⱼ[k, i])
                    for n ∈ 1:8
                        ∂ν╱∂g[n] +=
                            2 * (
                                real(∂ν╱∂μ[k, i]) * ∂Reμ╱∂g[n] -
                                imag(∂ν╱∂μ[k, i]) * ∂Imμ╱∂g[n]
                            )
                    end
                end

                # They also reach it through the harmonics at the pixel's rotor, which enter
                # the synthesis.  Pull ∂ν╱∂μ back to the synthesized series: through μ, at
                # each output time, by the transpose of its matrix, and then through ι by
                # the spline's transpose.  The result replaces the pixel's series in the
                # block.
                for i ∈ 1:Nᵗ′
                    τ[i] = t′[i] * κ⁻¹ + α
                    mixing_matrix!(M, κ⁻¹, b⁰ + τ[i] * b¹, ð²α, dc; e, z)
                    mul!(view(∂ν╱∂ι, :, i), transpose(M), view(∂ν╱∂μ, :, i))
                end
                fill!(∂ν╱∂σ, 0)
                for k ∈ 1:Nᵈ
                    @views spline_transpose!(
                        ∂ν╱∂σ[k, :], ∂ν╱∂ι[k, :], cubic_spline_cache, t, τ; ∂ν╱∂σ̈, ζ
                    )
                end
                permutedims!(view(∂ν╱∂σb, :, q, :), ∂ν╱∂σ, (2, 1))
            end
        end

        # Through σ.  The synthesis is σⱼ[k, h] = Σₘ Yₘ data[m, h, k], so for each spin
        # weight ∂ν╱∂Yₘ = Σ ∂ν╱∂σ[k, h] data[m, h, k].  These replace the harmonics in the
        # block.
        synthesis_transpose_block!(∂ν╱∂Yb, data, ∂ν╱∂σb, spins, N₀)

        # Assemble the gradient of each pixel's ν,
        # ∂ν╱∂θ = ∂ν╱∂g ∂g/∂θ + 2 Re Σ ∂ν╱∂Yₘ ∂Yₘ/∂θ, one chunk of the parameters at a time.
        # The derivatives ∂g/∂θ are the partials of the chunk's geometry.  For ∂Yₘ/∂θ the
        # harmonics are evaluated again, at the chunk's dual rotor.  Storing those for the
        # whole block would take P + 1 times the memory of the harmonics themselves.
        foreach_chunk(nb, buffers) do b, chunk
            (; calc, Yᵈ) = b
            for q ∈ chunk
                j = block[q]
                for (ps, G) ∈ zip(chunks, geometries)
                    g = (
                        G.κ⁻¹ₚ[j],
                        G.αₚ[j],
                        reim(G.ðt′╱2κₚ[1, j])...,
                        reim(G.ðt′╱2κₚ[2, j])...,
                        reim(G.ð²αₚ[j])...,
                    )
                    # The chunk's W slots of ∂ν╱∂θ for this pixel.
                    ∂ν╱∂θᶜ = zero(ForwardDiff.Partials{W,T})
                    for n ∈ 1:8
                        ∂ν╱∂θᶜ += ∂ν╱∂gᵇ[n, q] * ForwardDiff.partials(g[n])
                    end
                    sYlm_rows!(Yᵈ, calc, G.Rₚ[j])
                    for s ∈ unique_spins
                        col = s - first(spin_range) + 1
                        for m ∈ (s ^ 2 + 1):N₀
                            ∂ReY╱∂θ, ∂ImY╱∂θ = reim_partials(Yᵈ[m, col])
                            ∂ν╱∂Yₘ = ∂ν╱∂Yb[m, q, col]
                            ∂ν╱∂θᶜ += 2 * (∂ReY╱∂θ * real(∂ν╱∂Yₘ) - ∂ImY╱∂θ * imag(∂ν╱∂Yₘ))
                        end
                    end
                    for (a, p) ∈ enumerate(ps)
                        ∂ν╱∂θ[p, j] = ∂ν╱∂θᶜ[a]
                    end
                end
            end
        end
    end
    Φ = sum(ν)
    ∂Φ╱∂θ = vec(sum(∂ν╱∂θ; dims=2))
    if !isnothing(Π)
        Φ = sum(u[k] * sum(v[i] * sum(abs2, view(Π[k], :, i)) for i ∈ 1:Nᵗ′) for k ∈ 1:Nᵈ)
    end
    if plan.normalize
        Φ /= plan.norm²
        ∂Φ╱∂θ ./= plan.norm²
    end
    return Φ, ∂Φ╱∂θ
end

"""
    tanhc_sqrt(x)

Return `tanh(√x)/√x`, the factor that turns a rapidity vector into a velocity, for `x ≥ 0`.
Near zero it is a Taylor series in `x`, so that it is smooth there and automatic
differentiation never meets the infinite derivative of `√x` at zero.  The series is used
below `eps^(1/5)` in the precision of `x`, where the first term it omits, about `x⁵/113`,
is below roundoff, so the switch is exact to roundoff in any floating-point type.
"""
function tanhc_sqrt(x::T) where {T}
    PT = primal_float(T)
    if x < eps(PT)^(one(PT) / 5)
        return 1 + x * (-1 // 3 + x * (2 // 15 + x * (-17 // 315 + x * (62 // 2835))))
    end
    s = sqrt(x)
    return tanh(s) / s
end

"""
    bms_from_parameters(θ)

Return the boost velocity `v⃗`, the frame rotation `R`, and the supertranslation modes `α`
that the parameter vector `θ` of [`objective_and_gradient!`](@ref) describes, in the forms
that [`transform!`](@ref) and [`transform_objective`](@ref) take.  `θ[1:3]` is the
rotation vector, with `R = exp(QuatVec(θ[1:3]...)/2)`; `θ[4:6]` is the rapidity vector,
with `v⃗` of speed `tanh(|θ[4:6]|)` along it; and `θ[7:end]` holds the real degrees of
freedom of `α`, which must number a perfect square `(ℓ_α+1)²`, or zero for no
supertranslation.  These are stored in the
order of the modes: at the position of mode `(ℓ, m)` is the real part of `α_{ℓ,m}` if
`m ≥ 0`, and the imaginary part of `α_{ℓ,|m|}` if `m < 0`; the reality condition then fixes
the modes of negative `m`.  Generic over the element type, so `θ` may hold dual numbers.
"""
function bms_from_parameters(θ::AbstractVector)
    Nα = length(θ) - 6
    L = isqrt(max(Nα, 0))
    (Nα ≥ 0 && L^2 == Nα) || throw(
        ArgumentError(
            "θ has $(length(θ)) entries; it needs 6 plus a perfect square number of " *
            "supertranslation parameters",
        ),
    )
    r⃗ = QuatVec(θ[1], θ[2], θ[3])
    f = tanhc_sqrt(θ[4]^2 + θ[5]^2 + θ[6]^2)
    v⃗ = QuatVec(f * θ[4], f * θ[5], f * θ[6])
    a = view(θ, 7:length(θ))
    α = zeros(Complex{eltype(θ)}, max(Nα, 1))  # six parameters: no supertranslation
    for ℓ ∈ 0:(L - 1)
        i₀ = ℓ^2 + ℓ + 1
        α[i₀] = a[i₀]
        for m ∈ 1:ℓ
            αₘ = complex(a[i₀ + m], a[i₀ - m])
            α[i₀ + m] = αₘ
            α[i₀ - m] = (-1)^m * conj(αₘ)
        end
    end
    return v⃗, exp(r⃗ / 2), α
end

"""
    parameters_from_bms(v⃗, R, α)

The parameter vector `θ` of [`objective_and_gradient!`](@ref) for the boost velocity `v⃗`,
frame rotation `R`, and supertranslation modes `α`; the inverse of
[`bms_from_parameters`](@ref), up to the sign of `R`, which describes the same rotation;
the rotation vector returned has length at most π.
The modes are made real first, as everywhere else, by averaging each with its partner.
"""
function parameters_from_bms(v⃗, R, α::AbstractVector)
    v⃗ = as_quatvec(v⃗)
    β² = abs2vec(v⃗)
    β² < 1 || throw(ArgumentError("The boost speed must be less than 1"))
    # atanh(β)/β, the inverse of `tanhc_sqrt`, by its series near zero, which is used below
    # `eps^(1/4)`, where the first term it omits, β⁸/9, is below roundoff.
    PT = primal_float(typeof(β²))
    g = if β² < eps(PT)^(one(PT) / 4)
        1 + β² * (1 // 3 + β² * (1 // 5 + β² * (1 // 7)))
    else
        atanh(sqrt(β²)) / sqrt(β²)
    end
    # `±R` are the same rotation; the sign with a nonnegative scalar part gives the
    # shortest rotation vector, so that `-1` gives zero rather than a vector of length 2π.
    R = rotor_argument(R)
    r⃗ = 2 * log(R.w < 0 ? -R : R)
    L = isqrt(length(α))
    L^2 == length(α) || throw(ArgumentError("α must have a perfect-square length"))
    a = zeros(real(eltype(float.(α))), length(α))
    αʳ = impose_reality(collect(complex.(float.(α))), L - 1, 1)
    for ℓ ∈ 0:(L - 1)
        i₀ = ℓ^2 + ℓ + 1
        a[i₀] = real(αʳ[i₀])
        for m ∈ 1:ℓ
            a[i₀ + m] = real(αʳ[i₀ + m])
            a[i₀ - m] = imag(αʳ[i₀ + m])
        end
    end
    η⃗ = g .* vec(v⃗)
    return [vec(r⃗)..., η⃗..., a...]
end

"""
    spline_transpose!(∂ν╱∂σ, ∂ν╱∂ι, cubic_spline_cache, t, τ; ∂ν╱∂σ̈, ζ)

Take the backward step through the interpolation, for one component of one pixel.  If `S`
is the linear map from the samples `σ` at the times `t` to the values `ι` of their natural
cubic spline at the times `τ` (length `Nᵗ′`), as `transform_pixel!` evaluates them, this
adds `Sᵀ ∂ν╱∂ι` to `∂ν╱∂σ` (length `Nᵗ`), which turns derivatives with respect to `ι` into
derivatives with respect to `σ`.  `S` is real, so it acts alike on real and imaginary
parts.  The interval of each time is chosen exactly as in `transform_pixel!`, and the
weights of the knot values and second derivatives come from `spline_eval` itself, so the
transpose matches the evaluation to roundoff.  The second derivatives solve a symmetric
tridiagonal system, so their transpose is another solve with the same Thomas factors.

The keywords `∂ν╱∂σ̈` and `ζ` are work vectors of length `Nᵗ`, for callers that want to
avoid allocating them.  They hold the derivatives with respect to the spline's second
derivatives, and the solution of the tridiagonal system, which the gradient page calls
``ζ``.
"""
function spline_transpose!(
    ∂ν╱∂σ,
    ∂ν╱∂ι,
    cubic_spline_cache::CubicSplineCache,
    t,
    τ;
    ∂ν╱∂σ̈=zeros(eltype(∂ν╱∂σ), length(t)),
    ζ=similar(∂ν╱∂σ̈),
)
    Nᵗ = length(t)
    h, h⁻¹, l, u⁻¹ = cubic_spline_cache.h,
    cubic_spline_cache.h⁻¹, cubic_spline_cache.l,
    cubic_spline_cache.u⁻¹
    fill!(∂ν╱∂σ̈, 0)
    for j′ ∈ eachindex(τ, ∂ν╱∂ι)
        j = clamp(searchsortedlast(t, τ[j′]), 1, Nᵗ - 1)
        Δ = τ[j′] - t[j]
        𝟙, 𝟘 = one(eltype(h)), zero(eltype(h))
        ∂ν╱∂σ[j] += spline_eval(𝟙, 𝟘, 𝟘, 𝟘, h[j], h⁻¹[j], Δ) * ∂ν╱∂ι[j′]
        ∂ν╱∂σ[j + 1] += spline_eval(𝟘, 𝟙, 𝟘, 𝟘, h[j], h⁻¹[j], Δ) * ∂ν╱∂ι[j′]
        ∂ν╱∂σ̈[j] += spline_eval(𝟘, 𝟘, 𝟙, 𝟘, h[j], h⁻¹[j], Δ) * ∂ν╱∂ι[j′]
        ∂ν╱∂σ̈[j + 1] += spline_eval(𝟘, 𝟘, 𝟘, 𝟙, h[j], h⁻¹[j], Δ) * ∂ν╱∂ι[j′]
    end
    # The end second derivatives are fixed at zero, so only the interior ones have
    # derivatives; solve A ζ = ∂ν╱∂σ̈ on the interior with the Thomas factors, as the forward
    # sweep of `transform_pixel!` does, since A is symmetric.
    ζ[2] = ∂ν╱∂σ̈[2] * u⁻¹[1]
    for j ∈ 3:(Nᵗ - 1)
        ζ[j] = (∂ν╱∂σ̈[j] - h[j - 1] * ζ[j - 1]) * u⁻¹[j - 1]
    end
    for j ∈ (Nᵗ - 2):-1:2
        ζ[j] -= l[j - 1] * ζ[j + 1]
    end
    # Transpose of r[j] = 6(h⁻¹[j](σ[j+1] − σ[j]) − h⁻¹[j-1](σ[j] − σ[j-1])), j = 2..Nᵗ-1.
    for j ∈ 2:(Nᵗ - 1)
        ∂ν╱∂σ[j + 1] += 6 * h⁻¹[j] * ζ[j]
        ∂ν╱∂σ[j] -= 6 * (h⁻¹[j] + h⁻¹[j - 1]) * ζ[j]
        ∂ν╱∂σ[j - 1] += 6 * h⁻¹[j - 1] * ζ[j]
    end
    return ∂ν╱∂σ
end

"""
    synthesis_transpose_block!(∂ν╱∂Y, data, ∂ν╱∂σ, spins, N₀)

Take the backward step through the synthesis, for a block of pixels.  This applies the
transpose of [`synthesize_block!`](@ref) to `∂ν╱∂σ`, the derivatives with respect to the
block's synthesized time series, which gives the derivatives with respect to the harmonics:
```
∂ν╱∂Y[m, q, s - s₀ + 1] = Σₖ Σₕ data[m, h, k] ∂ν╱∂σ[h, q, k],
```
where the sum over `k` runs over the components of spin weight `s`.  This is done for each
spin weight `s` among the components' spin weights `spins`, whose least is `s₀`, and each
`m` from `s^2 + 1` to `N₀`.  There is no conjugation, as in every backward step here, and
the other entries of `∂ν╱∂Y` are left alone.  As in `synthesize_block!`, arrays of the same
BLAS element type use one matrix product for each component, and other types use parallel
tasks over the pixels.
"""
function synthesis_transpose_block!(
    ∂ν╱∂Y::StridedArray{T,3}, data::StridedArray{T,3}, ∂ν╱∂σ::StridedArray{T,3}, spins, N₀
) where {T<:LinearAlgebra.BlasComplex}
    s₀ = minimum(spins)
    for s ∈ unique(spins)
        m = (s ^ 2 + 1):N₀
        ∂ν╱∂Yˢ = view(∂ν╱∂Y, m, :, s - s₀ + 1)
        for (n, k) ∈ enumerate(findall(==(s), spins))
            mul!(∂ν╱∂Yˢ, view(data, m, :, k), view(∂ν╱∂σ, :, :, k), true, n > 1)
        end
    end
    return ∂ν╱∂Y
end
function synthesis_transpose_block!(∂ν╱∂Y, data, ∂ν╱∂σ, spins, N₀)
    s₀ = minimum(spins)
    OhMyThreads.tforeach(axes(∂ν╱∂σ, 2)) do q
        for s ∈ unique(spins)
            col = s - s₀ + 1
            fill!(view(∂ν╱∂Y, (s ^ 2 + 1):N₀, q, col), 0)
            for k ∈ findall(==(s), spins), h ∈ axes(∂ν╱∂σ, 1)
                ∂ν╱∂σₕₖ = ∂ν╱∂σ[h, q, k]
                @inbounds @simd for m ∈ (s ^ 2 + 1):N₀
                    ∂ν╱∂Y[m, q, col] = muladd(data[m, h, k], ∂ν╱∂σₕₖ, ∂ν╱∂Y[m, q, col])
                end
            end
        end
    end
    return ∂ν╱∂Y
end

"""
    mixing_matrix!(M, κ⁻¹, ðt′╱2κ, ð²α, dc; e, z)

Fill the `Nᵈ × Nᵈ` matrix `M` with the linear part of [`mix_components!`](@ref) at these
scalars, so that `mix_components!(x, …) == M * x + mix_components!(0, …)`: the mixing is
complex-linear in the components, with an offset only for the shear and the strain.  Column
`k` is the image of the `k`-th unit vector less the image of zero.  The keywords `e` and
`z` are work vectors of length `Nᵈ`.  Returns `M`.
"""
function mixing_matrix!(
    M::AbstractMatrix{Complex{T}},
    κ⁻¹,
    ðt′╱2κ,
    ð²α,
    dc;
    e=zeros(Complex{T}, size(M, 1)),
    z=zeros(Complex{T}, size(M, 1)),
) where {T}
    Nᵈ = size(M, 1)
    fill!(z, 0)
    mix_components!(z, κ⁻¹, ðt′╱2κ, ð²α, dc)
    for k ∈ 1:Nᵈ
        fill!(e, 0)
        e[k] = 1
        mix_components!(e, κ⁻¹, ðt′╱2κ, ð²α, dc)
        @views M[:, k] .= e .- z
    end
    return M
end
