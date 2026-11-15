# The analytic gradient of the objective.  The expensive part of each pixel's transformation
# is differentiated backwards by hand; the cheap per-pixel geometry is differentiated with
# ForwardDiff's dual numbers, with a handful of partials.

@doc raw"""
    objective_and_gradient!(∇, θ, data, t, dc, plan; ℓₘₐₓ₀=nothing)

Return the objective of [`transform_objective`](@ref) at the BMS transformation described
by the parameter vector `θ`, and store its gradient with respect to `θ` in `∇`.  The cost is
a few times that of the objective alone, whatever the number of parameters, where
ForwardDiff's costs about one evaluation per parameter.

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
    θᵈ = [seeded(ParameterTag, T(θ[p]), p, Val(P)) for p ∈ 1:P]
    v⃗, R, αᵢₙ = bms_from_parameters(θᵈ)
    Φ, Ωθ = gradient_core(data, t, v⃗, R, αᵢₙ, dc, plan, ℓₘₐₓ₀, T, Val(P))
    ∇ .= Ωθ
    return Φ
end

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

# The objective and its gradient with respect to the P parameters whose derivatives the
# dual numbers in v⃗, R, and αᵢₙ hold.  The geometry of each pixel (its rest-frame rotor,
# the conformal factor, the supertranslation and its eth derivatives there, and the
# null-rotation parameter) is computed with those dual numbers, which is cheap: a few
# operations per pixel and parameter.  The expensive part, synthesizing each pixel's time
# series from the modes, interpolating, and mixing, runs on plain numbers, and its gradient
# is propagated backwards by hand: through the mixing by its matrix, through the spline by
# its transpose, and through the synthesis by one more product with the modes, as costly as
# the synthesis itself.  The synthesis and its transpose are matrix products over a block
# of pixels.  The dependence on the eight real scalars of a pixel is found by running the
# interpolation and mixing once more with dual numbers in those scalars.  The derivation is
# on the documentation page "The gradient of the objective".
function gradient_core(
    data::Array{Complex{T1},3},
    t::AbstractVector{T2},
    v⃗,
    R,
    αᵢₙ,
    dc::DataComponents{C,I},
    plan::ObjectivePlan{TQ},
    ℓₘₐₓ₀,
    ::Type{T},
    ::Val{P},
) where {T1,T2,C,I,TQ,T,P}
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(data, t, v⃗, R, αᵢₙ, dc, ℓₘₐₓ₀)
    validate_plan(plan, dc, N₀)
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    R′ₚ, w, t′, v, u, c = plan.R′ₚ, plan.w, plan.t′, plan.v, plan.c, plan.target
    Nᵖ = length(R′ₚ)
    Nᵗ′ = length(t′)

    # The geometry, with derivatives with respect to the parameters, exactly as the
    # objective computes it.
    c_α = dc.conventions.c_α
    γ = lorentz_factor(v⃗)
    vˣ, vʸ, vᶻ = vec(v⃗)
    # At ℐ⁺, where `I = 1`, this is the frame rotor `Λ = Boost(v⃗) * R`, which gives the new
    # legs 𝐞′ = Λ 𝐞 Λ̃ (see "Frames and Components" in the docs).  At ℐ⁻, where `I = -1`, it is
    # the parity conjugate `ℂconj(Λ) = Boost(-v⃗) * R`, which realizes the frame rotor on the
    # antipodal labeling of the sphere.
    Λℐ = Boost(I * v⃗) * R
    cubic_spline_cache = CubicSplineCache(t)
    Rₚ, αₚ, ðαₚ, ð²αₚ = supertranslation_on_grid(R′ₚ, Λℐ, αᵢₙ, c_α)
    Tₚ = eltype(Rₚ)
    validate_t′(t′, t, αₚ, Rₚ, v⃗, I)
    ðt′╱2κₚ = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)
    κ⁻¹ₚ = [γ * (1 - I * v⃗dotk̂(Rₚ[j], vˣ, vʸ, vᶻ)) for j ∈ 1:Nᵖ]

    unique_spins = unique(spins)
    spin_range = minimum(spins):maximum(spins)
    Dg = ForwardDiff.Dual{PixelTag,T,8}

    # The geometry without its derivatives.  `Rotor{T}` keeps the rotors' components
    # exactly, without renormalizing them, so harmonics evaluated at these rotors are the
    # values of the harmonics evaluated with derivatives.
    Rₚ⁰ = [Rotor{T}(valueof.(components(r))...) for r ∈ Rₚ]
    κ⁻¹ₚ⁰ = valueof.(κ⁻¹ₚ)
    αₚ⁰ = valueof.(αₚ)
    ðt′╱2κₚ⁰ = valueof.(ðt′╱2κₚ)
    ð²αₚ⁰ = valueof.(ð²αₚ)
    scalars(j) = (κ⁻¹ₚ⁰[j], αₚ⁰[j], ðt′╱2κₚ⁰[1, j], ðt′╱2κₚ⁰[2, j], ð²αₚ⁰[j])

    # With ℓ_fit, the objective compares the projections of the transformed waveform and of
    # the target onto the modes up to ℓ_fit, which no single pixel determines, so a first
    # pass computes the transformed waveform's projection, exactly as the objective does.
    # Its difference from the target's, R̂, is what each pixel's adjoint is made from in
    # the main pass.
    R̂ = if plan.ℓ_fit ≥ 0
        geometry⁰ = (; Rₚ=Rₚ⁰, κ⁻¹ₚ=κ⁻¹ₚ⁰, αₚ=αₚ⁰, ðt′╱2κₚ=ðt′╱2κₚ⁰, ð²αₚ=ð²αₚ⁰)
        Pproj = transformed_projection(
            data, t, t′, cubic_spline_cache, geometry⁰, dc, spins, N₀, plan.Ȳ, T
        )
        [Pproj[k] .- plan.T̂[k] for k ∈ 1:Nᵈ]
    else
        nothing
    end

    # The main pass takes the pixels a block at a time.  The time series of a block are
    # synthesized from the harmonics without derivatives.  With ℓ_fit, the cotangents of
    # the block's transformed series come from R̂ by one more matrix product per
    # component, into `ΩD′ᵇ`.
    blocks = synthesis_blocks(
        data,
        Rₚ⁰,
        spins;
        extra_bytes_per_pixel=isnothing(R̂) ? 0 : Nᵗ′ * Nᵈ * sizeof(Complex{T}),
    )
    B = maximum(length, blocks)
    ΩD′ᵇ = Array{Complex{T}}(undef, isnothing(R̂) ? (0, 0, 0) : (B, Nᵗ′, Nᵈ))
    conjR̂ = isnothing(R̂) ? nothing : [conj.(R̂[k]) for k ∈ 1:Nᵈ]
    Yᵇ = Array{Complex{T}}(undef, Nᵐ, B, length(spin_range))  # harmonics, then cotangents
    Xᵇ = Array{Complex{T}}(undef, Nᵗ, B, Nᵈ)  # synthesized series, then cotangents
    calcs⁰ = [sYlmCalculator(one(eltype(Rₚ⁰)), ℓₘₐₓ, spin_range) for _ ∈ 1:nthreads()]
    # One set of buffers for each parallel task.  Cotangents are in the canonical
    # representation of the gradient page: a name such as Ωx holds ∂/∂(Re x) − i ∂/∂(Im x)
    # of the objective, so that a change δx changes it by Re Σ Ωx δx.  Each is a covector,
    # and each backward step applies the transpose of its forward step's Jacobian, with no
    # conjugation.
    buffers = [
        (;
            calc=sYlmCalculator(one(Tₚ), ℓₘₐₓ, spin_range),
            Yᵈ=Matrix{Complex{basetype(Tₚ)}}(undef, Nᵐ, length(spin_range)),
            xⱼ=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ),
            ẍⱼ=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ),
            d′ⱼ=Matrix{Complex{Dg}}(undef, Nᵈ, Nᵗ′),
            Ωd′=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′),  # cotangent of the transformed series
            Ωy=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′),  # cotangent of the interpolated series
            Ωx=Matrix{Complex{T}}(undef, Nᵈ, Nᵗ),  # cotangent of the synthesized series
            M=Matrix{Complex{T}}(undef, Nᵈ, Nᵈ),
            e=Vector{Complex{T}}(undef, Nᵈ),
            z=Vector{Complex{T}}(undef, Nᵈ),
            Ωẍ=Vector{Complex{T}}(undef, Nᵗ),
            ẑ=Vector{Complex{T}}(undef, Nᵗ),
            τ=Vector{T}(undef, Nᵗ′),
            Ωg=Vector{T}(undef, 8),  # cotangent of the eight real numbers
        ) for _ ∈ 1:nthreads()
    ]
    # Each pixel's share of the objective and of its gradient (the cotangent of θ, ∂Φ/∂θᵖ,
    # a covector), summed in pixel order at the end.
    Φₚ = zeros(T, Nᵖ)
    Ωθₚ = zeros(T, P, Nᵖ)

    for block ∈ blocks
        nb = length(block)
        Yb, Xb = view(Yᵇ, :, 1:nb, :), view(Xᵇ, :, 1:nb, :)
        harmonics_block!(Yb, calcs⁰, Rₚ⁰, block)
        synthesize_block!(Xb, data, Yb, spins, N₀)
        if !isnothing(R̂)
            # ΩD′ᵇ[q, i, k] = Σₙ Ȳ[k][n, block[q]] conj(R̂[k][n, i])
            for k ∈ 1:Nᵈ
                mul!(view(ΩD′ᵇ, 1:nb, :, k), transpose(view(plan.Ȳ[k], :, block)), conjR̂[k])
            end
        end

        foreach_chunk(nb, buffers) do b, chunk
            (; xⱼ, ẍⱼ, d′ⱼ, Ωd′, Ωy, Ωx, M, e, z, Ωẍ, ẑ, τ, Ωg) = b
            for q ∈ chunk
                j = block[q]
                permutedims!(xⱼ, view(Xb, :, q, :), (2, 1))

                # Interpolate and mix with dual numbers in the pixel's eight real scalars.
                κ⁻¹, α, ð₀, ð₁, ð² = scalars(j)
                g = (
                    seeded(PixelTag, T(κ⁻¹), 1, Val(8)),
                    seeded(PixelTag, T(α), 2, Val(8)),
                    complex(
                        seeded(PixelTag, T(real(ð₀)), 3, Val(8)),
                        seeded(PixelTag, T(imag(ð₀)), 4, Val(8)),
                    ),
                    complex(
                        seeded(PixelTag, T(real(ð₁)), 5, Val(8)),
                        seeded(PixelTag, T(imag(ð₁)), 6, Val(8)),
                    ),
                    complex(
                        seeded(PixelTag, T(real(ð²)), 7, Val(8)),
                        seeded(PixelTag, T(imag(ð²)), 8, Val(8)),
                    ),
                )
                transform_pixel!(d′ⱼ, xⱼ, ẍⱼ, cubic_spline_cache, t, t′, g..., dc)

                # The cotangent of the transformed series.  For the full objective it is
                # the pixel's own weighted residual, conjugated; with ℓ_fit it is the
                # conjugated residual of the projections, pulled back to the pixel through
                # the weighted harmonics, as computed for the whole block above.
                @inbounds for i ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                    if isnothing(R̂)
                        r = valueof(d′ⱼ[k, i]) - c[j, i, k]
                        Φₚ[j] += w[j] * v[i] * u[k] * abs2(r)
                        Ωd′[k, i] = 2 * w[j] * v[i] * u[k] * conj(r)
                    else
                        Ωd′[k, i] = 2 * v[i] * u[k] * ΩD′ᵇ[q, i, k]
                    end
                end

                # The parameters reach the pixel's scalars through the geometry.
                # Each is Re Σ Ωd′ ∂d′/∂g, contracting the cotangent with a column of the
                # forward-mode Jacobian of the interpolation and mixing.
                fill!(Ωg, 0)
                @inbounds for i ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                    ∂re, ∂im = reim_partials(d′ⱼ[k, i])
                    for n ∈ 1:8
                        Ωg[n] += real(Ωd′[k, i]) * ∂re[n] - imag(Ωd′[k, i]) * ∂im[n]
                    end
                end
                for (n, gₙ) ∈ enumerate((
                    κ⁻¹ₚ[j],
                    αₚ[j],
                    real(ðt′╱2κₚ[1, j]),
                    imag(ðt′╱2κₚ[1, j]),
                    real(ðt′╱2κₚ[2, j]),
                    imag(ðt′╱2κₚ[2, j]),
                    real(ð²αₚ[j]),
                    imag(ð²αₚ[j]),
                ))
                    ∂gₙ = ForwardDiff.partials(gₙ)
                    for p ∈ 1:P
                        Ωθₚ[p, j] += Ωg[n] * ∂gₙ[p]
                    end
                end

                # They also reach it through the harmonics at the pixel's rotor, which enter
                # the synthesis.  Pull the cotangent back to the synthesized time series:
                # through the mixing at each output time by the transpose of its matrix,
                # then through the spline by its transpose.  The result replaces the
                # pixel's time series in the block.
                for i ∈ 1:Nᵗ′
                    τ[i] = t′[i] * κ⁻¹ + α
                    mixing_matrix!(M, κ⁻¹, ð₀ + τ[i] * ð₁, ð², dc; e, z)
                    mul!(view(Ωy, :, i), transpose(M), view(Ωd′, :, i))
                end
                fill!(Ωx, 0)
                for k ∈ 1:Nᵈ
                    @views spline_transpose!(
                        Ωx[k, :], Ωy[k, :], cubic_spline_cache, t, τ; Ωẍ, z=ẑ
                    )
                end
                permutedims!(view(Xb, :, q, :), Ωx, (2, 1))
            end
        end

        # Through the synthesis, xⱼ[k, h] = Σₘ Yₘ data[m, h, k], so for each spin weight
        # the cotangent of the harmonic Yₘ is ΩYₘ = Σ data[m, h, k] Ωx[k, h].  These
        # replace the harmonics in the block.
        synthesis_transpose_block!(Yb, data, Xb, spins, N₀)

        # The parameters' share is Re Σ ΩYₘ ∂Yₘ/∂θ, for which the harmonics are evaluated
        # again, this time with their derivatives.  Storing those for the whole block would
        # take P + 1 times the memory of the harmonics themselves.
        foreach_chunk(nb, buffers) do b, chunk
            (; calc, Yᵈ) = b
            for q ∈ chunk
                j = block[q]
                sYlm_rows!(Yᵈ, calc, Rₚ[j])
                for s ∈ unique_spins
                    col = s - first(spin_range) + 1
                    for m ∈ (s ^ 2 + 1):N₀
                        ∂re, ∂im = reim_partials(Yᵈ[m, col])
                        ΩYᵐ = Yb[m, q, col]
                        for p ∈ 1:P
                            Ωθₚ[p, j] += ∂re[p] * real(ΩYᵐ) - ∂im[p] * imag(ΩYᵐ)
                        end
                    end
                end
            end
        end
    end
    Φ = sum(Φₚ)
    Ωθ = vec(sum(Ωθₚ; dims=2))
    if !isnothing(R̂)
        Φ = sum(u[k] * sum(v[i] * sum(abs2, view(R̂[k], :, i)) for i ∈ 1:Nᵗ′) for k ∈ 1:Nᵈ)
    end
    if plan.normalize
        Φ /= plan.norm²
        Ωθ ./= plan.norm²
    end
    return Φ, Ωθ
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
    spline_transpose!(Ωx, Ωy, cubic_spline_cache, t, τ; Ωẍ, z)

Add to `Ωx` (length `Nᵗ`) the transpose of natural-cubic-spline interpolation applied to
`Ωy`: if `S` is the linear map from the values `x` at the knots `t` to the values `y` of
the spline at the times `τ` (length `Nᵗ′`), as `transform_pixel!` evaluates them, this
adds `Sᵀ Ωy`, which pulls a cotangent of `y` back to one of `x`.  `S` is real, so it acts
alike on real and imaginary parts.  The
interval of each time is chosen exactly as in `transform_pixel!`, and the weights of the
knot values and second derivatives come from `spline_eval` itself, so the transpose matches
the evaluation to roundoff.  The second derivatives solve a symmetric tridiagonal system,
so their transpose is another solve with the same Thomas factors.  The keywords `Ωẍ` and
`z` are work vectors of length `Nᵗ`, for callers that want to avoid allocating them.
"""
function spline_transpose!(
    Ωx,
    Ωy,
    cubic_spline_cache::CubicSplineCache,
    t,
    τ;
    Ωẍ=zeros(eltype(Ωx), length(t)),
    z=similar(Ωẍ),
)
    Nᵗ = length(t)
    h, h⁻¹, l, u⁻¹ = cubic_spline_cache.h,
    cubic_spline_cache.h⁻¹, cubic_spline_cache.l,
    cubic_spline_cache.u⁻¹
    fill!(Ωẍ, 0)  # adjoint of the second derivatives
    for j′ ∈ eachindex(τ, Ωy)
        j = clamp(searchsortedlast(t, τ[j′]), 1, Nᵗ - 1)
        Δ = τ[j′] - t[j]
        𝟙, 𝟘 = one(eltype(h)), zero(eltype(h))
        Ωx[j] += spline_eval(𝟙, 𝟘, 𝟘, 𝟘, h[j], h⁻¹[j], Δ) * Ωy[j′]
        Ωx[j + 1] += spline_eval(𝟘, 𝟙, 𝟘, 𝟘, h[j], h⁻¹[j], Δ) * Ωy[j′]
        Ωẍ[j] += spline_eval(𝟘, 𝟘, 𝟙, 𝟘, h[j], h⁻¹[j], Δ) * Ωy[j′]
        Ωẍ[j + 1] += spline_eval(𝟘, 𝟘, 𝟘, 𝟙, h[j], h⁻¹[j], Δ) * Ωy[j′]
    end
    # The end second derivatives are fixed at zero, so only the interior ones have
    # adjoints; solve A z = Ωẍ on the interior with the Thomas factors, as the forward
    # sweep of `transform_pixel!` does, since A is symmetric.
    z[2] = Ωẍ[2] * u⁻¹[1]
    for j ∈ 3:(Nᵗ - 1)
        z[j] = (Ωẍ[j] - h[j - 1] * z[j - 1]) * u⁻¹[j - 1]
    end
    for j ∈ (Nᵗ - 2):-1:2
        z[j] -= l[j - 1] * z[j + 1]
    end
    # Transpose of r[j] = 6(h⁻¹[j](x[j+1] − x[j]) − h⁻¹[j-1](x[j] − x[j-1])), j = 2..Nᵗ-1.
    for j ∈ 2:(Nᵗ - 1)
        Ωx[j + 1] += 6 * h⁻¹[j] * z[j]
        Ωx[j] -= 6 * (h⁻¹[j] + h⁻¹[j - 1]) * z[j]
        Ωx[j - 1] += 6 * h⁻¹[j - 1] * z[j]
    end
    return Ωx
end

@doc raw"""
    synthesis_transpose_block!(ΩY, data, ΩX, spins, N₀)

Apply the transpose of [`synthesize_block!`](@ref) to the cotangents `ΩX` of a block's
synthesized time series, which pulls them back to cotangents of the harmonics:
```math
ΩY[m, q, s - s₀ + 1] = ∑_{k:\, s_k = s}\; ∑_h \mathrm{data}[m, h, k]\, ΩX[h, q, k],
```
for each spin weight ``s`` among the components' spin weights `spins`, whose least is
``s₀``, and each ``m`` from ``s^2 + 1`` to `N₀`.  There is no conjugation, as for every
cotangent here, and the other entries of `ΩY` are left alone.  As in `synthesize_block!`,
data and cotangents of the same BLAS element type use one matrix product for each
component, and other types use parallel tasks over the pixels.
"""
function synthesis_transpose_block!(
    ΩY::StridedArray{T,3}, data::StridedArray{T,3}, ΩX::StridedArray{T,3}, spins, N₀
) where {T<:LinearAlgebra.BlasComplex}
    s₀ = minimum(spins)
    for s ∈ unique(spins)
        m = (s ^ 2 + 1):N₀
        ΩYˢ = view(ΩY, m, :, s - s₀ + 1)
        for (n, k) ∈ enumerate(findall(==(s), spins))
            mul!(ΩYˢ, view(data, m, :, k), view(ΩX, :, :, k), true, n > 1)
        end
    end
    return ΩY
end
function synthesis_transpose_block!(ΩY, data, ΩX, spins, N₀)
    s₀ = minimum(spins)
    OhMyThreads.tforeach(axes(ΩX, 2)) do q
        for s ∈ unique(spins)
            col = s - s₀ + 1
            fill!(view(ΩY, (s ^ 2 + 1):N₀, q, col), 0)
            for k ∈ findall(==(s), spins), h ∈ axes(ΩX, 1)
                ΩXₕₖ = ΩX[h, q, k]
                @inbounds @simd for m ∈ (s ^ 2 + 1):N₀
                    ΩY[m, q, col] = muladd(data[m, h, k], ΩXₕₖ, ΩY[m, q, col])
                end
            end
        end
    end
    return ΩY
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
