# The analytic gradient of the objective.  The expensive part of each pixel's transformation
# is differentiated backwards by hand; the cheap per-pixel geometry is differentiated with
# ForwardDiff's dual numbers, with a handful of partials.

@doc raw"""
    objective_and_gradient!(∇, θ, data, t, dc, plan; ℓₘₐₓ₀=nothing)

Return the objective of [`transform_objective`](@ref) at the BMS transformation described
by the parameter vector `θ`, and store its gradient with respect to `θ` in `∇`.  The cost is
about twice that of the objective alone, whatever the number of parameters, where
ForwardDiff's costs about one evaluation per parameter.

The parameters are
```math
θ = (u_x, u_y, u_z,\; b_x, b_y, b_z,\; a_1, a_2, …, a_{(ℓ_α+1)^2}),
```
read by [`bms_from_parameters`](@ref): the rotation vector ``𝐮`` gives the frame rotation
``R = \exp(𝐮/2)``, the rapidity vector ``𝐛`` gives the boost velocity
``v⃗ = \tanh(|𝐛|)\, 𝐛/|𝐛|``, and the ``a`` are the real degrees of freedom of the
supertranslation's modes up to ``ℓ_α``.  Every ``θ`` is a valid transformation, the
identity is ``θ = 0``, and nothing is singular there, so an optimizer may start at zero and
move freely.  Use [`parameters_from_bms`](@ref) to start elsewhere.

The arguments `data`, `t`, `dc`, and `plan` are as for `transform_objective`, including
plans with `ℓ_fit`, which cost one more pass over the pixels.  The shape of `∇` must match
`θ`.
"""
function objective_and_gradient!(
    ∇::AbstractVector,
    θ::AbstractVector{<:Real},
    data::Array{Complex{T1},3},
    t::Vector{T2},
    dc::DataComponents,
    plan::ObjectivePlan{TQ};
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {T1<:Real,T2<:Real,TQ<:AbstractFloat}
    P = length(θ)
    length(∇) == P || throw(ArgumentError("∇ has length $(length(∇)), but θ has length $P"))
    T = promote_type(float(eltype(θ)), T1, T2, TQ)
    θᵈ = [seeded(ParameterTag, T(θ[p]), p, Val(P)) for p ∈ 1:P]
    v⃗, R, αᵢₙ = bms_from_parameters(θᵈ)
    J, g = gradient_core(data, t, v⃗, R, αᵢₙ, dc, plan, ℓₘₐₓ₀, T, Val(P))
    ∇ .= g
    return J
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
# dual numbers in v⃗, R, and αᵢₙ hold.  The geometry of each pixel (its rest-frame rotor, the conformal
# factor, the supertranslation and its eth derivatives there, and the null-rotation
# parameter) is computed with those dual numbers, which is cheap: a few operations per
# pixel and parameter.  The expensive part, synthesizing each pixel's time series from the
# modes, interpolating, and mixing, runs on plain numbers, and its gradient is propagated
# backwards by hand: through the mixing by its matrix, through the spline by its transpose,
# and through the synthesis by one more contraction with the modes, as costly as the
# synthesis itself.  The dependence on the eight real scalars of a pixel is found by
# running the interpolation and mixing once more with dual numbers in those scalars.  The
# derivation is on the documentation page "The gradient of the objective".
function gradient_core(
    data::Array{Complex{T1},3},
    t::Vector{T2},
    v⃗,
    R,
    αᵢₙ,
    dc::DataComponents{C,I},
    plan::ObjectivePlan{TQ},
    ℓₘₐₓ₀,
    ::Type{T},
    ::Val{P},
) where {T1,T2,C,I,TQ,T,P}
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    spins == plan.spins || throw(
        ArgumentError(
            "`dc` has components of spin weights $spins, but the plan was built " *
            "for spin weights $(plan.spins)",
        ),
    )
    R′ₚ, w, t′, v, c, target = plan.R′ₚ, plan.w, plan.t′, plan.v, plan.c, plan.target
    Nᵖ = length(R′ₚ)
    Nᵗ′ = length(t′)

    # The geometry, with derivatives with respect to the parameters, exactly as the
    # objective computes it.
    c_α = dc.conventions.c_α
    γ = lorentz_factor(v⃗)
    vˣ, vʸ, vᶻ = vec(v⃗)
    Λ = Boost(I * v⃗) * R
    cubic_spline_cache = CubicSplineCache(t)
    Rₚ, αₚ, ðαₚ, ð²αₚ = supertranslation_on_grid(R′ₚ, Λ, αᵢₙ, c_α)
    Tₚ = eltype(Rₚ)
    validate_t′(t′, t, αₚ, Rₚ, v⃗, I)
    ðt′╱2κₚ = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)
    κ⁻¹ₚ = [γ * (1 - I * v⃗dotk̂(Rₚ[i], vˣ, vʸ, vᶻ)) for i ∈ 1:Nᵖ]

    unique_spins = unique(spins)
    spin_components = [findall(==(s), spins) for s ∈ unique_spins]
    spin_range = minimum(spins):maximum(spins)
    Dσ = ForwardDiff.Dual{PixelTag,T,8}

    # Fill `xᵢ` with pixel i's time series, synthesized from the modes with plain numbers,
    # leaving the harmonics there, with their derivatives, in `Yᵈ`.
    function synthesize!(xᵢ, Y, Yᵈ, calc, i)
        sYlm_rows!(Yᵈ, calc, Rₚ[i])
        Y .= valueof.(Yᵈ)
        for (s, ks) ∈ zip(unique_spins, spin_components)
            Yˢ = view(Y, :, s - first(spin_range) + 1)
            for k ∈ ks
                @inbounds for j ∈ 1:Nᵗ
                    xᵢ[k, j] = pixel_dot(Yˢ, view(data, :, j, k), s ^ 2 + 1, N₀)
                end
            end
        end
        return xᵢ
    end
    scalars(i) = (
        valueof(κ⁻¹ₚ[i]),
        valueof(αₚ[i]),
        valueof(ðt′╱2κₚ[1, i]),
        valueof(ðt′╱2κₚ[2, i]),
        valueof(ð²αₚ[i]),
    )
    chunks = OhMyThreads.index_chunks(1:Nᵖ; n=nthreads())

    # With ℓ_fit, the objective compares the projections of the transformed waveform and of
    # the target onto the modes up to ℓ_fit, which no single pixel determines, so a first
    # pass computes the transformed waveform's projection.  Its difference from the
    # target's, R̂, is what each pixel's adjoint is made from in the second pass.
    R̂ = if plan.ℓ_fit ≥ 0
        Pproj = OhMyThreads.tmapreduce(+, chunks; chunking=false) do chunk
            calc = sYlmCalculator(one(Tₚ), ℓₘₐₓ, spin_range)
            Yᵈ = Matrix{Complex{basetype(Tₚ)}}(undef, Nᵐ, length(spin_range))
            Y = Matrix{Complex{T}}(undef, Nᵐ, length(spin_range))
            xᵢ = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ)
            ẍᵢ = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ)
            d′ᵢ = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′)
            Pc = [zeros(Complex{T}, size(plan.Ȳ[k], 1), Nᵗ′) for k ∈ 1:Nᵈ]
            for i ∈ chunk
                synthesize!(xᵢ, Y, Yᵈ, calc, i)
                transform_pixel!(d′ᵢ, xᵢ, ẍᵢ, cubic_spline_cache, t, t′, scalars(i)..., dc)
                for k ∈ 1:Nᵈ, j ∈ 1:Nᵗ′
                    @inbounds @simd for n ∈ axes(plan.Ȳ[k], 1)
                        Pc[k][n, j] = muladd(plan.Ȳ[k][n, i], d′ᵢ[k, j], Pc[k][n, j])
                    end
                end
            end
            return Pc
        end
        [Pproj[k] .- plan.T̂[k] for k ∈ 1:Nᵈ]
    else
        nothing
    end

    J, ∇ = OhMyThreads.tmapreduce(
        ((J₁, g₁), (J₂, g₂)) -> (J₁ + J₂, g₁ .+ g₂), chunks; chunking=false
    ) do chunk
        calc = sYlmCalculator(one(Tₚ), ℓₘₐₓ, spin_range)
        Yᵈ = Matrix{Complex{basetype(Tₚ)}}(undef, Nᵐ, length(spin_range))
        Y = Matrix{Complex{T}}(undef, Nᵐ, length(spin_range))
        G = Matrix{Complex{T}}(undef, Nᵐ, length(spin_range))
        xᵢ = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ)
        ẍᵢ = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ)
        d′ᵢ = Matrix{Complex{Dσ}}(undef, Nᵈ, Nᵗ′)
        # A name such as ∇x holds 2∂J/∂x*: the gradient with respect to the real and
        # imaginary parts of x, packed into one complex number.
        ∇d′ = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′)  # adjoint of the transformed series
        ∇y = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ′)
        ∇x = Matrix{Complex{T}}(undef, Nᵈ, Nᵗ)  # adjoint of the synthesized series
        M = Matrix{Complex{T}}(undef, Nᵈ, Nᵈ)
        e, z = Vector{Complex{T}}(undef, Nᵈ), Vector{Complex{T}}(undef, Nᵈ)
        ∇ẍ, ẑ = Vector{Complex{T}}(undef, Nᵗ), Vector{Complex{T}}(undef, Nᵗ)
        τ = Vector{T}(undef, Nᵗ′)
        ∇σ = Vector{T}(undef, 8)
        Jc = zero(T)
        gc = zeros(T, P)
        for i ∈ chunk
            synthesize!(xᵢ, Y, Yᵈ, calc, i)

            # Interpolate and mix with dual numbers in the pixel's eight real scalars.
            κ⁻¹, α, ð₀, ð₁, ð² = scalars(i)
            σ = (
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
            transform_pixel!(d′ᵢ, xᵢ, ẍᵢ, cubic_spline_cache, t, t′, σ..., dc)

            # The adjoint of the transformed series, so that a change δd′ changes J by
            # Re Σ conj(∇d′) δd′.  For the full objective it is the pixel's own weighted
            # residual; with ℓ_fit it is the projection's residual pulled back to the
            # pixel through the weighted harmonics.
            @inbounds for j ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                if isnothing(R̂)
                    r = valueof(d′ᵢ[k, j]) - target[i, j, k]
                    Jc += w[i] * v[j] * c[k] * abs2(r)
                    ∇d′[k, j] = 2 * w[i] * v[j] * c[k] * r
                else
                    acc = zero(Complex{T})
                    for n ∈ axes(plan.Ȳ[k], 1)
                        acc = muladd(R̂[k][n, j], conj(plan.Ȳ[k][n, i]), acc)
                    end
                    ∇d′[k, j] = 2 * v[j] * c[k] * acc
                end
            end

            # The parameters reach the pixel's scalars through the geometry.
            fill!(∇σ, 0)
            @inbounds for j ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                ∂re, ∂im = reim_partials(d′ᵢ[k, j])
                for n ∈ 1:8
                    ∇σ[n] += real(∇d′[k, j]) * ∂re[n] + imag(∇d′[k, j]) * ∂im[n]
                end
            end
            for (n, q) ∈ enumerate((
                κ⁻¹ₚ[i],
                αₚ[i],
                real(ðt′╱2κₚ[1, i]),
                imag(ðt′╱2κₚ[1, i]),
                real(ðt′╱2κₚ[2, i]),
                imag(ðt′╱2κₚ[2, i]),
                real(ð²αₚ[i]),
                imag(ð²αₚ[i]),
            ))
                ∂q = ForwardDiff.partials(q)
                for p ∈ 1:P
                    gc[p] += ∇σ[n] * ∂q[p]
                end
            end

            # They also reach it through the harmonics at the pixel's rotor, which enter
            # the synthesis.  Propagate the adjoint back to the synthesized time series:
            # through the mixing of each output time by its matrix, then through the
            # spline by its transpose.
            for j ∈ 1:Nᵗ′
                τ[j] = t′[j] * κ⁻¹ + α
                mixing_matrix!(M, κ⁻¹, ð₀ + τ[j] * ð₁, ð², dc; e, z)
                mul!(view(∇y, :, j), M', view(∇d′, :, j))
            end
            fill!(∇x, 0)
            for k ∈ 1:Nᵈ
                @views spline_transpose!(
                    ∇x[k, :], ∇y[k, :], cubic_spline_cache, t, τ; ∇ẍ, z=ẑ
                )
            end
            # Through the synthesis, xᵢ[k, j] = Σₘ Yₘ data[m, j, k], so for each spin weight
            # the derivative of J with respect to Yₘ is Gₘ = Σ data[m, j, k] conj(∇x[k, j]).
            fill!(G, 0)
            for (s, ks) ∈ zip(unique_spins, spin_components)
                col = s - first(spin_range) + 1
                for k ∈ ks, j ∈ 1:Nᵗ
                    conj∇x = conj(∇x[k, j])
                    @inbounds @simd for m ∈ (s ^ 2 + 1):N₀
                        G[m, col] = muladd(data[m, j, k], conj∇x, G[m, col])
                    end
                end
                for m ∈ (s ^ 2 + 1):N₀
                    ∂re, ∂im = reim_partials(Yᵈ[m, col])
                    Gᵐ = G[m, col]
                    for p ∈ 1:P
                        gc[p] += ∂re[p] * real(Gᵐ) - ∂im[p] * imag(Gᵐ)
                    end
                end
            end
        end
        return (Jc, gc)
    end
    if !isnothing(R̂)
        J = sum(c[k] * sum(v[j] * sum(abs2, view(R̂[k], :, j)) for j ∈ 1:Nᵗ′) for k ∈ 1:Nᵈ)
    end
    if plan.normalize
        J /= plan.norm²
        ∇ ./= plan.norm²
    end
    return J, ∇
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
freedom of `α`, which must number a perfect square `(ℓ_α+1)²`.  These are stored in the
order of the modes: at the position of mode `(ℓ, m)` is the real part of `α_{ℓ,m}` if
`m ≥ 0`, and the imaginary part of `α_{ℓ,|m|}` if `m < 0`; the reality condition then fixes
the modes of negative `m`.  Generic over the element type, so `θ` may hold dual numbers.
"""
function bms_from_parameters(θ::AbstractVector)
    Nα = length(θ) - 6
    L = isqrt(max(Nα, 0))
    (Nα ≥ 1 && L^2 == Nα) || throw(
        ArgumentError(
            "θ has $(length(θ)) entries; it needs 6 plus a perfect square number of " *
            "supertranslation parameters",
        ),
    )
    𝐮 = QuatVec(θ[1], θ[2], θ[3])
    f = tanhc_sqrt(θ[4]^2 + θ[5]^2 + θ[6]^2)
    v⃗ = QuatVec(f * θ[4], f * θ[5], f * θ[6])
    a = view(θ, 7:length(θ))
    α = zeros(Complex{eltype(θ)}, Nα)
    for ℓ ∈ 0:(L - 1)
        i₀ = ℓ^2 + ℓ + 1
        α[i₀] = a[i₀]
        for m ∈ 1:ℓ
            αₘ = complex(a[i₀ + m], a[i₀ - m])
            α[i₀ + m] = αₘ
            α[i₀ - m] = (-1)^m * conj(αₘ)
        end
    end
    return v⃗, exp(𝐮 / 2), α
end

"""
    parameters_from_bms(v⃗, R, α)

The parameter vector `θ` of [`objective_and_gradient!`](@ref) for the boost velocity `v⃗`,
frame rotation `R`, and supertranslation modes `α`; the inverse of
[`bms_from_parameters`](@ref), up to the sign of `R`, which describes the same rotation.
The modes are made real first, as everywhere else, by averaging each with its partner.
"""
function parameters_from_bms(v⃗, R, α::AbstractVector)
    v⃗ = QuatVec(v⃗)
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
    𝐮 = 2 * log(rotor_argument(R))
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
    𝐛 = g .* vec(v⃗)
    return [vec(𝐮)..., 𝐛..., a...]
end

"""
    spline_transpose!(∇x, ∇y, cubic_spline_cache, t, τ; ∇ẍ, z)

Add to `∇x` (length `Nᵗ`) the transpose of natural-cubic-spline interpolation applied to
`∇y`: if `S` is the linear map from the values `x` at the knots `t` to the values `y` of
the spline at the times `τ` (length `Nᵗ′`), as `transform_pixel!` evaluates them, this
adds `Sᵀ ∇y`.
Interpolation is real-linear, so for complex values the transpose is the adjoint.  The
interval of each time is chosen exactly as in `transform_pixel!`, and the weights of the
knot values and second derivatives come from `spline_eval` itself, so the transpose matches
the evaluation to roundoff.  The second derivatives solve a symmetric tridiagonal system,
so their transpose is another solve with the same Thomas factors.  The keywords `∇ẍ` and
`z` are work vectors of length `Nᵗ`, for callers that want to avoid allocating them.
"""
function spline_transpose!(
    ∇x,
    ∇y,
    cubic_spline_cache::CubicSplineCache,
    t,
    τ;
    ∇ẍ=zeros(eltype(∇x), length(t)),
    z=similar(∇ẍ),
)
    Nᵗ = length(t)
    h, h⁻¹, l, u⁻¹ = cubic_spline_cache.h,
    cubic_spline_cache.h⁻¹, cubic_spline_cache.l,
    cubic_spline_cache.u⁻¹
    fill!(∇ẍ, 0)  # adjoint of the second derivatives
    for j′ ∈ eachindex(τ, ∇y)
        j = clamp(searchsortedlast(t, τ[j′]), 1, Nᵗ - 1)
        Δ = τ[j′] - t[j]
        𝟙, 𝟘 = one(eltype(h)), zero(eltype(h))
        ∇x[j] += spline_eval(𝟙, 𝟘, 𝟘, 𝟘, h[j], h⁻¹[j], Δ) * ∇y[j′]
        ∇x[j + 1] += spline_eval(𝟘, 𝟙, 𝟘, 𝟘, h[j], h⁻¹[j], Δ) * ∇y[j′]
        ∇ẍ[j] += spline_eval(𝟘, 𝟘, 𝟙, 𝟘, h[j], h⁻¹[j], Δ) * ∇y[j′]
        ∇ẍ[j + 1] += spline_eval(𝟘, 𝟘, 𝟘, 𝟙, h[j], h⁻¹[j], Δ) * ∇y[j′]
    end
    # The end second derivatives are fixed at zero, so only the interior ones have
    # adjoints; solve A z = ∇ẍ on the interior with the Thomas factors, as the forward
    # sweep of `transform_pixel!` does, since A is symmetric.
    z[2] = ∇ẍ[2] * u⁻¹[1]
    for j ∈ 3:(Nᵗ - 1)
        z[j] = (∇ẍ[j] - h[j - 1] * z[j - 1]) * u⁻¹[j - 1]
    end
    for j ∈ (Nᵗ - 2):-1:2
        z[j] -= l[j - 1] * z[j + 1]
    end
    # Transpose of r[j] = 6(h⁻¹[j](x[j+1] − x[j]) − h⁻¹[j-1](x[j] − x[j-1])), j = 2..Nᵗ-1.
    for j ∈ 2:(Nᵗ - 1)
        ∇x[j + 1] += 6 * h⁻¹[j] * z[j]
        ∇x[j] -= 6 * (h⁻¹[j] + h⁻¹[j - 1]) * z[j]
        ∇x[j - 1] += 6 * h⁻¹[j - 1] * z[j]
    end
    return ∇x
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
