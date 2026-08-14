"""
    validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)

Check the shared input contract of [`transform!`](@ref) and [`transform_objective`](@ref):
the boost speed is subluminal; `t` has at least 4 samples; `data` is a 3-dimensional array
whose first dimension is a perfect square and whose second and third dimensions match those
of `t` and `dc`, respectively; `αᵢₙ` fits within the band limit; and any declared input band
limit `ℓₘₐₓ₀` is in range and honored by `data`.  This function throws an `ArgumentError` on
any violation.  If there is none, it returns `(Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀)`, where `N₀` is the
number of modes within the declared input band limit (`Nᵐ` if none was declared).
"""
function validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    # These are explicit throws (not `@assert`) so they survive even when asserts are
    # disabled.
    check(cond, msg) = cond || throw(ArgumentError(msg))
    check(absvec(v⃗) < 1, "Input `v⃗` has magnitude $(absvec(v⃗)), but expected less than 1")
    check(length(t) ≥ 4, "Input `t` has only $(length(t)) samples, but expected at least 4")
    check(ndims(data) == 3, "Input `data` has $(ndims(data)) dimensions, but expected 3")
    Nᵐ, Nᵗ, Nᵈ = size(data)
    L = isqrt(Nᵐ)
    check(L^2 == Nᵐ, "Input `data` has $Nᵐ modes, which is not a perfect square")
    check(Nᵗ == length(t), "Input `data` has $Nᵗ samples, but input `t` has $(length(t))")
    check(
        Nᵈ == ncomponents(dc),
        "Input `data` has $Nᵈ components, but `dc` has $(ncomponents(dc))",
    )
    check(
        length(αᵢₙ) ≤ Nᵐ, "Input `αᵢₙ` has $(length(αᵢₙ)) modes, but expected at most $Nᵐ"
    )
    ℓₘₐₓ = L - 1
    if !isnothing(ℓₘₐₓ₀)
        check(0 ≤ ℓₘₐₓ₀ ≤ ℓₘₐₓ, "ℓₘₐₓ₀=$ℓₘₐₓ₀ must be between 0 and ℓₘₐₓ=$ℓₘₐₓ")
        check(
            all(iszero, @view data[((ℓₘₐₓ₀ + 1) ^ 2 + 1):end, :, :]),
            "`data` contains nonzero modes above the declared input band limit " *
            "ℓₘₐₓ₀=$ℓₘₐₓ₀",
        )
    end
    N₀ = isnothing(ℓₘₐₓ₀) ? Nᵐ : (ℓₘₐₓ₀ + 1)^2
    return (Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀)
end

"""
    transform_pixel!(d′ᵢ, dᵢ, d̈ᵢ, cubic_spline_cache, t, t′,
                     κ⁻¹ᵢ, αₚᵢ, ðt′╱2κₚ₀ᵢ, ðt′╱2κₚ₁ᵢ, ð²αₚᵢ, dc)

Fill `d′ᵢ` (`Nᵈ × Nᵗ′`) with the BMS-transformed time series of a single pixel, given the
pixel's rest-frame time series `dᵢ` (`Nᵈ × Nᵗ`, sampled at `t`) and a scratch buffer `d̈ᵢ`
of the same size as `dᵢ` for the spline second derivatives.  This is the per-pixel kernel
shared by [`transform!`](@ref), [`transform_objective`](@ref), and
[`pixel_waveform`](@ref): a Thomas forward sweep for the natural cubic spline, followed by
a fused backward sweep that interpolates each output time `t′[j′] κ⁻¹ᵢ + αₚᵢ` and applies
the component-mixing transformation laws via [`mix_components!`](@ref).

The scalar arguments are this pixel's inverse conformal factor `κ⁻¹ᵢ`, supertranslation
value `αₚᵢ`, the two parts of the null-rotation parameter from [`compute_ðt′╱2κ`](@ref)
(`ðt′╱2κₚ₀ᵢ` constant, `ðt′╱2κₚ₁ᵢ` the coefficient of time), and the second
eth-derivative `ð²αₚᵢ`.

Everything is generic over element types: `dᵢ`, `d̈ᵢ`, and `d′ᵢ` may hold complex dual
numbers, and the scalars may be dual, while `cubic_spline_cache`, `t`, and `t′` remain
primal.  The output grid `t′` may have any length ≥ 1 — in particular it need *not* match
`length(t)`; the caller guarantees (via [`validate_t′`](@ref)) that every
`t′[j′] κ⁻¹ᵢ + αₚᵢ` lies within the span of `t` (up to roundoff; a τ ≈ 0⁻ excursion below
`t[1]` is extrapolated with error O(τ³)).  Returns `d′ᵢ`.
"""
function transform_pixel!(
    d′ᵢ,
    dᵢ,
    d̈ᵢ,
    cubic_spline_cache::CubicSplineCache,
    t,
    t′,
    κ⁻¹ᵢ,
    αₚᵢ,
    ðt′╱2κₚ₀ᵢ,
    ðt′╱2κₚ₁ᵢ,
    ð²αₚᵢ,
    dc::DataComponents,
)
    Nᵈ, Nᵗ = size(dᵢ)
    Nᵗ′ = length(t′)

    # `d̈` forward sweep (Thomas algorithm, natural BC: d̈[1]=d̈[Nᵗ]=0)
    @inbounds let
        @simd ivdep for k ∈ 1:Nᵈ
            d̈ᵢ[k, 1] = 0
        end
        @simd ivdep for k ∈ 1:Nᵈ
            r =
                6 * (
                    cubic_spline_cache.h⁻¹[2] * (dᵢ[k, 3] - dᵢ[k, 2]) -
                    cubic_spline_cache.h⁻¹[1] * (dᵢ[k, 2] - dᵢ[k, 1])
                )
            d̈ᵢ[k, 2] = r * cubic_spline_cache.u⁻¹[1]
        end
        for j ∈ 3:(Nᵗ - 1)
            @simd ivdep for k ∈ 1:Nᵈ
                r =
                    6 * (
                        cubic_spline_cache.h⁻¹[j] * (dᵢ[k, j + 1] - dᵢ[k, j]) -
                        cubic_spline_cache.h⁻¹[j - 1] * (dᵢ[k, j] - dᵢ[k, j - 1])
                    )
                d̈ᵢ[k, j] =
                    (r - cubic_spline_cache.h[j - 1] * d̈ᵢ[k, j - 1]) *
                    cubic_spline_cache.u⁻¹[j - 1]
            end
        end
        @simd ivdep for k ∈ 1:Nᵈ
            d̈ᵢ[k, Nᵗ] = 0
        end
    end

    # `d̈` backward sweep combined with interpolation and application of the BMS
    # transformation laws
    @inbounds let
        j′ = Nᵗ′
        tᵢⱼ′ = t′[j′] * κ⁻¹ᵢ + αₚᵢ  # original-frame time for output index j′
        for j ∈ (Nᵗ - 1):-1:1
            # Backward sweep step: d̈ᵢ[j] = z[j] − l[j-1]·d̈ᵢ[j+1]
            # (j=Nᵗ-1 and j=1 are natural-BC endpoints; no update needed)
            if 2 ≤ j ≤ Nᵗ-2
                @simd ivdep for k ∈ 1:Nᵈ
                    d̈ᵢ[k, j] -= cubic_spline_cache.l[j - 1] * d̈ᵢ[k, j + 1]
                end
            end

            # Evaluate all output times that fall in interval [t[j], t[j+1]].
            # The `|| j == 1` catches any j′ whose tᵢⱼ′ landed infinitesimally
            # below t[1] due to floating-point roundoff — τ will be ≈ 0⁻, and
            # the cubic extrapolation error is O(τ³), i.e., negligible.
            while j′ ≥ 1 && (tᵢⱼ′ ≥ t[j] || j == 1)
                let τ = tᵢⱼ′ - t[j]
                    @simd ivdep for k ∈ 1:Nᵈ
                        d′ᵢ[k, j′] = spline_eval(
                            dᵢ[k, j],
                            dᵢ[k, j + 1],
                            d̈ᵢ[k, j],
                            d̈ᵢ[k, j + 1],
                            cubic_spline_cache.h[j],
                            cubic_spline_cache.h⁻¹[j],
                            τ,
                        )
                    end
                end
                ðt′╱2κᵢⱼ = ðt′╱2κₚ₀ᵢ + tᵢⱼ′ * ðt′╱2κₚ₁ᵢ
                @views mix_components!(d′ᵢ[:, j′], κ⁻¹ᵢ, ðt′╱2κᵢⱼ, ð²αₚᵢ, dc)
                j′ -= 1
                if j′ ≥ 1
                    tᵢⱼ′ = t′[j′] * κ⁻¹ᵢ + αₚᵢ
                end
            end
        end
    end

    return d′ᵢ
end

"""
    transform!(data, t, v⃗, R, αᵢₙ, dc; t′=nothing, ℓₘₐₓ₀=nothing)

Transform the mode weights in `data` — sampled at times `t` — from the rest frame to the
BMS-transformed frame.  The BMS transformation is specified by the boost velocity `v⃗`, the
overall rotation `R`, the supertranslation `αᵢₙ`, and the `DataComponents` descriptor `dc`.
The transformation is performed in-place, modifying the input `data` array.

The complex `data` array is expected to have dimensions `(Nᵐ, Nᵗ, Nᵈ)`, where `Nᵐ` is the
number of modes, `Nᵗ` is the number of time samples, and `Nᵈ` is the number of data
components (e.g., strain and/or Newman-Penrose Weyl components).  The modes are expected to
be ordered by increasing `ℓ`, then by increasing `m` within each `ℓ`.  For data with spin
weight ``s ≠ 0``, the modes with ``ℓ < |s|`` are expected to be present, but will be
ignored.  The maximum `ℓ` value is determined by the size of the first dimension of `data`
as `ℓₘₐₓ = √Nᵐ - 1`.  The array must have complex type, with the underlying real type being
at least as wide as the types of the other inputs.

The `t` array is expected to have length `Nᵗ`, matching the second dimension of `data`.  The
`v⃗` and `R` inputs are expected to be of types `QuatVec` and `Rotor`.

The `αᵢₙ` array must be a complex vector, and must be ordered as described above for the
final dimension of `data`, though it may have a smaller `ℓₘₐₓ`.  (A larger `ℓₘₐₓ` cannot be
allowed, because the result would have higher angular dependence than `data`, which is not
possible since we are transforming `data` in-place.)  This array is expected to represent a
supertranslation in the rest frame, which is a real-valued function (with spin weight 0) on
the sphere.  The reality condition is that the mode weights satisfy ``α_{ℓ,-m} = (-1)^m
ᾱ_{ℓ,m}``.  This function will automatically impose this condition by averaging each mode
with its complex-conjugate partner.  This is done in a copy of the array for simplicity,
rather than being done in place.

The `dc` argument is a [`DataComponents`](@ref) value specifying which field components are
stored in `data`, in the order they appear along its third dimension.  Because of the
hierarchical nature of the BMS transformation, any Weyl component ``ψᵢ`` must be accompanied
by all higher-index components ``ψⱼ`` for ``j > i``.  Note that `DataComponents` includes a
sign indicating whether `data` represents data on ``ℐ⁺`` if `ℐ` = +1 or ``ℐ⁻`` if `ℐ` = -1.

The `dc` argument also contains the [`Conventions`](@ref) the data are expressed in, and
the transformation laws applied are the ones *native to that convention* — see the
["Convention dependence"](@ref convention_dependence_fields) section of the "BMS Action on
Fields" documentation page.  With the default conventions (SXS) every factor is the
identity and compiles away.

The optional keyword `t′` supplies the output time grid explicitly, instead of the default
grid constructed by [`compute_t′`](@ref).  It must be a real vector of the same length as
`t`, strictly increasing, and must lie within the range for which every pixel of the sphere
maps into the span of `t` (as validated by [`validate_t′`](@ref); an out-of-range grid
throws an `ArgumentError`).  This is chiefly useful in optimization loops over BMS
parameters, where a fixed output grid makes the results directly comparable across
iterations.

The optional keyword `ℓₘₐₓ₀` declares the *input* band limit: a promise that every mode
with ``ℓ > ℓₘₐₓ₀`` in `data` is exactly zero (as when band-limited data have been
zero-padded to a larger array, per the pad-first workflow of [Choosing
``ℓ_\\mathrm{max}``](@ref)).  The first (synthesis) stage then processes only the nonzero
modes, reducing its cost by roughly the factor `(ℓₘₐₓ₀+1)²/Nᵐ`.  The promise is verified,
and an `ArgumentError` is thrown if any declared-zero mode is nonzero.

Returns `(data, t′)`, where `t′` is the output time grid (the supplied one, if given).

The transformation is differentiable end-to-end with ForwardDiff with respect to the BMS
parameters, provided the output grid is held fixed with `t′` and ForwardDiff is loaded
(activating a package extension); see the [Automatic Differentiation](@ref)
documentation page.  Note, however, that differentiating through `transform!` requires a
dual-typed copy of the entire `data` array; when the goal is to optimize the BMS
parameters against a fixed target waveform, the memory-lean
[`transform_objective`](@ref) computes the L² objective directly, leaving `data` primal
and untouched.

"""
function transform!(
    data::Array{Complex{T1}},
    t::Vector{T2},
    v⃗::QuatVec{T3},
    R::Rotor{T4},
    αᵢₙ::Vector{Complex{T5}},
    dc::DataComponents{C,I};
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real,C,I}
    # Use this `let` block to ensure that we don't accidentally use `T` below, because that
    # could lead to type instability.
    let T = promote_type(T1, T2, T3, T4, T5)
        if T != T1
            throw(
                ArgumentError(
                    "\nInput `data` type $T1 does not match common input type $T.\n" *
                    "Because `transform!` modifies `data` in place, its type must be\n" *
                    "compatible with all the other input types:\n" *
                    "  - `t` has element type $T2\n" *
                    "  - `v⃗` has element type $T3\n" *
                    "  - `R` has element type $T4\n" *
                    "  - `α` has element type $T5\n",
                ),
            )
        end
    end

    # Check that the input data has the expected dimensions and properties.
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    Nᵖ = Nᵐ
    block_size = max(1, min(Nᵗ, cachesize_L2 ÷ (Nᵐ * sizeof(Complex{T1}))))

    # The time-law sign in t′ = κ(t − c_α α) comes from the data conventions; the remaining
    # convention factors (dyad rescaling of the mixing parameter and the F_σ/F_h shift
    # factors) are computed in `mix_components!`, which reads them from `dc` itself.
    c_α = dc.conventions.c_α

    ###
    ### Stage 0: Precompute various quantities needed for the transformation
    ###

    β = absvec(v⃗)
    γ = 1 / √(1 - β^2)
    vˣ, vʸ, vᶻ = vec(v⃗)
    Λ = Boost(I * v⃗) * R

    # Compute uniformly spaced rotors that are simple to produce, but close to ideal for
    # sampling the sphere.  Use spin weight 0 to accommodate all fields on the same grid.
    # Normally, when using `SSHTDirect`, we would just let it choose the grid for us, but
    # this would usually depend on the spin weight, and we need the same grid for all
    # components.  Moreover, we actually want the grid to be uniformly spaced in the
    # transformed frame, which means that we have to evaluate on a non-uniform grid in the
    # rest frame.
    # The grid is parameter-independent, so build it (and, below, the analysis
    # factorizations) at the primal float type — see `primal_float`; this keeps AD dual
    # numbers out of the constant `qr`/`lu` factorizations, whose Householder steps would
    # otherwise turn the (identically zero) perturbations into NaNs.
    T4′ = primal_float(T4)
    R′ₚ = golden_ratio_spiral_rotors(0, ℓₘₐₓ, T4′)

    # That uniformly spaced grid will be as seen in the transformed frame; here we compute
    # the corresponding rotors in the rest frame, on which we will evaluate the input data.
    # This is the boosted or distorted grid.
    Tₚ = promote_type(Rotor{T4}, T3)
    Rₚ = similar(R′ₚ, Tₚ)
    # `ℐ = +1` (ℐ⁺) vs `-1` (ℐ⁻) selects the past-vs-future-cone direction map;
    # see the `aberration` docstring and the "Future and past null infinity" conventions.
    Polyester.@batch for i ∈ eachindex(Rₚ)
        Rₚ[i] = aberration(R′ₚ[i], Λ)
    end

    # Calculate the LU factorization of the tridiagonal matrix for cubic spline
    task_cubic_spline_cache = OhMyThreads.@spawn CubicSplineCache(t)

    # Impose the reality condition on the input supertranslation, and pad with zeros up to
    # `ℓₘₐₓ` if necessary.  We always impose reality on the modes, rather than just taking
    # the real part of the result after evaluation, because we also need ðα and ð²α, which
    # need to be consistent with the reality condition.  This is done in a separate thread
    # to overlap with the computation of the SSHTs.
    task_α = OhMyThreads.@spawn impose_reality(αᵢₙ, ℓₘₐₓ, c_α)

    # Construct the set of spin-spherical-harmonic transforms, for each spin weight.  Here
    # we use `OffsetVector` so that they can be indexed by their spin weight.
    𝒯 = OffsetVector(
        OhMyThreads.tmap(
            s -> ₛ𝐘(s, ℓₘₐₓ, basetype(Tₚ), Rₚ),
            Matrix{Complex{basetype(Tₚ)}},
            -2:2;
            chunking=false,
        ),
        -3,
    )
    task_augmented_lu = OhMyThreads.@spawn begin
        # Build augmented square analysis matrices.  See the documentation page "Augmented
        # direct SSHT" for details.
        OffsetVector(
            map(-2:2) do s
                ₛY = ₛ𝐘(s, ℓₘₐₓ, T4′, R′ₚ)
                if s == 0
                    lu(ₛY)
                else
                    F = qr(ₛY)
                    Q = F.Q * Matrix{Bool}(LinearAlgebra.I, Nᵐ, Nᵐ)  # full Nᵐ×Nᵐ unitary
                    Q⊥ = Q[:, (Nᵐ - s ^ 2 + 1):end]  # Nᵐ × s² null-space columns
                    lu([Q⊥ ₛY])  # Nᵐ × Nᵐ, square
                end
            end,
            -3,
        )
    end

    # NOTE: From this point on, `α` will represent the corrected version that accounts for
    # `c_α`.  That is, we can now interpret `α` as being involved in the time translation as
    # t' = κ(t - α), rather than trying to keep that factor of c_α around.
    α = fetch(task_α)

    # Evaluate α on the boosted grid.  Make a copy because the 𝒯 act in place.
    task_αₚ = OhMyThreads.@spawn real.(𝒯[0] * copy(α))

    # Compute ðα, which is needed for ðt′/2κ in the Weyl transformation laws.
    # ð returns a full Nᵐ-element vector, but spin-1 modes start at ℓ=1, so skip the first
    # 1² = 1 leading zero entry before passing to 𝒯[1] (which expects Nᵐ − 1² modes).
    task_ðαₚ = OhMyThreads.@spawn 𝒯[1] * (ð(0, 0, ℓₘₐₓ, T5) * α)[2:end]

    # Compute ð²α, which is needed for σ or h data.
    # Same reasoning: spin-2 modes start at ℓ=2, so skip the first 2² = 4 leading zeros.
    task_ð²αₚ = OhMyThreads.@spawn 𝒯[2] * (ð(1, 0, ℓₘₐₓ, T5) * ð(0, 0, ℓₘₐₓ, T5) * α)[5:end]

    # Compute t′ — unless the caller supplied a grid, in which case just validate it here,
    # synchronously, so that a bad grid throws a plain ArgumentError rather than a
    # TaskFailedException.
    αₚ = fetch(task_αₚ)  # αₚ is also needed elsewhere, so fetch it before the task
    task_t′ = if isnothing(t′)
        # `compute_t′` returns the vector and the crossover time tᵪ, but we don't need the
        # latter, so just take the first element of the tuple.
        OhMyThreads.@spawn first(compute_t′(t, αₚ, Rₚ, v⃗, I))
    else
        # `transform!` writes the result back into `data` in place, so — unlike
        # `transform_objective` and `pixel_waveform` — its output grid must have exactly
        # as many samples as the input grid.
        if length(t′) != length(t)
            throw(
                ArgumentError(
                    "Supplied `t′` has $(length(t′)) samples, but `t` has $(length(t)).  " *
                    "Because `transform!` writes the result back into `data` in place, " *
                    "its output grid must have the same length as `t`.  For an output " *
                    "grid of a different length, see `transform_objective` and " *
                    "`pixel_waveform`.",
                ),
            )
        end
        validate_t′(t′, t, αₚ, Rₚ, v⃗, I)
        nothing
    end

    # Compute ðt′/2κ parts.
    ðαₚ = fetch(task_ðαₚ)
    task_ðt′╱2κₚ = OhMyThreads.@spawn compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)

    ###
    ### Stage 1: Evaluate all input data on the distorted grid
    ###

    # For each field component d, apply its SWSH synthesis matrix 𝒯[s] (Nᵐ×Nᵐ) to convert
    # mode weights → pixel values in-place.  The view `view(data, :, :, d)` is a fully
    # contiguous (Nᵐ × Nᵗ) matrix (leading dimension Nᵐ), so BLAS operates at full bandwidth
    # for this component slice.
    #
    # Each Julia task processes a contiguous column-chunk of time samples and further
    # subdivides it into (Nᵐ × block_size) sub-blocks sized to fit in L2 cache.  The
    # sub-block loop is necessary because each task's column-count (≈ Nᵗ/nthreads()) is
    # far larger than block_size — without it the task-local buffer would require
    # ~Nᵐ × Nᵗ/nthreads() × 16 B ≈ 69 MB, comparable to a full out-of-place copy.
    let
        cols = axes(data, 2)
        # If the caller declares an input band limit ℓₘₐₓ₀ (asserting that all higher
        # modes are zero), the synthesis GEMM only needs the first N₀ − s² columns of
        # 𝒯[s] — a column-contiguous view, so BLAS still runs at full speed — cutting
        # stage 1's cost from O(Nᵐ²Nᵗ) to O(NᵐN₀Nᵗ).
        for k ∈ 1:Nᵈ
            s = spin_weight(C[k])
            valid_modes = (s ^ 2 + 1):N₀  # skip leading ℓ < |s| entries
            data_k = view(data, :, :, k)  # (Nᵐ × Nᵗ), fully contiguous
            workspace = Matrix{Complex{T1}}(undef, length(valid_modes), block_size)
            𝒯ₛ = view(𝒯[s], :, 1:length(valid_modes))
            for sub_start ∈ cols[begin:block_size:end]
                sub = sub_start:min(sub_start + block_size - 1, cols[end])
                workspace_view = view(workspace, :, 1:length(sub))
                copyto!(workspace_view, view(data_k, valid_modes, sub))
                mul!(view(data_k, :, sub), 𝒯ₛ, workspace_view)
            end
        end
    end

    ###
    ### Stage 2: Interpolate to new slices and apply transformation laws
    ###

    cubic_spline_cache = fetch(task_cubic_spline_cache)
    # Both branches produce the same eltype that `compute_t′` would, keeping the loop below
    # type-stable; `convert` is a no-op when the supplied grid already has that type.
    t′ = let Tt = promote_type(eltype(t), eltype(αₚ), typeof(γ))
        isnothing(task_t′) ? convert(Vector{Tt}, t′) : fetch(task_t′)::Vector{Tt}
    end
    ðt′╱2κₚ = fetch(task_ðt′╱2κₚ)
    ð²αₚ = fetch(task_ð²αₚ)

    OhMyThreads.@tasks for i ∈ 1:Nᵖ
        OhMyThreads.@set scheduler = :static
        OhMyThreads.@set ntasks = nthreads()
        OhMyThreads.@local begin
            dᵢ = Matrix{Complex{T1}}(undef, Nᵈ, Nᵗ)
            d̈ᵢ = Matrix{Complex{T1}}(undef, Nᵈ, Nᵗ)
            d′ᵢ = Matrix{Complex{T1}}(undef, Nᵈ, Nᵗ)
        end

        κ⁻¹ᵢ = γ * (1 - I * v⃗dotk̂(Rₚ[i], vˣ, vʸ, vᶻ))

        # Copy pixel time series into the dᵢ buffer.  Note that tests comparing this
        # `permutedims!` approach to `LinearAlgebra.copy_transpose!` and to `.= transpose`
        # show this to be fastest and least allocating by up to ~2x, depending on Nᵈ.
        data_view = view(data, i, :, :)
        permutedims!(dᵢ, data_view, (2, 1))

        transform_pixel!(
            d′ᵢ,
            dᵢ,
            d̈ᵢ,
            cubic_spline_cache,
            t,
            t′,
            κ⁻¹ᵢ,
            αₚ[i],
            ðt′╱2κₚ[1, i],
            ðt′╱2κₚ[2, i],
            ð²αₚ[i],
            dc,
        )

        # Copy the transformed pixel time series back from the d′ᵢ buffer
        permutedims!(data_view, d′ᵢ, (2, 1))
    end  # OhMyThreads.@tasks

    ###
    ### Stage 3: Transform back to modes of the transformed data
    ###

    augmented_lu = fetch(task_augmented_lu)
    let
        for k ∈ 1:Nᵈ
            s = spin_weight(C[k])
            for sub_start ∈ 1:block_size:Nᵗ
                sub = sub_start:min(sub_start + block_size - 1, Nᵗ)
                ldiv!(augmented_lu[s], view(data, :, sub, k))
            end
        end
    end

    return data, t′
end

"""
    transform!(data, t, v⃗, R, αᵢₙ;
               data_components=nothing, ℐ=+1, conventions=Conventions(),
               t′=nothing, ℓₘₐₓ₀=nothing)

Keyword-argument convenience form.  See the main docstring for details.

The `data_components` argument may be a `DataComponents` value, a tuple of symbols such as
`(:ψ₄, :ψ₃)`, or a sequence of strings that indicate those symbols.  The strings are parsed
in a flexible way, so that, for example, `"psi4"`, `"Psi_4"`, and `"PSI₄"` all indicate the
same component `:ψ₄`.  Alternatively, if the argument is `nothing` (the default), the first
`Nᵈ` of `(:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀)` on ``ℐ⁺`` (or `(:h, :ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄)` on
``ℐ⁻``) will be chosen — though a warning will be issued.

The `ℐ` and `conventions` keywords are used only when `data_components` is *not* already a
`DataComponents` value (which contains its own).
"""
function transform!(
    data::Array{Complex{T1}},
    t::Vector{T2},
    v⃗::QuatVec{T3},
    R::Rotor{T4},
    αᵢₙ::Vector{Complex{T5}};
    data_components=nothing,
    ℐ::Int=+1,
    conventions::Conventions=Conventions(),
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real}
    Nᵈ = size(data, 3)
    dc = if data_components isa DataComponents
        data_components
    elseif isnothing(data_components)
        # The strain leads on both ends of null infinity; the Weyl tower orders differ.
        full_dc = ℐ == 1 ? (:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀) : (:h, :ψ₀, :ψ₁, :ψ₂, :ψ₃, :ψ₄)
        default_dc = full_dc[1:Nᵈ]
        @warn "Defaulting to data components $(default_dc).\n" *
            "Check that this is correct for your input data.\n" *
            "Consider passing a `DataComponents` value explicitly."
        DataComponents(default_dc...; ℐ, conventions)
    else
        DataComponents(data_components...; ℐ, conventions)
    end
    return transform!(data, t, v⃗, R, αᵢₙ, dc; t′, ℓₘₐₓ₀)
end

function transform!(
    data::Array{<:Complex},
    t::Vector{<:Real},
    v⃗::Vector{<:Real},
    R::Vector{<:Real},
    αᵢₙ::Vector{<:Complex};
    data_components=nothing,
    ℐ::Int=+1,
    conventions::Conventions=Conventions(),
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
)
    return transform!(
        data, t, QuatVec(v⃗), Rotor(R), αᵢₙ; data_components, ℐ, conventions, t′, ℓₘₐₓ₀
    )
end

"""
    transform!(data, t, g::BMS, dc::DataComponents; t′=nothing, ℓₘₐₓ₀=nothing)

Apply the BMS element `g` to `data` in place.  This convenience wrapper unpacks the boost
velocity, frame rotation, and supertranslation from `g` — via [`boost_velocity`](@ref),
[`frame_rotation`](@ref), and [`supertranslation`](@ref) — and forwards to the main
[`transform!`](@ref) method.  Because the accessors define `lorentz(g) = inv(Boost(v⃗) *
Lorentz(R))`, these parts have exactly the meaning the `v⃗` and `R` arguments have in the
main method, so `transform!(data, t, g, dc)` reproduces the action of `g` on the data.

The data are given on the null infinity singled out by `dc`'s `ℐ` type parameter, so `g` is
first re-expressed in that representation — and with the supertranslation sign given by the
data conventions' `c_α` — via the [`BMS`](@ref) conversion constructor, which accounts
exactly for whatever conventions (`A`, `I`) `g` was constructed with.  Returns `(data, t′)`,
as the main method does.
"""
# Generic-argument convenience: accept anything `QuatVec`/`Rotor` can convert — plain
# vectors, tuples, quaternions, slices of an optimizer's parameter vector, … — so callers
# never need to load Quaternionic or know its constructors.  The fully-typed main method
# is strictly more specific, so this cannot recurse.
function transform!(
    data::Array{<:Complex},
    t::Vector{<:Real},
    v⃗,
    R,
    αᵢₙ::Vector{<:Complex},
    dc::DataComponents;
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
)
    return transform!(data, t, QuatVec(v⃗), Rotor(R), αᵢₙ, dc; t′, ℓₘₐₓ₀)
end

function transform!(
    data::Array{<:Complex},
    t::Vector{<:Real},
    g::BMS,
    dc::DataComponents{C,I};
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {C,I}
    # Re-express `g` with the supertranslation sign matching the data conventions' `c_α`,
    # since the main method will interpret α through that sign.
    gᴵ = BMS(g; c_α=dc.conventions.c_α * 1, ℐ=I)
    return transform!(
        data, t, boost_velocity(gᴵ), frame_rotation(gᴵ), supertranslation(gᴵ), dc; t′, ℓₘₐₓ₀
    )
end

@doc raw"""
    pixel_waveform(data, t, dc, t′)

Synthesize the mode weights in `data` onto the standard pixel grid and interpolate to the
time samples `t′`, returning an array of pixel values of size `(Nᵖ, length(t′), Nᵈ)` with
`Nᵖ = Nᵐ`.  This is the *un*-transformed waveform, evaluated exactly where
[`transform_objective`](@ref) evaluates the transformed one: pixel `i` corresponds to the
direction of the `i`-th rotor of `golden_ratio_spiral_rotors(0, ℓₘₐₓ, T)`, and the values
are spline-interpolated from `t` to `t′` with the same natural cubic spline used by
[`transform!`](@ref).

The `data`, `t`, and `dc` arguments are exactly as in [`transform!`](@ref); in particular
`data` must have size `(Nᵐ, length(t), Nᵈ)`, with `Nᵐ` a perfect square setting
`ℓₘₐₓ = √Nᵐ - 1`.  If the target waveform is known at a smaller band limit than the data
to be transformed, zero-pad its modes up to the same `Nᵐ` before calling this function
(the pad-first workflow of [Choosing ``ℓ_\mathrm{max}``](@ref)) — the pixel grid, and
hence the pixel correspondence with `transform_objective`, is set by `Nᵐ`.

The output grid `t′` must be strictly increasing and contained in `[t[begin], t[end]]`,
but — unlike the `t′` keyword of `transform!` — it may have **any** length ≥ 1.  In
particular it may be much coarser than `t`, which makes the optimization objective
correspondingly cheaper.

This function is entirely independent of any BMS parameters, so in an optimization loop
it should be called **once**, up front, to build the `target` argument of
[`transform_objective`](@ref); every optimizer iteration then reuses the result.
"""
function pixel_waveform(
    data::Array{Complex{T1},3}, t::Vector{T2}, dc::DataComponents{C,I}, t′::Vector{<:Real}
) where {T1<:Real,T2<:Real,C,I}
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, _ = validate_transform_inputs(
        data, t, zero(QuatVec{T2}), Complex{T1}[], dc, nothing
    )
    Nᵖ = Nᵐ
    Nᵗ′ = length(t′)
    # Validate `t′` against the trivial transformation (no boost, no supertranslation),
    # for which the valid span is exactly [t[begin], t[end]].
    validate_t′(t′, t, zeros(T2, 1), [one(Rotor{T2})], zero(QuatVec{T2}), I)

    TP = promote_type(T1, T2, eltype(t′))
    TP′ = primal_float(TP)

    # Synthesize modes → pixels on the uniform grid, one GEMM per unique spin weight,
    # exactly as in stage 1 of `transform!` — but out of place, into `P`.
    R′ₚ = golden_ratio_spiral_rotors(0, ℓₘₐₓ, TP′)
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    P = Array{Complex{TP},3}(undef, Nᵖ, Nᵗ, Nᵈ)
    for s ∈ unique(spins)
        # ₛ𝐘 columns cover modes with ℓ ≥ |s|; the leading ℓ < |s| rows of `data` are
        # ignored, as in `transform!`.
        ₛY = ₛ𝐘(s, ℓₘₐₓ, TP′, R′ₚ)
        for k ∈ findall(==(s), spins)
            mul!(view(P, :, :, k), ₛY, view(data, (s ^ 2 + 1):Nᵐ, :, k))
        end
    end

    # Interpolate each pixel's time series from `t` to `t′` using the shared per-pixel
    # kernel with neutral transformation parameters (κ⁻¹ = 1, αₚ = ðt′/2κ = ð²α = 0),
    # for which every mixing law reduces to the identity.
    cubic_spline_cache = CubicSplineCache(t)
    out = Array{Complex{TP},3}(undef, Nᵖ, Nᵗ′, Nᵈ)
    OhMyThreads.@tasks for i ∈ 1:Nᵖ
        OhMyThreads.@set scheduler = :static
        OhMyThreads.@set ntasks = nthreads()
        OhMyThreads.@local begin
            dᵢ = Matrix{Complex{TP}}(undef, Nᵈ, Nᵗ)
            d̈ᵢ = Matrix{Complex{TP}}(undef, Nᵈ, Nᵗ)
            d′ᵢ = Matrix{Complex{TP}}(undef, Nᵈ, Nᵗ′)
        end
        permutedims!(dᵢ, view(P, i, :, :), (2, 1))
        transform_pixel!(
            d′ᵢ,
            dᵢ,
            d̈ᵢ,
            cubic_spline_cache,
            t,
            t′,
            one(TP),
            zero(TP),
            zero(Complex{TP}),
            zero(Complex{TP}),
            zero(Complex{TP}),
            dc,
        )
        permutedims!(view(out, i, :, :), d′ᵢ, (2, 1))
    end
    return out
end

@doc raw"""
    transform_objective(data, t, v⃗, R, αᵢₙ, dc, target; t′, weights=nothing, ℓₘₐₓ₀=nothing)

Compute the weighted L² difference between the BMS-transformed `data` and a fixed
`target`, as a single real scalar:

```math
\mathrm{obj} = ∑_i w_i ∑_j ∑_k \left| d′_{ijk} - \mathrm{target}_{ijk} \right|^2,
```

where ``d′_{ijk}`` is the transformed waveform at pixel `i`, output time `t′[j]`, and
component `k` — exactly the pixel values that [`transform!`](@ref) would compute in its
interpolation stage, before its final analysis back to modes.  This is the natural
objective for optimizing the BMS parameters `(v⃗, R, αᵢₙ)` to align `data` with a target
waveform, and it is designed for use with
[ForwardDiff](https://juliadiff.org/ForwardDiff.jl/): see the [Automatic
Differentiation](@ref) documentation page for a complete recipe.

The `data`, `t`, `v⃗`, `R`, `αᵢₙ`, `dc`, and `ℓₘₐₓ₀` arguments are exactly as in
[`transform!`](@ref) — except that `data` is **never modified**, and its element type
need *not* be widened to match the parameter types.  In particular, when `v⃗`, `R`, or
`αᵢₙ` carry ForwardDiff dual numbers, `data` stays primal and is simply read; no
dual-typed copy of the (typically very large) data array is ever made, no dual-typed
synthesis matrices are built, and the analysis stage of `transform!` is skipped entirely.
Dual numbers are confined to per-pixel scalars and a few task-local buffers of size
`O(Nᵈ × Nᵗ)`, so the memory overhead of differentiation is negligible, and `data` can be
reused across optimizer iterations without copying.  A `BMS`-valued transformation may be
passed in place of `(v⃗, R, αᵢₙ)`, and `v⃗` and `R` may be anything the `QuatVec` and
`Rotor` constructors accept — e.g., plain vectors, such as slices of an optimizer's
parameter vector.

The keyword `t′` is **required**: it is the fixed output time grid on which the
transformed waveform is compared to `target`, and holding it fixed across parameter
values is what makes objective values (and derivatives) comparable across optimizer
iterations.  Use [`compute_t′(t, β, δt)`](@ref compute_t′(::Any, ::Real, ::Real)) to
construct a grid valid for every transformation the optimizer may try.  Unlike the `t′`
keyword of `transform!`, it may have **any** length ≥ 1 — a coarse grid makes the
objective correspondingly cheaper — but it must be strictly increasing and remain within
the valid span ([`validate_t′`](@ref)) at every parameter value.

The `target` array must have size `(Nᵖ, length(t′), Nᵈ)`, with `Nᵖ = Nᵐ`, holding the
target waveform's *pixel* values on the standard grid — build it once with
[`pixel_waveform`](@ref) (zero-padding the target's modes up to `Nᵐ` first, if needed).
The optional `weights` vector supplies the per-pixel weights ``w_i`` (length `Nᵖ`,
default: uniform weight 1); any time- or component-dependent weighting can be folded into
the preprocessing that builds `target`.

Note that this is a *pixel-domain* L² norm.  Because the pixel grid has exactly as many
points as modes and the synthesis matrices are full rank, it vanishes exactly when the
mode-domain difference vanishes and defines an equivalent norm — but it is **not**
numerically equal to `sum(abs2, ...)` over mode weights, so objective values are
comparable to each other, not to mode-space norms.

Two fine points.  First, the spline-interval selection inside the interpolation is
discontinuous on a measure-zero set of parameter values (where a transformed time
`t′[j] κ⁻¹ + αₚ` crosses a knot of `t`); the objective itself is continuous there, and
derivatives are unaffected away from those isolated points.  Second, the result is
deterministic: for a fixed number of threads, repeated calls with identical inputs return
bitwise-identical values.

Returns a real scalar of the common promoted type of the inputs — a ForwardDiff dual when
the parameters carry duals.
"""
function transform_objective(
    data::Array{Complex{T1},3},
    t::Vector{T2},
    v⃗::QuatVec{T3},
    R::Rotor{T4},
    αᵢₙ::Vector{Complex{T5}},
    dc::DataComponents{C,I},
    target::AbstractArray{<:Complex,3};
    t′::Vector{<:Real},
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real,C,I}
    ###
    ### Stage A: Validation and parameter-independent setup
    ###

    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    Nᵖ = Nᵐ
    Nᵗ′ = length(t′)
    if size(target) != (Nᵖ, Nᵗ′, Nᵈ)
        throw(
            ArgumentError(
                "Input `target` has size $(size(target)), but expected " *
                "(Nᵖ, length(t′), Nᵈ) = $((Nᵖ, Nᵗ′, Nᵈ))",
            ),
        )
    end
    if !isnothing(weights) && length(weights) != Nᵖ
        throw(
            ArgumentError(
                "Input `weights` has length $(length(weights)), but expected Nᵖ=$Nᵖ"
            ),
        )
    end

    # TΘ is the (possibly dual) common type of the BMS parameters; TR is the common type
    # of everything, in which the accumulation runs.
    TΘ = promote_type(T3, T4, T5)
    TR = promote_type(T1, T2, TΘ)

    c_α = dc.conventions.c_α
    β = absvec(v⃗)
    γ = 1 / √(1 - β^2)
    vˣ, vʸ, vᶻ = vec(v⃗)
    Λ = Boost(I * v⃗) * R

    # The uniform output pixel grid is parameter-independent, so build it at the primal
    # float type (see `primal_float`), as in `transform!`.  The same goes for the ð
    # operator matrices, which depend only on the band limit.
    R′ₚ = golden_ratio_spiral_rotors(0, ℓₘₐₓ, primal_float(T4))
    T5′ = primal_float(T5)

    # As in `transform!`: from here on, `α` is the c_α-corrected supertranslation, so it
    # enters the time law as t′ = κ(t − α).  `ðα` and `ð²α` are full Nᵐ-length mode
    # vectors whose leading (ℓ < s) entries are ignored below.
    α = impose_reality(αᵢₙ, ℓₘₐₓ, c_α)
    ðα = ð(0, 0, ℓₘₐₓ, T5′) * α
    ð²α = ð(1, 0, ℓₘₐₓ, T5′) * ðα

    cubic_spline_cache = CubicSplineCache(t)

    ###
    ### Stage B: Per-pixel prepass — distorted grid, α and its derivatives on it
    ###

    # These per-pixel values must all exist before the main loop, because the grid
    # validation and `compute_ðt′╱2κ` consume them whole.
    Tₚ = promote_type(Rotor{T4}, T3)
    Rₚ = Vector{Tₚ}(undef, Nᵖ)
    αₚ = Vector{TΘ}(undef, Nᵖ)
    ðαₚ = Vector{Complex{TΘ}}(undef, Nᵖ)
    ð²αₚ = Vector{Complex{TΘ}}(undef, Nᵖ)
    OhMyThreads.@tasks for i ∈ 1:Nᵖ
        OhMyThreads.@set scheduler = :static
        OhMyThreads.@set ntasks = nthreads()
        OhMyThreads.@local storage = sYlm_prep(ℓₘₐₓ, 2, TΘ)
        Rₚ[i] = aberration(R′ₚ[i], Λ)
        # Evaluate α, ðα, and ð²α at this pixel by dotting each mode vector with the
        # matching ₛYₗₘ row (spins 0, 1, and 2).  Each `sYlm_values!` call overwrites the
        # storage, so each dot product completes before the next spin is computed.  The
        # sums run over the modes with ℓ ≥ s; `α` is real on the sphere by construction.
        Y = sYlm_values!(storage, Rₚ[i], 0)
        αₚ[i] = real(pixel_dot(Y, α, 1, Nᵐ))
        Y = sYlm_values!(storage, Rₚ[i], 1)
        ðαₚ[i] = pixel_dot(Y, ðα, 2, Nᵐ)
        Y = sYlm_values!(storage, Rₚ[i], 2)
        ð²αₚ[i] = pixel_dot(Y, ð²α, 5, Nᵐ)
    end

    # Validate the output grid at these parameter values (synchronously, so a bad grid
    # throws a plain ArgumentError), then assemble the null-rotation parameter parts.
    validate_t′(t′, t, αₚ, Rₚ, v⃗, I)
    ðt′╱2κₚ = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)

    ###
    ### Stage C: Fused per-pixel synthesis → interpolation/mixing → accumulation
    ###

    # Components sharing a spin weight reuse a single ₛYₗₘ evaluation per pixel.
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    unique_spins = unique(spins)
    spin_components = [findall(==(s), spins) for s ∈ unique_spins]

    return OhMyThreads.tmapreduce(
        +, OhMyThreads.index_chunks(1:Nᵖ; n=nthreads()); chunking=false
    ) do chunk
        # Task-local storage and buffers, allocated once per chunk; these are the only
        # dual-typed arrays of the whole computation.
        storage = sYlm_prep(ℓₘₐₓ, 2, TΘ)
        dᵢ = Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ)
        d̈ᵢ = Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ)
        d′ᵢ = Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ′)
        acc = zero(TR)
        for i ∈ chunk
            # Synthesize this pixel's time series directly from the (primal, read-only)
            # mode data: dᵢ[k, j] = Σₘ Yₘ data[m, j, k], over the modes with ℓ ≥ |s| and
            # within any declared input band limit (m ≤ N₀).  The inner loop over m walks
            # a contiguous column of `data`.
            for (s, ks) ∈ zip(unique_spins, spin_components)
                Y = sYlm_values!(storage, Rₚ[i], s)
                for k ∈ ks
                    @inbounds for j ∈ 1:Nᵗ
                        dᵢ[k, j] = pixel_dot(Y, view(data, :, j, k), s ^ 2 + 1, N₀)
                    end
                end
            end
            κ⁻¹ᵢ = γ * (1 - I * v⃗dotk̂(Rₚ[i], vˣ, vʸ, vᶻ))
            transform_pixel!(
                d′ᵢ,
                dᵢ,
                d̈ᵢ,
                cubic_spline_cache,
                t,
                t′,
                κ⁻¹ᵢ,
                αₚ[i],
                ðt′╱2κₚ[1, i],
                ðt′╱2κₚ[2, i],
                ð²αₚ[i],
                dc,
            )
            objᵢ = zero(TR)
            @inbounds for j ∈ 1:Nᵗ′, k ∈ 1:Nᵈ
                objᵢ += abs2(d′ᵢ[k, j] - target[i, j, k])
            end
            wᵢ = isnothing(weights) ? one(TR) : TR(weights[i])
            acc += wᵢ * objᵢ
        end
        acc
    end
end

"""
    pixel_dot(Y, f, m₁, m₂)

Contract the ₛYₗₘ row `Y` with the mode vector `f` over the index range `m₁:m₂` — no
conjugation, matching the synthesis convention `pixels = ₛ𝐘 modes`.  Both arguments may
hold dual numbers.
"""
@inline function pixel_dot(Y, f, m₁, m₂)
    acc = zero(promote_type(eltype(Y), eltype(f)))
    @inbounds @simd for m ∈ m₁:m₂
        acc = muladd(Y[m], f[m], acc)
    end
    return acc
end

# Generic-argument convenience, exactly as for `transform!`: accept anything
# `QuatVec`/`Rotor` can convert.  The fully-typed main method is strictly more specific,
# so this cannot recurse.
function transform_objective(
    data::Array{<:Complex},
    t::Vector{<:Real},
    v⃗,
    R,
    αᵢₙ::Vector{<:Complex},
    dc::DataComponents,
    target::AbstractArray{<:Complex,3};
    t′::Vector{<:Real},
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
)
    return transform_objective(
        data, t, QuatVec(v⃗), Rotor(R), αᵢₙ, dc, target; t′, weights, ℓₘₐₓ₀
    )
end

"""
    transform_objective(data, t, g::BMS, dc, target; t′, weights=nothing, ℓₘₐₓ₀=nothing)

Convenience method taking a BMS element `g` in place of `(v⃗, R, αᵢₙ)`, unpacked exactly
as in the corresponding [`transform!`](@ref) method — including re-expressing `g` in the
representation and supertranslation sign of `dc`.  See the main
[`transform_objective`](@ref) docstring for details.
"""
function transform_objective(
    data::Array{<:Complex},
    t::Vector{<:Real},
    g::BMS,
    dc::DataComponents{C,I},
    target::AbstractArray{<:Complex,3};
    t′::Vector{<:Real},
    weights::Union{Nothing,AbstractVector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {C,I}
    # Re-express `g` with the supertranslation sign matching the data conventions' `c_α`,
    # since the main method will interpret α through that sign.
    gᴵ = BMS(g; c_α=dc.conventions.c_α * 1, ℐ=I)
    return transform_objective(
        data,
        t,
        boost_velocity(gᴵ),
        frame_rotation(gᴵ),
        supertranslation(gᴵ),
        dc,
        target;
        t′,
        weights,
        ℓₘₐₓ₀,
    )
end
