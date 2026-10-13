"""
    validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)

Check the shared input contract of [`transform!`](@ref) and [`transform_objective`](@ref):
`v⃗` is finite and subluminal; `t` is finite, strictly increasing, and has at least 4
samples; `data` is a 3-dimensional array whose first dimension is a perfect square and whose
second and third dimensions match those of `t` and `dc`, respectively; the band limit is at
least 2; `αᵢₙ` is finite, has a perfect-square length, and fits within the band limit; and
any declared input band limit `ℓₘₐₓ₀` is in range and honored by `data`.

This function throws an `ArgumentError` on any violation.  Its callers run it before
spawning any task, so that a bad input never surfaces as a `TaskFailedException`.  If there
is no violation, it returns `(Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀)`, where `N₀` is the number of modes
within the declared input band limit (`Nᵐ` if none was declared).
"""
function validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    # These are explicit throws (not `@assert`) so they survive even when asserts are
    # disabled.
    check(cond, msg) = cond || throw(ArgumentError(msg))
    check(all(isfinite, vec(v⃗)), "Input `v⃗` = $v⃗ has a non-finite component")
    check(absvec(v⃗) < 1, "Input `v⃗` has magnitude $(absvec(v⃗)), but expected less than 1")
    check(length(t) ≥ 4, "Input `t` has only $(length(t)) samples, but expected at least 4")
    check(all(isfinite, t), "Input `t` has a non-finite element")
    check(issorted(t; lt=(≤)), "Input `t` must be strictly increasing")
    check(ndims(data) == 3, "Input `data` has $(ndims(data)) dimensions, but expected 3")
    Nᵐ, Nᵗ, Nᵈ = size(data)
    L = isqrt(Nᵐ)
    check(L^2 == Nᵐ, "Input `data` has $Nᵐ modes, which is not a perfect square")
    check(Nᵗ == length(t), "Input `data` has $Nᵗ samples, but input `t` has $(length(t))")
    check(
        Nᵈ == ncomponents(dc),
        "Input `data` has $Nᵈ components, but `dc` has $(ncomponents(dc))",
    )
    check(all(isfinite, αᵢₙ), "Input `αᵢₙ` has a non-finite element")
    check(
        isqrt(length(αᵢₙ))^2 == length(αᵢₙ),
        "Input `αᵢₙ` has $(length(αᵢₙ)) modes, which is not a perfect square",
    )
    check(
        length(αᵢₙ) ≤ Nᵐ, "Input `αᵢₙ` has $(length(αᵢₙ)) modes, but expected at most $Nᵐ"
    )
    ℓₘₐₓ = L - 1
    # Both functions evaluate ð²α, of spin weight 2, whatever the data components are, and
    # the strain and the shear have spin weight ±2 themselves.
    check(ℓₘₐₓ ≥ 2, "Input `data` has ℓₘₐₓ=$ℓₘₐₓ, but the transformation needs ℓₘₐₓ ≥ 2")
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
    widen_parameters(v⃗, R, αᵢₙ, TP)

Convert the BMS parameters `v⃗`, `R`, and `αᵢₙ` so that each element type is at least as
wide as the real type `TP`, and return them as a tuple.  [`transform!`](@ref) and
[`transform_objective`](@ref) pass the primal float type of the data and the times as `TP`,
so that, for example, `Float32` parameters with `Float64` data have their geometry computed
in `Float64`.  A parameter that is already wide enough is returned unchanged; in particular,
the promotion turns no primal parameter into a dual number.
"""
function widen_parameters(
    v⃗::QuatVec{T3}, R::Rotor{T4}, αᵢₙ::Vector{Complex{T5}}, ::Type{TP}
) where {T3,T4,T5,TP}
    return (
        QuatVec{promote_type(T3, TP)}(v⃗),
        Rotor{promote_type(T4, TP)}(R),
        convert(Vector{Complex{promote_type(T5, TP)}}, αᵢₙ),
    )
end

@doc raw"""
    rotor_argument(R)

The rotor meant by a rotation argument given as something other than a `Rotor`, for the
generic-argument methods of [`transform!`](@ref) and [`transform_objective`](@ref).  Three
components, as a `QuatVec`, a vector, or a tuple, are a rotation vector ``\vec{r}``: the
right-handed rotation by the angle ``|\vec{r}|`` about the axis ``\vec{r}``, which is the
rotor ``\exp(\vec{r}/2)``.  So the zero vector is the identity, and slices of an
optimizer's parameter vector can start there.  Anything else goes to the `Rotor`
constructor, which reads four components as ``(w, x, y, z)`` and normalizes them.
"""
rotor_argument(R::Rotor) = R
rotor_argument(r⃗::QuatVec) = exp(r⃗ / 2)
function rotor_argument(R::Union{AbstractVector,Tuple})
    return length(R) == 3 ? rotor_argument(QuatVec(R...)) : Rotor(R)
end
rotor_argument(R) = Rotor(R)

"""
    augmented_factorization(R′ₚ, ℓₘₐₓ, s)

The LU factorization of the square analysis matrix of spin weight `s` on the grid `R′ₚ`,
which has as many points as there are modes up to `ℓₘₐₓ`.  For `s = 0` that matrix is the
synthesis matrix itself.  Otherwise the synthesis matrix lacks the `s²` columns of the modes
with `ℓ < |s|`, and is augmented with `s²` columns spanning the null space of its
transpose, from a full QR factorization; see the documentation page "Augmented Direct
SSHT".  Solving with the result gives the modes, preceded by `s²` coefficients that measure
what no band-limited field of that spin weight can represent.
"""
function augmented_factorization(R′ₚ, ℓₘₐₓ, s)
    ₛY = sYlm_matrix(R′ₚ, ℓₘₐₓ, s)
    s == 0 && return lu(ₛY)
    Nᵐ = size(ₛY, 1)
    Q = qr(ₛY).Q * Matrix{Bool}(LinearAlgebra.I, Nᵐ, Nᵐ)  # full Nᵐ×Nᵐ unitary
    Q⊥ = Q[:, (Nᵐ - s ^ 2 + 1):end]  # Nᵐ × s² null-space columns
    return lu([Q⊥ ₛY])  # Nᵐ × Nᵐ, square
end

"""
    supertranslation_on_grid(R′ₚ, Λ, αᵢₙ, c_α)

Return the distorted grid `Rₚ`, the rotors in the rest frame of the pixels `R′ₚ` of the
transformed frame under the Lorentz transformation `Λ`, together with the supertranslation
`αₚ`, its eth derivative `ðαₚ`, and its second eth derivative `ð²αₚ` at those rotors.
Either `Λ` or `αᵢₙ` may hold dual numbers, and the types of the results follow from
theirs: the rotors are those that [`aberration`](@ref) returns, and the values have the
real type of `Λ` promoted with that of `αᵢₙ`.

The reality condition is imposed on the modes, and they are multiplied by `c_α`, so that
from here on `α` enters the time law as `t′ = κ(t − α)`.  Each pixel's values are three dot
products of the mode vectors with the harmonics of spin weights 0, 1, and 2 there, from a
calculator stepped pixel by pixel.  The calculator works at the supertranslation's own band
limit (at least 2, the largest spin weight needed), which is usually far below that of the
data, so this costs much less than synthesis matrices at the data's band limit would.
"""
function supertranslation_on_grid(R′ₚ, Λ, αᵢₙ, c_α)
    TΛ, Tᵅ = real(basetype(Λ)), real(eltype(αᵢₙ))
    TΘ = promote_type(TΛ, Tᵅ)
    Tₚ = Rotor{TΛ}
    T5′ = primal_float(Tᵅ)
    ℓᵅ = max(2, isqrt(length(αᵢₙ)) - 1)
    Nᵅ = (ℓᵅ + 1)^2
    # `ðα` and `ð²α` are full mode vectors whose leading (ℓ < s) entries are zero.
    α = impose_reality(αᵢₙ, ℓᵅ, c_α)
    ðα = ð(0, 0, ℓᵅ, T5′) * α
    ð²α = ð(1, 0, ℓᵅ, T5′) * ðα
    Nᵖ = length(R′ₚ)
    Rₚ = Vector{Tₚ}(undef, Nᵖ)
    αₚ = Vector{TΘ}(undef, Nᵖ)
    ðαₚ = Vector{Complex{TΘ}}(undef, Nᵖ)
    ð²αₚ = Vector{Complex{TΘ}}(undef, Nᵖ)
    OhMyThreads.@tasks for i ∈ 1:Nᵖ
        OhMyThreads.@set scheduler = :static
        OhMyThreads.@set ntasks = nthreads()
        OhMyThreads.@local begin
            calc = sYlmCalculator(one(Tₚ), ℓᵅ, 0:2)
            Y = Matrix{Complex{basetype(Tₚ)}}(undef, Nᵅ, 3)
        end
        # `ℐ = +1` (ℐ⁺) vs `-1` (ℐ⁻) selects the past-vs-future-cone direction map; see the
        # `aberration` docstring and the "Future and past null infinity" conventions.
        Rₚ[i] = aberration(R′ₚ[i], Λ)
        # The columns of `Y` hold spins 0, 1, and 2; the sums run over the modes with ℓ ≥ s,
        # and `α` is real on the sphere by construction.
        sYlm_rows!(Y, calc, Rₚ[i])
        αₚ[i] = real(pixel_dot(view(Y, :, 1), α, 1, Nᵅ))
        ðαₚ[i] = pixel_dot(view(Y, :, 2), ðα, 2, Nᵅ)
        ð²αₚ[i] = pixel_dot(view(Y, :, 3), ð²α, 5, Nᵅ)
    end
    return Rₚ, αₚ, ðαₚ, ð²αₚ
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
    transform!(data, t, v⃗, R, αᵢₙ, dc; t′=nothing, ℓₘₐₓ₀=nothing, plan=nothing)

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
`v⃗` and `R` inputs are expected to be of types `QuatVec` and `Rotor`; other methods convert
anything those constructors accept, such as plain vectors.  A rotation given as three
components is read as a rotation vector, and four as a rotor; see [`rotor_argument`](@ref).

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

The optional keyword `ℓₘₐₓ₀` declares the *input* band limit: a
promise that every mode with ``ℓ > ℓₘₐₓ₀`` in `data` is exactly zero
(as when band-limited data have been zero-padded to a larger array,
per the pad-first workflow of [Choosing ``ℓ_\\mathrm{max}``](@ref
choosing_ell_max)).  The first (synthesis) stage then processes only
the nonzero modes, reducing its cost by roughly the factor
`(ℓₘₐₓ₀+1)²/Nᵐ`.  The promise is verified, and an `ArgumentError` is
thrown if any declared-zero mode is nonzero.

The optional keyword `plan` supplies a [`TransformPlan`](@ref) for the band limit, precision,
and spin weights of `data`, so that the uniform grid and its analysis factorizations, which
depend on no BMS parameter, are not rebuilt on every call.

Returns `(data, t′)`, where `t′` is the output time grid (the supplied one, if given).

The transformation is differentiable end-to-end with ForwardDiff with respect to the BMS
parameters, provided the output grid is held fixed with `t′`; see the [Automatic
Differentiation](@ref) documentation page.  Note, however, that differentiating through `transform!` requires a
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
    plan::Union{Nothing,TransformPlan}=nothing,
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

    # Compute the geometry at no lower precision than the data: with Float32 parameters and
    # Float64 data, the grid, the harmonics, and the factorizations would otherwise all be
    # Float32.  The branch is decided by the types alone, so it compiles away.
    let TP = promote_type(primal_float(T1), primal_float(T2))
        widened = (promote_type(T3, TP), promote_type(T4, TP), promote_type(T5, TP))
        if widened !== (T3, T4, T5)
            return transform!(
                data, t, widen_parameters(v⃗, R, αᵢₙ, TP)..., dc; t′, ℓₘₐₓ₀, plan
            )
        end
    end

    # Check that the input data has the expected dimensions and properties.
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    Nᵖ = Nᵐ
    cachesize_L2[] > 0 || read_cache_size!()
    block_size = max(1, min(Nᵗ, cachesize_L2[] ÷ (Nᵐ * sizeof(Complex{T1}))))

    # The time-law sign in t′ = κ(t − c_α α) comes from the data conventions; the remaining
    # convention factors (dyad rescaling of the mixing parameter and the F_σ/F_h shift
    # factors) are computed in `mix_components!`, which reads them from `dc` itself.
    c_α = dc.conventions.c_α

    ###
    ### Stage 0: Precompute various quantities needed for the transformation
    ###

    γ = lorentz_factor(v⃗)
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
    spins = unique(Int[spin_weight(Val(c)) for c ∈ C])
    if isnothing(plan)
        R′ₚ = golden_ratio_spiral_rotors(0, ℓₘₐₓ, T4′)
        # The analysis factorization on the uniform grid, for each spin weight that the
        # data have, and no others; it overlaps with the setup below.
        task_augmented_lu = OhMyThreads.@spawn Dict(
            s => augmented_factorization(R′ₚ, ℓₘₐₓ, s) for s ∈ spins
        )
    else
        check_plan(plan, ℓₘₐₓ, T4′, spins)
        R′ₚ = plan.R′ₚ
        task_augmented_lu = nothing
    end

    # Calculate the LU factorization of the tridiagonal matrix for cubic spline
    task_cubic_spline_cache = OhMyThreads.@spawn CubicSplineCache(t)

    # That uniformly spaced grid will be as seen in the transformed frame; here we compute
    # the corresponding rotors in the rest frame, on which we will evaluate the input data
    # (the boosted or distorted grid), together with the supertranslation and its eth
    # derivatives there.  See `supertranslation_on_grid`.
    Rₚ, αₚ, ðαₚ, ð²αₚ = supertranslation_on_grid(R′ₚ, Λ, αᵢₙ, c_α)
    Tₚ = eltype(Rₚ)

    # The synthesis matrix on the distorted grid, for each spin weight that the data have.
    𝒯 = Dict(
        zip(
            spins,
            OhMyThreads.tmap(
                s -> sYlm_matrix(Rₚ, ℓₘₐₓ, s),
                Matrix{Complex{basetype(Tₚ)}},
                spins;
                chunking=false,
            ),
        ),
    )

    # Compute t′ — or, if the caller supplied a grid, validate it.  Either way this runs
    # synchronously, so that a bad grid, or a transformation that leaves no complete slice,
    # throws a plain ArgumentError rather than a TaskFailedException.  Both branches produce
    # the same element type, keeping the loop below type-stable; `convert` is a no-op when
    # the supplied grid already has that type.
    t′ = let Tt = promote_type(eltype(t), eltype(αₚ), typeof(γ))
        if isnothing(t′)
            compute_t′(t, αₚ, Rₚ, v⃗, I)::Vector{Tt}
        else
            # `transform!` writes the result back into `data` in place, so — unlike
            # `transform_objective` and `pixel_waveform` — its output grid must have
            # exactly as many samples as the input grid.
            if length(t′) != length(t)
                throw(
                    ArgumentError(
                        "Supplied `t′` has $(length(t′)) samples, but `t` has " *
                        "$(length(t)).  Because `transform!` writes the result back " *
                        "into `data` in place, its output grid must have the same " *
                        "length as `t`.  For an output grid of a different length, see " *
                        "`transform_objective` and `pixel_waveform`.",
                    ),
                )
            end
            validate_t′(t′, t, αₚ, Rₚ, v⃗, I)
            convert(Vector{Tt}, t′)
        end
    end

    # Compute ðt′/2κ parts.
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
    ðt′╱2κₚ = fetch(task_ðt′╱2κₚ)

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

    augmented_lu =
        isnothing(task_augmented_lu) ? plan.factorizations : fetch(task_augmented_lu)
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
               t′=nothing, ℓₘₐₓ₀=nothing, plan=nothing)

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
    plan::Union{Nothing,TransformPlan}=nothing,
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
    return transform!(data, t, v⃗, R, αᵢₙ, dc; t′, ℓₘₐₓ₀, plan)
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
    plan::Union{Nothing,TransformPlan}=nothing,
)
    return transform!(
        data,
        t,
        QuatVec(v⃗),
        rotor_argument(R),
        αᵢₙ;
        data_components,
        ℐ,
        conventions,
        t′,
        ℓₘₐₓ₀,
        plan,
    )
end

# Generic-argument convenience: accept anything `QuatVec`/`Rotor` can convert — plain
# vectors, tuples, quaternions, slices of an optimizer's parameter vector, … — so callers
# never need to load Quaternionic or know its constructors.  A three-component rotation is
# a rotation vector (see `rotor_argument`).  The fully-typed main method is strictly more
# specific, so this cannot recurse.
function transform!(
    data::Array{<:Complex},
    t::Vector{<:Real},
    v⃗,
    R,
    αᵢₙ::Vector{<:Complex},
    dc::DataComponents;
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
    plan::Union{Nothing,TransformPlan}=nothing,
)
    return transform!(data, t, QuatVec(v⃗), rotor_argument(R), αᵢₙ, dc; t′, ℓₘₐₓ₀, plan)
end

"""
    transform!(data, t, g::BMS, dc::DataComponents; t′=nothing, ℓₘₐₓ₀=nothing, plan=nothing)

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
function transform!(
    data::Array{<:Complex},
    t::Vector{<:Real},
    g::BMS,
    dc::DataComponents{C,I};
    t′::Union{Nothing,Vector{<:Real}}=nothing,
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
    plan::Union{Nothing,TransformPlan}=nothing,
) where {C,I}
    # Re-express `g` with the supertranslation sign matching the data conventions' `c_α`,
    # since the main method will interpret α through that sign.
    gᴵ = BMS(g; c_α=dc.conventions.c_α * 1, ℐ=I)
    return transform!(
        data,
        t,
        boost_velocity(gᴵ),
        frame_rotation(gᴵ),
        supertranslation(gᴵ),
        dc;
        t′,
        ℓₘₐₓ₀,
        plan,
    )
end

@doc raw"""
    pixel_waveform(data, t, dc, t′, [R′ₚ])

Synthesize the mode weights in `data` at the rotors `R′ₚ` and interpolate to the time
samples `t′`, returning an array of pixel values of size `(length(R′ₚ), length(t′), Nᵈ)`.
This is the *un*-transformed waveform, evaluated on a grid in the frame of the transformed
data: pixel `i` is the direction of `R′ₚ[i]`, and the values are spline-interpolated from
`t` to `t′` with the same natural cubic spline used by [`transform!`](@ref).  Without
`R′ₚ`, the grid is `golden_ratio_spiral_rotors(0, ℓₘₐₓ, T)`, the `Nᵐ` pixels on which
`transform!` evaluates the transformed waveform.

The `data`, `t`, and `dc` arguments are exactly as in
[`transform!`](@ref); in particular `data` must have size `(Nᵐ,
length(t), Nᵈ)`, with `Nᵐ` a perfect square setting `ℓₘₐₓ = √Nᵐ - 1`.

The output grid `t′` must be strictly increasing and contained in `[t[begin], t[end]]`,
but — unlike the `t′` keyword of `transform!` — it may have **any** length ≥ 1.  In
particular it may be much coarser than `t`.
"""
function pixel_waveform(
    data::Array{Complex{T1},3}, t::Vector{T2}, dc::DataComponents, t′::Vector{<:Real}
) where {T1<:Real,T2<:Real}
    ℓₘₐₓ = isqrt(size(data, 1)) - 1
    T′ = primal_float(promote_type(T1, T2, eltype(t′)))
    return pixel_waveform(data, t, dc, t′, golden_ratio_spiral_rotors(0, ℓₘₐₓ, T′))
end
function pixel_waveform(
    data::Array{Complex{T1},3},
    t::Vector{T2},
    dc::DataComponents{C,I},
    t′::Vector{<:Real},
    R′ₚ::AbstractVector{<:Rotor},
) where {T1<:Real,T2<:Real,C,I}
    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, _ = validate_transform_inputs(
        data, t, zero(QuatVec{T2}), Complex{T1}[], dc, nothing
    )
    Nᵖ = length(R′ₚ)
    Nᵗ′ = length(t′)
    # Validate `t′` against the trivial transformation (no boost, no supertranslation),
    # for which the valid span is exactly [t[begin], t[end]].
    validate_t′(t′, t, zeros(T2, 1), [one(Rotor{T2})], zero(QuatVec{T2}), I)

    TP = promote_type(T1, T2, eltype(t′), basetype(eltype(R′ₚ)))

    # Synthesize modes → pixels on the grid, one GEMM per unique spin weight, as in stage 1
    # of `transform!` — but out of place, into `P`.
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    P = Array{Complex{TP},3}(undef, Nᵖ, Nᵗ, Nᵈ)
    for s ∈ unique(spins)
        # `sYlm_matrix` columns cover modes with ℓ ≥ |s|; the leading ℓ < |s| rows of
        # `data` are ignored, as in `transform!`.
        ₛY = sYlm_matrix(R′ₚ, ℓₘₐₓ, s)
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
    transform_objective(data, t, v⃗, R, αᵢₙ, dc, plan; ℓₘₐₓ₀=nothing)

Compute the objective of the [`ObjectivePlan`](@ref) `plan`, the weighted ``L²``
difference between the BMS-transformed `data` and the plan's target, as a single real
scalar:

```math
Φ = \frac{1}{𝒩} ∑_j w_j ∑_i v_i ∑_k u_k
    \left| {}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j} \right|^2,
```

where ``{}^{k}d'^{i}_{j}`` is the transformed waveform at the plan's pixel `j`, output time
`plan.t′[i]`, and component `k`, and the target ``c``, the weights ``w``, ``v``, and ``u``,
and the normalization ``𝒩`` are the plan's.  With the plan's default grid and time
weights, ``Φ`` is a quadrature of the integral of ``|d' - c|²`` over the
sphere and the time window, exact over the sphere for band-limited fields.  This is the
natural objective for optimizing the BMS parameters `(v⃗, R, αᵢₙ)` to align `data` with a
target waveform, and it is designed for use with
[ForwardDiff](https://juliadiff.org/ForwardDiff.jl/): see the
[Automatic Differentiation](@ref) documentation page for a complete recipe.

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
parameter vector.  A rotation given as three components is read as a rotation vector, and
four as a rotor; see [`rotor_argument`](@ref).

The plan's output grid `t′` is held fixed across parameter values, which is what makes
objective values (and derivatives) comparable across optimizer iterations.  It must remain
within the valid span ([`validate_t′`](@ref)) at every parameter value; build it with
[`compute_t′(t, β, δt)`](@ref compute_t′(::Any, ::Real, ::Real)) for the largest boost and
supertranslation the optimizer may try.  `dc` must have the same components as the plan.

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
    plan::ObjectivePlan{TQ};
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
) where {T1<:Real,T2<:Real,T3<:Real,T4<:Real,T5<:Real,C,I,TQ<:AbstractFloat}
    # As in `transform!`, compute the geometry at no lower precision than the data, the
    # times, and the plan.  The branch is decided by the types alone, so it compiles away.
    let TP = promote_type(primal_float(T1), primal_float(T2), TQ)
        widened = (promote_type(T3, TP), promote_type(T4, TP), promote_type(T5, TP))
        if widened !== (T3, T4, T5)
            return transform_objective(
                data, t, widen_parameters(v⃗, R, αᵢₙ, TP)..., dc, plan; ℓₘₐₓ₀
            )
        end
    end

    ###
    ### Stage A: Validation and parameter-independent setup
    ###

    Nᵐ, Nᵗ, Nᵈ, ℓₘₐₓ, N₀ = validate_transform_inputs(data, t, v⃗, αᵢₙ, dc, ℓₘₐₓ₀)
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    if spins != plan.spins
        throw(
            ArgumentError(
                "`dc` has components of spin weights $spins, but the plan was built " *
                "for spin weights $(plan.spins)",
            ),
        )
    end
    R′ₚ, w, t′, v, c, target = plan.R′ₚ, plan.w, plan.t′, plan.v, plan.c, plan.target
    Nᵖ = length(R′ₚ)
    Nᵗ′ = length(t′)

    # TΘ is the (possibly dual) common type of the BMS parameters; TR is the common type
    # of everything, in which the accumulation runs.
    TΘ = promote_type(T3, T4, T5)
    TR = promote_type(T1, T2, TΘ, TQ)

    c_α = dc.conventions.c_α
    γ = lorentz_factor(v⃗)
    vˣ, vʸ, vᶻ = vec(v⃗)
    Λ = Boost(I * v⃗) * R

    cubic_spline_cache = CubicSplineCache(t)

    ###
    ### Stage B: Per-pixel prepass — distorted grid, α and its derivatives on it
    ###

    # These per-pixel values must all exist before the main loop, because the grid
    # validation and `compute_ðt′╱2κ` consume them whole.
    Rₚ, αₚ, ðαₚ, ð²αₚ = supertranslation_on_grid(R′ₚ, Λ, αᵢₙ, c_α)
    Tₚ = eltype(Rₚ)

    # Validate the output grid at these parameter values (synchronously, so a bad grid
    # throws a plain ArgumentError), then assemble the null-rotation parameter parts.
    validate_t′(t′, t, αₚ, Rₚ, v⃗, I)
    ðt′╱2κₚ = compute_ðt′╱2κ(Rₚ, v⃗, αₚ, ðαₚ, I)

    ###
    ### Stage C: Fused per-pixel synthesis → interpolation/mixing → accumulation
    ###

    # Components sharing a spin weight reuse a single ₛYₗₘ evaluation per pixel, and one
    # calculator evaluates every spin weight in the range from a single recurrence.
    unique_spins = unique(spins)
    spin_components = [findall(==(s), spins) for s ∈ unique_spins]
    spin_range = minimum(spins):maximum(spins)

    # Fill `d′ᵢ` with the transformed time series of pixel `i`.  The pixel's rest-frame
    # time series is synthesized directly from the (primal, read-only) mode data:
    # dᵢ[k, j] = Σₘ Yₘ data[m, j, k], over the modes with ℓ ≥ |s| and within any declared
    # input band limit (m ≤ N₀); the inner loop over m walks a contiguous column of `data`.
    function transformed_pixel!(d′ᵢ, dᵢ, d̈ᵢ, calc, Y, i)
        sYlm_rows!(Y, calc, Rₚ[i])
        for (s, ks) ∈ zip(unique_spins, spin_components)
            Yˢ = view(Y, :, s - first(spin_range) + 1)
            for k ∈ ks
                @inbounds for j ∈ 1:Nᵗ
                    dᵢ[k, j] = pixel_dot(Yˢ, view(data, :, j, k), s ^ 2 + 1, N₀)
                end
            end
        end
        κ⁻¹ᵢ = γ * (1 - I * v⃗dotk̂(Rₚ[i], vˣ, vʸ, vᶻ))
        return transform_pixel!(
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
    end
    # Task-local storage and buffers, allocated once per chunk; these are the only
    # dual-typed arrays of the whole computation.
    function buffers()
        return (
            sYlmCalculator(one(Tₚ), ℓₘₐₓ, spin_range),
            Matrix{Complex{basetype(Tₚ)}}(undef, Nᵐ, length(spin_range)),
            Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ),
            Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ),
            Matrix{Complex{TR}}(undef, Nᵈ, Nᵗ′),
        )
    end
    chunks = OhMyThreads.index_chunks(1:Nᵖ; n=nthreads())

    Φ = if plan.ℓ_fit < 0
        OhMyThreads.tmapreduce(+, chunks; chunking=false) do chunk
            calc, Y, dᵢ, d̈ᵢ, d′ᵢ = buffers()
            acc = zero(TR)
            for i ∈ chunk
                transformed_pixel!(d′ᵢ, dᵢ, d̈ᵢ, calc, Y, i)
                objᵢ = zero(TR)
                @inbounds for j ∈ 1:Nᵗ′
                    objᵢⱼ = zero(TR)
                    for k ∈ 1:Nᵈ
                        objᵢⱼ += c[k] * abs2(d′ᵢ[k, j] - target[i, j, k])
                    end
                    objᵢ += v[j] * objᵢⱼ
                end
                acc += w[i] * objᵢ
            end
            acc
        end
    else
        # The band-limited objective: project the transformed waveform onto the modes with
        # ℓ ≤ ℓ_fit by the quadrature, Pₖ[n, j] = Σᵢ wᵢ conj(ₛYₙ(Rᵢ)) d′ᵢⱼₖ, whose weighted
        # conjugate harmonics the plan holds as the columns of Ȳ[k], and compare with the
        # target's projection.  The partial projections of the chunks add up.
        P = OhMyThreads.tmapreduce(+, chunks; chunking=false) do chunk
            calc, Y, dᵢ, d̈ᵢ, d′ᵢ = buffers()
            Pc = [zeros(Complex{TR}, size(plan.Ȳ[k], 1), Nᵗ′) for k ∈ 1:Nᵈ]
            for i ∈ chunk
                transformed_pixel!(d′ᵢ, dᵢ, d̈ᵢ, calc, Y, i)
                @inbounds for k ∈ 1:Nᵈ
                    Ȳₖ, Pₖ = plan.Ȳ[k], Pc[k]
                    for j ∈ 1:Nᵗ′
                        d = d′ᵢ[k, j]
                        for n ∈ axes(Ȳₖ, 1)
                            Pₖ[n, j] = muladd(Ȳₖ[n, i], d, Pₖ[n, j])
                        end
                    end
                end
            end
            Pc
        end
        sum(
            c[k] * sum(
                v[j] * sum(abs2, view(P[k], :, j) .- view(plan.T̂[k], :, j)) for j ∈ 1:Nᵗ′
            ) for k ∈ 1:Nᵈ
        )
    end
    return plan.normalize ? Φ / plan.norm² : Φ
end

"""
    sYlm_rows!(Y, calc, R)

Evaluate the spin-weighted spherical harmonics at the rotor `R` with the calculator `calc`,
and store them in the columns of `Y`, one column for each spin weight that `calc` serves,
in ascending order.  Each column is the ₛYₗₘ row of the synthesis matrix in the flat mode
ordering from ``ℓ = 0`` that the data use, so it has `(ℓₘₐₓ + 1)^2` entries, of which those
with ``ℓ < |s|`` are zero.  The element type of `R` must be that of `calc`.
"""
function sYlm_rows!(Y::AbstractMatrix, calc, R)
    set_R!(calc, R)
    for (ℓ, Yˡ) ∈ calc
        for (j, s) ∈ enumerate(SphericalFunctions.spins(calc)), m ∈ (-ℓ):ℓ
            Y[ℓ ^ 2 + ℓ + m + 1, j] = Yˡ[s, m]
        end
    end
    return Y
end

"""
    pixel_dot(Y, f, m₁, m₂)

Contract the ₛYₗₘ row `Y` with the mode vector `f` over the index range `m₁:m₂` — no
conjugation, matching the synthesis convention `pixels = sYlm_matrix(R⃗, ℓₘₐₓ, s) * modes`.
Both arguments may hold dual numbers.
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
    plan::ObjectivePlan;
    ℓₘₐₓ₀::Union{Nothing,Int}=nothing,
)
    return transform_objective(data, t, QuatVec(v⃗), rotor_argument(R), αᵢₙ, dc, plan; ℓₘₐₓ₀)
end

"""
    transform_objective(data, t, g::BMS, dc, plan; ℓₘₐₓ₀=nothing)

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
    plan::ObjectivePlan;
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
        plan;
        ℓₘₐₓ₀,
    )
end
