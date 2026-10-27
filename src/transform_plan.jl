"""
    transform_grid(ℓₘₐₓ, T)

The rotors, of element type `T`, of the pixels on which [`transform!`](@ref) evaluates the
transformed data and analyzes them back into modes of band limit `ℓₘₐₓ`.  There are as many
pixels as modes, `(ℓₘₐₓ + 1)^2`, and one grid serves every spin weight, because the
transformation laws mix the components pixel by pixel.  This is also the default grid of
[`pixel_waveform`](@ref).

The pixels are the Leja points of spin weight 0 from SphericalFunctions' `leja_rotors`,
chosen so that the matrix of harmonics on them is well conditioned.  The augmented analysis
of nonzero spin weights (see [`augmented_factorization`](@ref Scri.augmented_factorization))
inherits that conditioning: at ``ℓₘₐₓ = 32``, the condition numbers of the analysis matrices
of spin weights 0, 1, and 2 are about 70, 110, and 80, where on the golden-ratio spiral they
were about 10⁵, 7 × 10⁴, and 10⁴.  Choosing the points costs less than the factorizations
that `transform!` makes on them, about 0.06 seconds at ``ℓₘₐₓ = 32``.
"""
transform_grid(ℓₘₐₓ, ::Type{T}) where {T} = leja_rotors(0, ℓₘₐₓ, T)

@doc raw"""
    TransformPlan(ℓₘₐₓ, dc, T=Float64)
    TransformPlan(data, dc)

The part of [`transform!`](@ref) that depends on no BMS parameter, computed once for reuse
across calls: the grid [`transform_grid`](@ref Scri.transform_grid)`(ℓₘₐₓ, T)` on which the
transformed data are analyzed back into modes, and the factorization of the analysis matrix
on that grid for each spin weight among the components of `dc`.

Pass it as `transform!(…; plan)` when transforming many times at one band limit, precision,
and set of spin weights, as in a scan over parameters.  Each call then skips those
factorizations, which dominate the cost of `transform!` at high ``ℓ_\mathrm{max}``.  `T` is
the primal floating-point type of the data (`Float64` for `ComplexF64` data, also when the
data carry ForwardDiff duals); the second form takes `ℓₘₐₓ` and `T` from `data`.  A plan
built for some components serves any components whose spin weights it covers.
"""
struct TransformPlan{T<:AbstractFloat,F}
    ℓₘₐₓ::Int
    R′ₚ::Vector{Rotor{T}}
    factorizations::Dict{Int,F}
end

function TransformPlan(
    ℓₘₐₓ::Integer, dc::DataComponents{C}, ::Type{T}=Float64
) where {C,T<:AbstractFloat}
    ℓₘₐₓ ≥ 2 || throw(ArgumentError("ℓₘₐₓ=$ℓₘₐₓ, but the transformation needs ℓₘₐₓ ≥ 2"))
    spins = unique(Int[spin_weight(Val(c)) for c ∈ C])
    R′ₚ = transform_grid(ℓₘₐₓ, T)
    factorizations = OhMyThreads.tmap(s -> augmented_factorization(R′ₚ, ℓₘₐₓ, s), spins)
    return TransformPlan{T,eltype(factorizations)}(
        ℓₘₐₓ, R′ₚ, Dict(zip(spins, factorizations))
    )
end

function TransformPlan(data::AbstractArray{Complex{T}}, dc::DataComponents) where {T<:Real}
    return TransformPlan(isqrt(size(data, 1)) - 1, dc, primal_float(T))
end

function Base.show(io::IO, p::TransformPlan{T}) where {T}
    return print(
        io,
        "TransformPlan{$T}(ℓₘₐₓ=$(p.ℓₘₐₓ), spin weights $(sort(collect(keys(p.factorizations)))))",
    )
end

"""
    check_plan(plan, ℓₘₐₓ, T, spins)

Throw an `ArgumentError` unless the [`TransformPlan`](@ref) `plan` was built for band limit
`ℓₘₐₓ` and primal type `T`, and covers every spin weight in `spins`.
"""
function check_plan(plan::TransformPlan{TP}, ℓₘₐₓ, ::Type{T}, spins) where {TP,T}
    plan.ℓₘₐₓ == ℓₘₐₓ || throw(
        ArgumentError(
            "The plan was built for ℓₘₐₓ=$(plan.ℓₘₐₓ), but the data have ℓₘₐₓ=$ℓₘₐₓ"
        ),
    )
    TP == T || throw(
        ArgumentError("The plan was built for element type $TP, but the data need $T")
    )
    missing = setdiff(spins, keys(plan.factorizations))
    isempty(missing) || throw(
        ArgumentError("The plan has no factorization for spin weights $(sort(missing))")
    )
    return plan
end
