@doc raw"""
    ObjectivePlan(target, t, dc, t′; grid, weights, time_weights, component_weights,
                  normalize=false, ℓ_fit=nothing)

Everything about the objective of [`transform_objective`](@ref) that does not depend on the
BMS parameters, computed once and reused at every evaluation:

```math
Φ = \frac{1}{𝒩} ∑_j w_j ∑_i v_i ∑_k u_k
    \left| μ^{ki}_{j} - c^{ki}_{j} \right|^2,
```

where ``μ^{ki}_{j}`` is the transformed waveform at pixel `j`, output time `t′[i]`, and
component `k`, and ``c^{ki}_{j}`` is the target there.  The notation is that of the page
[The gradient of the objective](@ref objective_gradient).

The arguments `target`, `t`, and `dc` are the target waveform's mode weights, its time
samples, and its [`DataComponents`](@ref), as for [`transform!`](@ref); `t′` is the fixed
output grid on which the transformed data are compared with the target.  Use
[`compute_t′(t, β, δt)`](@ref compute_t′(::Any, ::Real, ::Real)) to build a grid valid for
every transformation an optimizer may try.  It may have any length ≥ 1.

The keywords are:

- `grid`: a pair `(R′ₚ, w)` of rotors and sphere weights ``w_j``.  The default is
  [`sphere_quadrature`](@ref Scri.sphere_quadrature)`(ℓₘₐₓ)`, at the target's band limit,
  which makes the sum over pixels an exact integral over the sphere when the transformed
  data ``μ`` are band-limited to that ``ℓₘₐₓ`` too, so that ``Φ`` is then exactly the
  ``L²`` norm.  A boosted field is not band-limited, so in general the sum is a
  spectrally accurate approximation of that norm rather than an exact one; a finer grid,
  such as `sphere_quadrature(L)` with a larger `L`, can be passed to reduce the error.
  With the default grid, data of a higher band limit than the target's are rejected,
  because the quadrature would not be exact even for them; zero-pad the target's modes
  to the data's band limit first.
- `weights`: per-pixel factors multiplying ``w_j``, for example to down-weight a region of
  the sphere.  The default is one.
- `time_weights`: the ``v_i``.  The default is the trapezoid rule on `t′`, so that the sum
  over time is an integral even on a non-uniform grid; a single time has weight 1.
- `component_weights`: the ``u_k``, one per component of `dc`.  The default is one.
- `normalize`: if `true`, ``𝒩`` is the same sum with ``μ`` set to zero, so that ``Φ`` is
  relative to the target's norm; otherwise ``𝒩 = 1``.
- `ℓ_fit`: if given, compare only the modes with ``ℓ ≤ ℓ_\mathrm{fit}`` of the difference.
  This leaves out the modes that the data and the target were padded with to make room
  for the transformation, and any high-``ℓ`` noise in either:

  ```math
  Φ = \frac{1}{𝒩} ∑_i v_i ∑_k u_k ∑_{ℓm:\, |s_k| ≤ ℓ ≤ ℓ_\mathrm{fit}}
      \left| Π^{ki}{}_{ℓm} \right|^2,
  \qquad
  Π^{ki}{}_{ℓm}
  = ∑_j w_j\, \overline{{}_{s_k}Y'^{ℓm}_{j}} \left( μ^{ki}_{j} - c^{ki}_{j} \right),
  ```

  where ``{}_{s}Y'^{ℓm}_{j}`` is the harmonic evaluated at pixel `j` of the grid.

  The projection is the quadrature's, so with the default grid it is exact when ``μ`` is
  band-limited to ``2ℓₘₐₓ - ℓ_\mathrm{fit}``; `ℓ_fit` must then be at most the target's
  band limit.  It must be at least the largest ``|s|`` among the components.  The
  projection takes one more pass over the pixels, so the objective and its gradient cost
  about 1.3 to 1.4 times as much with `ℓ_fit` as without it.

All weights must be finite and nonnegative.  The target's pixel values are computed with
[`pixel_waveform`](@ref) on the grid.
"""
struct ObjectivePlan{T<:AbstractFloat}
    R′ₚ::Vector{Rotor{T}}
    w::Vector{T}
    t′::Vector{T}
    v::Vector{T}
    u::Vector{T}
    c::Array{Complex{T},3}
    spins::Vector{Int}
    norm²::T
    normalize::Bool
    ℓ_fit::Int                     # -1 for the full objective
    ℓ_grid::Int                    # band limit of the default grid; -1 for a grid passed in
    Ȳ::Vector{Matrix{Complex{T}}}  # per component, wⱼ conj(ₛYₙ(R′ⱼ)) at [n, j]
    ĉ::Vector{Matrix{Complex{T}}}  # per component, the projection of `c` at [n, i]
end

function ObjectivePlan(
    target::Array{Complex{T1},3},
    t::AbstractVector{T2},
    dc::DataComponents{C},
    t′::AbstractVector{T3};
    grid=nothing,
    weights=nothing,
    time_weights=nothing,
    component_weights=nothing,
    normalize::Bool=false,
    ℓ_fit::Union{Nothing,Integer}=nothing,
) where {T1<:Real,T2<:Real,T3<:Real,C}
    T = float(promote_type(T1, T2, T3))
    T <: AbstractFloat ||
        throw(ArgumentError("An `ObjectivePlan` holds no derivatives; got element type $T"))
    check(cond, msg) = cond || throw(ArgumentError(msg))
    good(x) = all(y -> isfinite(y) && y ≥ 0, x)
    ℓₘₐₓ = isqrt(size(target, 1)) - 1
    R′ₚ, w = if isnothing(grid)
        sphere_quadrature(max(ℓₘₐₓ, 0), T)
    else
        Rg, wg = grid
        check(
            length(Rg) == length(wg),
            "`grid` has $(length(Rg)) rotors but $(length(wg)) weights",
        )
        convert(Vector{Rotor{T}}, collect(Rg)), convert(Vector{T}, collect(wg))
    end
    if !isnothing(weights)
        check(
            length(weights) == length(w),
            "`weights` has length $(length(weights)), " *
            "but the grid has $(length(w)) pixels",
        )
        w = w .* weights
    end
    check(good(w), "Pixel weights must be finite and nonnegative")
    t′ᵀ = convert(Vector{T}, t′)
    c = convert(Array{Complex{T},3}, pixel_waveform(target, t, dc, t′ᵀ, R′ₚ))
    v = isnothing(time_weights) ? trapezoid_weights(t′ᵀ) : convert(Vector{T}, time_weights)
    check(
        length(v) == length(t′ᵀ),
        "`time_weights` has length $(length(v)), " * "but `t′` has $(length(t′ᵀ)) samples",
    )
    check(good(v), "Time weights must be finite and nonnegative")
    u = if isnothing(component_weights)
        ones(T, length(C))
    else
        convert(Vector{T}, component_weights)
    end
    check(
        length(u) == length(C),
        "`component_weights` has length $(length(u)), " *
        "but `dc` has $(length(C)) components",
    )
    check(good(u), "Component weights must be finite and nonnegative")
    spins = Int[spin_weight(Val(s)) for s ∈ C]
    Ȳ, ĉ = Matrix{Complex{T}}[], Matrix{Complex{T}}[]
    norm² = zero(T)
    if isnothing(ℓ_fit)
        for k ∈ eachindex(u), i ∈ eachindex(v), j ∈ eachindex(w)
            norm² += w[j] * v[i] * u[k] * abs2(c[j, i, k])
        end
    else
        check(
            ℓ_fit ≥ maximum(abs, spins),
            "ℓ_fit=$ℓ_fit is below the largest |s| among the components",
        )
        check(
            !isnothing(grid) || ℓ_fit ≤ ℓₘₐₓ,
            "ℓ_fit=$ℓ_fit exceeds the band limit ℓₘₐₓ=$ℓₘₐₓ of the default grid",
        )
        for (k, s) ∈ enumerate(spins)
            push!(Ȳ, permutedims(conj.(sYlm_matrix(R′ₚ, ℓ_fit, s)) .* w))
            push!(ĉ, Ȳ[k] * view(c, :, :, k))
            for i ∈ eachindex(v)
                norm² += v[i] * u[k] * sum(abs2, view(ĉ[k], :, i))
            end
        end
    end
    check(
        !normalize || norm² > 0, "Cannot normalize by a target whose weighted norm is zero"
    )
    return ObjectivePlan{T}(
        R′ₚ,
        w,
        t′ᵀ,
        v,
        u,
        c,
        spins,
        norm²,
        normalize,
        isnothing(ℓ_fit) ? -1 : Int(ℓ_fit),
        isnothing(grid) ? max(ℓₘₐₓ, 0) : -1,
        Ȳ,
        ĉ,
    )
end

"""
    validate_plan(plan, dc, N₀)

Throw an `ArgumentError` unless `plan` can be used with data of the components `dc` whose
nonzero modes are the first `N₀`.  The components must have the plan's spin weights.  With
the default grid, the data's band limit must also not exceed the target's, because the grid
integrates exactly only up to that band limit; past it the objective is silently wrong.
"""
function validate_plan(plan::ObjectivePlan, dc::DataComponents{C}, N₀) where {C}
    spins = Int[spin_weight(Val(c)) for c ∈ C]
    spins == plan.spins || throw(
        ArgumentError(
            "`dc` has components of spin weights $spins, but the plan was built " *
            "for spin weights $(plan.spins)",
        ),
    )
    ℓ_data = isqrt(N₀) - 1
    plan.ℓ_grid < 0 ||
        ℓ_data ≤ plan.ℓ_grid ||
        throw(
            ArgumentError(
                "The data have band limit ℓₘₐₓ=$ℓ_data, but the plan's default grid is exact " *
                "only up to the target's band limit ℓₘₐₓ=$(plan.ℓ_grid).  Pad the target " *
                "to the data's ℓₘₐₓ before building the plan, pass `ℓₘₐₓ₀` if the data's " *
                "higher modes are zero, or pass a finer `grid` to the plan.",
            ),
        )
    return nothing
end

function Base.show(io::IO, p::ObjectivePlan{T}) where {T}
    return print(
        io,
        "ObjectivePlan{$T}($(length(p.R′ₚ)) pixels, $(length(p.t′)) times, ",
        "$(length(p.u)) components$(p.ℓ_fit ≥ 0 ? ", ℓ_fit=$(p.ℓ_fit)" : "")",
        "$(p.normalize ? ", normalized" : ""))",
    )
end

"""
    trapezoid_weights(t)

The weights of the trapezoid rule on the samples `t`, so that `sum(trapezoid_weights(t) .*
f.(t))` approximates the integral of `f` from `t[begin]` to `t[end]`, on any increasing
grid.  A single sample gets weight 1, so that the sum is just the value there.
"""
function trapezoid_weights(t::AbstractVector{T}) where {T<:Real}
    N = length(t)
    N == 1 && return ones(T, 1)
    v = Vector{T}(undef, N)
    v[1] = (t[2] - t[1]) / 2
    v[N] = (t[N] - t[N - 1]) / 2
    for j ∈ 2:(N - 1)
        v[j] = (t[j + 1] - t[j - 1]) / 2
    end
    return v
end
