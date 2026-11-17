# The part of the comparison that runs on `main`: finite differences of an objective built
# from `transform!`.  Run by `run.jl`; see there.
#
# On `main`, `transform!` chooses its own output times, which move with the transformation.
# So an objective that compares the transformed waveform with a fixed target has to resample
# the transformed modes onto the target's times at every evaluation.  Three versions of the
# objective are timed:
#
#   - with no resampling at all, which is not a usable objective, but bounds the cost of one
#     from below;
#   - with a carefully written resampling, done the way `transform!` itself interpolates;
#   - with a plainly written resampling, as one might write it in ten minutes.
#
# Every timed loop is inside a function whose arguments have concrete types.  A loop over
# untyped global variables once made these timings five times too long.

using Scri
using LinearAlgebra: Tridiagonal, lu
include("problem.jl")
const QuatVec, Rotor = Scri.QuatVec, Scri.Rotor

# The parameterization of `bms_from_parameters` on the branch, which `main` does not have: a
# rotation vector, a rapidity vector, and the real degrees of freedom of the
# supertranslation's modes.
function bms(θ)
    r⃗ = QuatVec(θ[1], θ[2], θ[3])
    η = sqrt(θ[4]^2 + θ[5]^2 + θ[6]^2)
    f = η == 0 ? 1.0 : tanh(η) / η
    v⃗ = QuatVec(f * θ[4], f * θ[5], f * θ[6])
    a = view(θ, 7:length(θ))
    L = isqrt(length(a))
    α = zeros(ComplexF64, length(a))
    for ℓ ∈ 0:(L - 1)
        i₀ = ℓ^2 + ℓ + 1
        α[i₀] = a[i₀]
        for m ∈ 1:ℓ
            αₘ = complex(a[i₀ + m], a[i₀ - m])
            α[i₀ + m] = αₘ
            α[i₀ - m] = (-1)^m * conj(αₘ)
        end
    end
    return v⃗, Rotor(exp(r⃗ / 2)), α
end

# Resample the modes `d` (one row for each mode, sampled at the times `tk`) onto the times
# `tq`, the way `transform!` itself interpolates on `main`: one natural cubic spline for
# each mode, the Thomas sweeps vectorized across the modes, the query times marched through
# in order, and the modes divided among the threads.  `d̈` is a work array like `d`.
function resample!(out, d̈, d, tk, tq)
    Nᵐ, Nᵏ = size(d)
    c = Scri.CubicSplineCache(tk)
    h, h⁻¹, l, u⁻¹ = c.h, c.h⁻¹, c.l, c.u⁻¹
    nchunks = Threads.nthreads()
    Threads.@threads :static for chunk ∈ 1:nchunks
        ms = (1 + ((chunk - 1) * Nᵐ) ÷ nchunks):((chunk * Nᵐ) ÷ nchunks)
        @inbounds begin
            @simd ivdep for m ∈ ms
                d̈[m, 1] = 0
                d̈[m, Nᵏ] = 0
            end
            @simd ivdep for m ∈ ms
                r = 6 * (h⁻¹[2] * (d[m, 3] - d[m, 2]) - h⁻¹[1] * (d[m, 2] - d[m, 1]))
                d̈[m, 2] = r * u⁻¹[1]
            end
            for j ∈ 3:(Nᵏ - 1)
                @simd ivdep for m ∈ ms
                    r =
                        6 * (
                            h⁻¹[j] * (d[m, j + 1] - d[m, j]) -
                            h⁻¹[j - 1] * (d[m, j] - d[m, j - 1])
                        )
                    d̈[m, j] = (r - h[j - 1] * d̈[m, j - 1]) * u⁻¹[j - 1]
                end
            end
            q = length(tq)
            for j ∈ (Nᵏ - 1):-1:1
                if 2 ≤ j ≤ Nᵏ - 2
                    @simd ivdep for m ∈ ms
                        d̈[m, j] -= l[j - 1] * d̈[m, j + 1]
                    end
                end
                while q ≥ 1 && (tq[q] ≥ tk[j] || j == 1)
                    τ = tq[q] - tk[j]
                    @simd ivdep for m ∈ ms
                        out[m, q] = Scri.spline_eval(
                            d[m, j], d[m, j + 1], d̈[m, j], d̈[m, j + 1], h[j], h⁻¹[j], τ
                        )
                    end
                    q -= 1
                end
            end
        end
    end
    return out
end

# The same resampling written plainly: one mode at a time on one thread, a tridiagonal solve
# from LinearAlgebra, and a binary search for every query time.
function resample_plain!(out, d, tk, tq)
    Nᵐ, Nᵏ = size(d)
    h = diff(tk)
    F = lu(Tridiagonal(h[2:(end - 1)], 2 .* (h[1:(end - 1)] .+ h[2:end]), h[2:(end - 1)]))
    for m ∈ 1:Nᵐ
        y = d[m, :]
        Δ = diff(y) ./ h
        r = 6 .* (Δ[2:end] .- Δ[1:(end - 1)])
        ÿ = [0; complex.(F \ real.(r), F \ imag.(r)); 0]
        for (q, τq) ∈ enumerate(tq)
            j = clamp(searchsortedlast(tk, τq), 1, Nᵏ - 1)
            out[m, q] = Scri.spline_eval(
                y[j], y[j + 1], ÿ[j], ÿ[j + 1], h[j], inv(h[j]), τq - tk[j]
            )
        end
    end
    return out
end

function norm²(a, b)
    s = 0.0
    @inbounds @simd for i ∈ eachindex(a, b)
        s += abs2(a[i] - b[i])
    end
    return s
end

function forward_differences!(f::F, ∇, θ; h=1e-6) where {F}
    f₀ = f(θ)
    θ′ = copy(θ)
    for p ∈ eachindex(θ)
        θ′[p] = θ[p] + h
        ∇[p] = (f(θ′) - f₀) / h
        θ′[p] = θ[p]
    end
    return f₀
end
function central_differences!(f::F, ∇, θ; h=1e-5) where {F}
    f₀ = f(θ)
    θ′ = copy(θ)
    for p ∈ eachindex(θ)
        θ′[p] = θ[p] + h
        f₊ = f(θ′)
        θ′[p] = θ[p] - h
        f₋ = f(θ′)
        ∇[p] = (f₊ - f₋) / 2h
        θ′[p] = θ[p]
    end
    return f₀
end

# Time one evaluation of the objective `f`, and then its gradient by forward and by central
# differences.  A gradient that would take longer than `limit` seconds is not run; its time
# is estimated as the number of evaluations times the time of one.
function time_variant(io, ℓₘₐₓ, components, name, f::F, θ; budget, limit) where {F}
    P = length(θ)
    ∇ = similar(θ)
    t₁, n = least_time(() -> f(θ); budget)
    record(io, ℓₘₐₓ, components, P, "$name: one evaluation", t₁, n, "measured")
    for (kind, differences!, N) ∈ (
        ("forward", forward_differences!, P + 1), ("central", central_differences!, 2P + 1)
    )
        if N * t₁ > limit
            record(io, ℓₘₐₓ, components, P, "$name: $kind", N * t₁, 0, "estimated")
        else
            tₙ, n = least_time(() -> differences!(f, ∇, θ); budget)
            record(io, ℓₘₐₓ, components, P, "$name: $kind", tₙ, n, "measured")
        end
    end
    return nothing
end

function time_main(io, ℓₘₐₓ, Nᵗ, components, ℓα; budget, limit)
    (; t, data, target, θ) = problem(ℓₘₐₓ, Nᵗ, components, ℓα)
    Nᵈ = length(components)
    dc = Scri.DataComponents(components...)
    work, out, d̈ = similar(data), similar(data), similar(data)
    # The target's times are fixed, and inside the window of output times for every
    # transformation near this one.
    copyto!(work, data)
    _, t′₀ = transform!(work, t, bms(θ)..., dc)
    tq = collect(LinRange(t′₀[begin] + 1, t′₀[end] - 1, Nᵗ))
    function objective_without_resampling(θ)
        copyto!(work, data)
        d′, _ = transform!(work, t, bms(θ)..., dc)
        return norm²(d′, target)
    end
    function objective_careful(θ)
        copyto!(work, data)
        d′, t′ = transform!(work, t, bms(θ)..., dc)
        for k ∈ 1:Nᵈ
            @views resample!(out[:, :, k], d̈[:, :, k], d′[:, :, k], t′, tq)
        end
        return norm²(out, target)
    end
    function objective_plain(θ)
        copyto!(work, data)
        d′, t′ = transform!(work, t, bms(θ)..., dc)
        for k ∈ 1:Nᵈ
            @views resample_plain!(out[:, :, k], d′[:, :, k], t′, tq)
        end
        return norm²(out, target)
    end
    # The two resamplings must agree, or the comparison means nothing.
    careful, plain = objective_careful(θ), objective_plain(θ)
    isapprox(careful, plain; rtol=1e-10) ||
        error("the two resamplings disagree: $careful and $plain")
    kw = (; budget, limit)
    time_variant(
        io, ℓₘₐₓ, components, "no resampling", objective_without_resampling, θ; kw...
    )
    time_variant(io, ℓₘₐₓ, components, "careful resampling", objective_careful, θ; kw...)
    time_variant(io, ℓₘₐₓ, components, "plain resampling", objective_plain, θ; kw...)
    return nothing
end

function main(args)
    (; file, Nᵗ, ℓs, ℓαs, budget, limit) = parse_options(args)
    println("main, threads = $(Threads.nthreads()), Nᵗ = $Nᵗ")
    open(file, "w") do io
        for ℓₘₐₓ ∈ ℓs, components ∈ component_sets, ℓα ∈ ℓαs
            ℓα ≤ ℓₘₐₓ || continue
            time_main(io, ℓₘₐₓ, Nᵗ, components, ℓα; budget, limit)
        end
    end
end
main(ARGS)
