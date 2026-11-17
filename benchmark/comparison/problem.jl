# The problems that the comparison times, and the helpers that time them.  This file is
# included by `branch.jl` and `main.jl`, which run on different versions of Scri, so it uses
# only Base and the standard libraries.  Both scripts then time exactly the same data and
# the same transformation.

using Printf: @printf, @sprintf
using Random: Xoshiro, randn

# The sets of components that are compared: the strain alone, the strain with three Weyl
# components, and the strain with all five.
const component_sets = ((:h,), (:h, :ψ₄, :ψ₃, :ψ₂), (:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀))

"""
    problem(ℓₘₐₓ, Nᵗ, components, ℓα)

Random mode data and a random target for the given band limit, number of time samples, and
components, with a parameter vector `θ` that describes a small but generic BMS
transformation whose supertranslation has modes up to `ℓα`.  The time step is 0.5, as for
a typical simulation.  The same arguments always give the same problem.
"""
function problem(ℓₘₐₓ, Nᵗ, components, ℓα)
    rng = Xoshiro(20261116)
    Nᵐ = (ℓₘₐₓ + 1)^2
    Nᵈ = length(components)
    t = collect(LinRange(0.0, (Nᵗ - 1) / 2, Nᵗ))
    data = randn(rng, ComplexF64, Nᵐ, Nᵗ, Nᵈ)
    target = randn(rng, ComplexF64, Nᵐ, Nᵗ, Nᵈ)
    θ = [0.1, -0.2, 0.3, 0.02, 0.01, -0.03, (0.02 .* randn(rng, (ℓα + 1)^2))...]
    return (; t, data, target, θ)
end

"""
    least_time(f; budget=10.0, nmax=50)

Call `f()` once to compile it, and then return the least of several timings of it, in
seconds, with the number of timings taken.  There are at least three timings, and more
until `budget` seconds have passed or there are `nmax` of them.  If the first call alone
takes longer than the budget, there is only one more.
"""
function least_time(f::F; budget=10.0, nmax=50) where {F}
    first_call = @elapsed f()
    nmin = first_call > budget ? 1 : 3
    best, n, start = Inf, 0, time()
    while n < nmin || (n < nmax && time() - start < budget)
        best = min(best, @elapsed f())
        n += 1
    end
    return best, n
end

# One line of results: tab-separated, to be read back by `run.jl`.  `kind` is "measured",
# "estimated", or "skipped".
function record(io, ℓₘₐₓ, components, P, name, seconds, samples, kind)
    println(io, join((ℓₘₐₓ, join(components, ","), P, name, seconds, samples, kind), '\t'))
    flush(io)
    @printf(
        "  ℓₘₐₓ = %2d │ %-18s │ P = %2d │ %-28s │ %10.4f s │ %s\n",
        ℓₘₐₓ,
        join(components, " "),
        P,
        name,
        seconds,
        kind == "measured" ? "least of $samples" : kind
    )
    return flush(stdout)
end

# The options that `run.jl` passes to both scripts: the results file, the number of time
# samples, the band limits, the supertranslations' band limits, the budget for each
# measurement, and the longest that one measurement may be expected to take, all in seconds.
function parse_options(args)
    file, Nᵗ, ℓs, ℓαs, budget, limit = args
    return (;
        file,
        Nᵗ=parse(Int, Nᵗ),
        ℓs=parse.(Int, split(ℓs, ',')),
        ℓαs=parse.(Int, split(ℓαs, ',')),
        budget=parse(Float64, budget),
        limit=parse(Float64, limit),
    )
end
