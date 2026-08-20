# Benchmark `transform!` across component counts and band limits.
#
# Run from the repository root with threads enabled:
#
#     julia --project=. -t auto benchmark/benchmarks.jl
#
# Optional arguments:  benchmarks.jl [Nᵗ] [ℓ_values...]
# e.g., `julia --project=. -t auto benchmark/benchmarks.jl 2000 8 16 32`.
#
# Times are best-of-3 wall-clock seconds per call, with a warm-up call excluded (the
# warm-up also covers compilation for shapes outside the precompile workload).  Keep the
# machine otherwise idle for stable numbers, and record results before and after any
# performance change.

using Printf: @printf
using Random: Xoshiro
using Quaternionic: QuatVec, Rotor
using Scri

const COMPONENT_SETS = Dict(
    1 => (:h,),
    2 => (:h, :ψ₄),
    4 => (:h, :ψ₄, :ψ₃, :ψ₂),
    6 => (:h, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
    8 => (:σ, :h, :News, :ψ₄, :ψ₃, :ψ₂, :ψ₁, :ψ₀),
)

function bench(ℓₘₐₓ::Int, Nᵈ::Int, Nᵗ::Int; samples::Int=3)
    rng = Xoshiro(1234)
    Nᵐ = (ℓₘₐₓ + 1)^2
    t = collect(LinRange(-500.0, 500.0, Nᵗ))
    v⃗ = QuatVec(1e-4, -2e-4, 3e-4)
    R = Rotor(1.0, 1e-2, -2e-2, 3e-2)
    α = 1e-2 * randn(rng, ComplexF64, 9)
    dc = DataComponents(COMPONENT_SETS[Nᵈ]...)
    mkdata() = randn(rng, ComplexF64, Nᵐ, Nᵗ, Nᵈ)
    transform!(mkdata(), copy(t), v⃗, R, copy(α), dc)  # warm-up / compile
    best_time = Inf
    allocs = 0
    for _ ∈ 1:samples
        data = mkdata()
        stats = @timed transform!(data, copy(t), v⃗, R, copy(α), dc)
        best_time = min(best_time, stats.time)
        allocs = stats.bytes
    end
    return best_time, allocs
end

function main()
    Nᵗ = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 10_000
    ℓs = length(ARGS) ≥ 2 ? parse.(Int, ARGS[2:end]) : [4, 8, 12, 16, 24, 32]
    Nᵈs = sort(collect(keys(COMPONENT_SETS)))
    println("transform! benchmarks: Nᵗ = $Nᵗ, threads = $(Threads.nthreads())")
    @printf("%6s │ %4s │ %12s │ %12s\n", "ℓₘₐₓ", "Nᵈ", "time [s]", "alloc [MiB]")
    println("───────┼──────┼──────────────┼─────────────")
    for ℓₘₐₓ ∈ ℓs, Nᵈ ∈ Nᵈs
        time, bytes = bench(ℓₘₐₓ, Nᵈ, Nᵗ)
        @printf("%6d │ %4d │ %12.4f │ %12.1f\n", ℓₘₐₓ, Nᵈ, time, bytes / 2^20)
    end
end

main()
