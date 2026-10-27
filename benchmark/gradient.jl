# Benchmark `objective_and_gradient!` against the objective alone and ForwardDiff.
#
# Run from the repository root with threads enabled:
#
#     julia --project=test -t auto benchmark/gradient.jl
#
# Times are best-of-5 wall-clock seconds (best-of-3 for ForwardDiff), after a warm-up.

using Printf: @printf
using Random: Xoshiro, randn
import ForwardDiff
using Scri

function bench(ℓₘₐₓ, Nᵗ, Nα, ℓ_fit)
    rng = Xoshiro(1234)
    N = (ℓₘₐₓ + 1)^2
    t = collect(LinRange(-200.0, 200.0, Nᵗ))
    dc = DataComponents(:h, :ψ₄)
    data = randn(rng, ComplexF64, N, Nᵗ, 2)
    target = randn(rng, ComplexF64, N, Nᵗ, 2)
    plan = ObjectivePlan(target, t, dc, compute_t′(t, 0.1, 1.0); ℓ_fit)
    θ = [0.1, -0.2, 0.3, 0.02, 0.01, -0.03, (0.05 .* randn(rng, Nα))...]
    J(θ) = transform_objective(data, t, bms_from_parameters(θ)..., dc, plan)
    ∇ = similar(θ)
    J(θ), objective_and_gradient!(∇, θ, data, t, dc, plan), ForwardDiff.gradient(J, θ)
    tJ = minimum(@elapsed(J(θ)) for _ ∈ 1:5)
    tG = minimum(@elapsed(objective_and_gradient!(∇, θ, data, t, dc, plan)) for _ ∈ 1:5)
    tF = minimum(@elapsed(ForwardDiff.gradient(J, θ)) for _ ∈ 1:3)
    return length(θ), tJ, tG, tF
end

println("objective_and_gradient! benchmarks, threads = $(Threads.nthreads())")
@printf(
    "%5s │ %5s │ %5s │ %3s │ %10s │ %10s │ %6s │ %10s │ %6s\n",
    "ℓₘₐₓ",
    "ℓ_fit",
    "Nᵗ",
    "P",
    "J [s]",
    "∇J [s]",
    "ratio",
    "FD [s]",
    "ratio"
)
for (ℓₘₐₓ, Nᵗ) ∈ ((8, 400), (16, 400)), ℓ_fit ∈ (nothing, ℓₘₐₓ ÷ 2), Nα ∈ (1, 16)
    P, tJ, tG, tF = bench(ℓₘₐₓ, Nᵗ, Nα, ℓ_fit)
    @printf(
        "%5d │ %5s │ %5d │ %3d │ %10.4f │ %10.4f │ %6.2f │ %10.4f │ %6.2f\n",
        ℓₘₐₓ,
        isnothing(ℓ_fit) ? "—" : string(ℓ_fit),
        Nᵗ,
        P,
        tJ,
        tG,
        tG / tJ,
        tF,
        tF / tJ
    )
end
