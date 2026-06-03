#!/bin/bash -l
# -*- mode: julia -*-
#SBATCH --output=%x_%j.log

#=
# Run as an executable, the shebang at the top of this file starts this as a bash script.

# The `SBATCH` directives allow us to tell `slurm`` how to run this job.

# This figures out where the *original* script was located, even though `slurm` copies it
# to a different location when running.  This is important because the script needs to know
# where its `Project.toml` is, and where to save output files.
export THIS_COMMAND=$(scontrol show job $SLURM_JOBID | grep "^   Command=" | head -n 1 | cut -d "=" -f 2-)

# The following line replaces the `bash` command with this file as a julia script.  All
# the comments at the top (including this whole block) are ignored by `julia`.
exec julia ${JULIA_SYSIMAGE:-} --threads=auto "${BASH_SOURCE[0]}" "$@"

=#

"""
# Spline errors

The purpose of this script is to demonstrate the fundamental errors that arise from using
splines to interpolate waveforms during the BMS transformation.  We generate a simple test
waveform with a known analytic form, apply a boost to it (no rotation or supertranslation),
and look at the "energy" in the transformed modes as functions of time.  Near zero, we see
clean exponential convergence to zero with increasing ℓ, but as we approach the critical
time samples at ±1/2β, we see a plateau of large errors that does not significantly improve
with increasing ℓ.  This is a fundamental feature of using interpolants with finite
smoothness, even when using extremely high-precision numerics.

To run this script, change to this directory, set up the project by running

    julia --project=.. -e 'using Pkg; Pkg.instantiate()'

and then run something like

    julia --threads=auto spline_errors.jl Float64 16

The raw data and plot will be saved to this directory.

To submit to `slurm`, run something like

    sbatch \
      --ntasks=1 \
      --cpus-per-task=192 \
      --mem-per-cpu=2G \
      -t 36:00:00 \
      spline_errors.jl BigFloat 35

"""

# We have to do this to deal with slurm's trickery involving copying the script
# to a different location upon submission, and running from there.
if "THIS_COMMAND" ∈ keys(ENV)
    THIS_COMMAND = ENV["THIS_COMMAND"]
    THIS_FILE = split(lstrip(THIS_COMMAND), " ")[1]
else
    THIS_COMMAND = PROGRAM_FILE * " " * join(ARGS, " ")
    THIS_FILE = abspath(@__FILE__)
end
THIS_DIR = dirname(THIS_FILE)

using Pkg
Pkg.activate(THIS_DIR)

using Random
using Quaternionic
using Scri
using Plots
using LaTeXStrings
using DoubleFloats
using JLD2
import ArgParse: ArgParse

ENV["GKSwstype"] = "100"  # suppress GR GUI window; still writes PNG files

function generate_parameters(T, ℓₘₐₓ; Nᵗ, t₁, t₂, β)
    rng = Random.Xoshiro(123)

    Nᵐ = (ℓₘₐₓ+1)^2
    ℓₘₐₓ₀ = min(8, ℓₘₐₓ)
    Nᵐ₀ = (ℓₘₐₓ₀ + 1)^2
    data_components = (:σ,)
    Nᵈ = length(data_components)
    dc = Scri.DataComponents(data_components...)

    t = collect(LinRange{BigFloat}(t₁, t₂, Nᵗ));

    δt = t[2] - t[1]
    Ωₙ = π / δt  # Nyquist frequency
    Ω = Ωₙ / 8ℓₘₐₓ₀
    v = ∛Ω
    v⃗ = β * normalize(randn(rng, QuatVec{BigFloat}))
    R = one(Rotor{T})
    αᵢₙ = zeros(Complex{T}, Nᵐ)

    # (ℓ, m) for each mode index
    ℓᵢ = [isqrt(i-1) for i ∈ 1:Nᵐ];
    mᵢ = [(i-1) - ℓᵢ[i]*(ℓᵢ[i]+1) for i ∈ 1:Nᵐ];

    c = zeros(Complex{T}, Nᵐ, Nᵈ)
    for (d, comp) ∈ enumerate(data_components)
        s = Scri.spin_weight(Val(comp))
        c[:, d] .= v .^ ℓᵢ
        c[ℓᵢ .< abs(s), d] .= 0
        c[ℓᵢ .> ℓₘₐₓ₀, d] .= 0
    end

    phases = @. exp(-im * mᵢ * Ω * t')  # Nᵐ × Nᵗ; note the *transpose* of t
    data = zeros(Complex{T}, Nᵐ, Nᵗ, Nᵈ)
    for d ∈ 1:Nᵈ
        @. data[:, :, d] = c[:, d] * phases
    end

    M₄ = maximum(ℓ -> v^ℓ * (ℓ * Ω)^4, 2:ℓₘₐₓ)
    ϵₜ = M₄ * δt^4 * ℓₘₐₓ^BigFloat(-3//2)

    (; data, t=T.(t), v⃗=QuatVec{T}(v⃗), R, αᵢₙ, M₄=T(M₄), ϵₜ=T(ϵₜ))
end

function compute_example(T, ℓₘₐₓ; Nᵗ=10_001, t₁=-5_000, t₂=5_000, β::Rational{Int}=1//1000)
    parameters = generate_parameters(T, ℓₘₐₓ; Nᵗ, t₁, t₂, β)
    data, t, v⃗, R, αᵢₙ, M₄, ϵₜ = parameters
    dc = Scri.DataComponents(:σ)

    N = floor(Int, 1/2β)  # Time samples away from 0 where plateau begins
    t₋, t₊ = let
        i₀ = argmin(abs.(t))
        t[i₀-N], t[i₀+N]
    end

    t1 = time_ns()
    data′, t′ = Scri.transform!(copy(data), t, v⃗, R, αᵢₙ; data_components=(:σ,))
    t2 = time_ns()
    println("Transformation took $((t2-t1)/1e9) seconds")
    JLD2.@save "spline_errors_$T.jld2" data t v⃗ data′ t′ M₄ ϵₜ t₋ t₊
    (; data, t, v⃗, data′, t′, M₄, ϵₜ, t₋, t₊)
end

function plot_example(data::Array{Complex{T}}, t, v⃗, data′, t′, M₄, ϵₜ, t₋, t₊, β) where {T}
    ℓₘₐₓ = isqrt(size(data, 1)) - 1
    beta = Float64(round(β, sigdigits=3))
    dc = Scri.DataComponents(:σ)
    plt = Scri.diagnostics(t′, data′, dc)[:σ]
    hline!(plt, [ℓₘₐₓ * eps(T)], color=:grey, ls=:dot, z_order=:back, label="")
    hline!(plt, [ϵₜ], color=:purple, z_order=:back, label=L"|\sigma^{(4)}| δt^4 ℓ_\mathrm{max}^{-3/2}")
    vline!(plt, [t₋, t₊], color=:black, z_order=:back, label=L"N = \lfloor 1/(2β) \rfloor")
    plot!(plt, legend_position=(0.82,0.99), title=L"\sigma' \quad(\beta=%$beta c)", dpi=300)
    savefig(plt, "spline_errors_$T.png")
end

function main(args=ARGS)
    println("Running on $(Threads.nthreads()) threads")

    s = ArgParse.ArgParseSettings(
        description="Demonstrate spline interpolation errors in BMS transformations"
    )
    ArgParse.@add_arg_table! s begin
        "T"
            help = "Floating-point precision type (e.g. Float64, Double64, BigFloat)"
            required = true
        "lmax"
            help = "Maximum harmonic mode ℓₘₐₓ"
            arg_type = Int
            required = true
        "--Nt"
            help = "Number of time samples"
            arg_type = Int
            default = 10_001
        "--t1"
            help = "Start time"
            arg_type = Float64
            default = -5_000.0
        "--t2"
            help = "End time"
            arg_type = Float64
            default = 5_000.0
        "--beta"
            help = "Boost velocity as rational p/q (e.g., 1/1000)"
            default = "1/1000"
    end

    parsed = ArgParse.parse_args(args, s)

    T = eval(Meta.parse(parsed["T"]))
    ℓₘₐₓ = parsed["lmax"]
    Nᵗ = parsed["Nt"]
    t₁ = parsed["t1"]
    t₂ = parsed["t2"]
    β = let parts = split(parsed["beta"], '/')
        length(parts) == 2 ? Rational(parse(Int, parts[1]), parse(Int, parts[2])) : Rational(parse(Int, parts[1]))
    end

    result = compute_example(T, ℓₘₐₓ; Nᵗ, t₁, t₂, β)
    plot_example(result..., β)
    return 0
end

@isdefined(var"@main") ? (@main) : exit(main(ARGS))
