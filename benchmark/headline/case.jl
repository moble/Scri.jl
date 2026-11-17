# One case of the headline benchmark, run in a process of its own, so that the peak memory
# of the process belongs to this case alone.  Run by `run.jl`; see there.
#
#     julia --project=<environment> -t <threads> case.jl <kind> <ℓₘₐₓ> <components> <ℓα>
#         <ℓ_fit> <Nᵗ> <budget> <limit> <file>
#
# The kinds are
#
#   - `transform`: one call of `transform!`, which works in place;
#   - `transform_plan`: the same with a `TransformPlan`, which only the branch has;
#   - `gradient`: one call of `objective_and_gradient!`, which only the branch has; and
#   - `differences`: the gradient of the same kind of objective by finite differences (the
#     forward kind, with one evaluation for each parameter and one more) through
#     `transform!`, with the transformed modes resampled onto the target's times by a
#     plainly written cubic spline, which is what `main` offers.
#
# Each timed call is measured by Chairmarks, one call to a sample.  The first call at full
# size is made the way a user would make it, and the process's peak memory is read just
# before and just after it.  Everything is compiled beforehand on a short problem, and the
# time of that first, compiling call is recorded too.
#
# Every timed call is inside a function whose arguments have concrete types.

const process_start = time_ns()
using Scri
const load_seconds = (time_ns() - process_start) / 1e9
using Chairmarks: Chairmarks, Sample
using LinearAlgebra: BLAS, Tridiagonal, lu
using Random: Xoshiro, randn
using Statistics: median
const QuatVec, Rotor = Scri.QuatVec, Scri.Rotor

# The number of time samples of the short problem that compiles everything.
const Nᵗ_short = 400

# The parameters of a small but generic BMS transformation whose supertranslation has modes
# up to `ℓα`: a rotation vector, a rapidity vector, and the real degrees of freedom of the
# supertranslation's modes, as for `bms_from_parameters` on the branch.
function parameters(ℓα)
    return [
        0.1, -0.2, 0.3, 0.02, 0.01, -0.03, (0.02 .* randn(Xoshiro(20261116), (ℓα + 1)^2))...
    ]
end

# Time samples with a step of 0.5, as for a typical simulation, and random mode weights.
times(Nᵗ) = collect(LinRange(0.0, (Nᵗ - 1) / 2, Nᵗ))
modes(seed, ℓₘₐₓ, Nᵗ, Nᵈ) = randn(Xoshiro(seed), ComplexF64, (ℓₘₐₓ + 1)^2, Nᵗ, Nᵈ)

# The boost velocity, rotor, and supertranslation modes for the parameters `θ`.  `main` has
# no `bms_from_parameters`, so both versions use this.
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

# One timed call of `f(setup())`, measured by Chairmarks, with the time of `setup` left out.
function sample(setup::S, f::F) where {S,F}
    return only(Chairmarks.benchmark(setup, f; evals=1, samples=0).samples)
end
sample(f::F) where {F} = sample(nothing, f)

# Add timed calls to `samples` until there are at least three of them (two, if a call takes
# longer than 100 seconds), and then until `budget` seconds have gone by or there are 30.
function repeat_samples!(samples, setup::S, f::F; budget) where {S,F}
    start = time()
    while length(samples) < 30
        n = length(samples)
        enough = n ≥ 3 || (n == 2 && minimum(s.time for s ∈ samples) > 100)
        enough && time() - start ≥ budget && break
        push!(samples, sample(setup, f))
    end
    return samples
end

# The lines of the results file, which `run.jl` reads back: tab-separated, the first field
# saying what kind of line it is.
info(io, name, value) = (println(io, join(("info", name, value), '\t')); flush(io))
memory(io, stage) = info(io, "memory $stage", Sys.maxrss())
function record(io, name, samples::Vector{Sample}; how="measured")
    seconds = [s.time for s ∈ samples]
    row = (
        "time",
        name,
        minimum(seconds),
        median(seconds),
        length(seconds),
        median(s.bytes for s ∈ samples),
        median(s.gc_fraction for s ∈ samples),
        how,
    )
    println(io, join(row, '\t'))
    flush(io)
    println(
        "    $name: $(round(minimum(seconds); sigdigits=4)) s, least of $(length(seconds))"
    )
    return flush(stdout)
end
function record_estimate(io, name, seconds)
    println(io, join(("time", name, seconds, seconds, 0, 0, 0, "estimated"), '\t'))
    flush(io)
    println("    $name: about $(round(seconds; sigdigits=4)) s, estimated")
    return flush(stdout)
end

################################################################################
# A plain transformation, with or without a `TransformPlan`.

apply!(data, t, v⃗, R, α, dc, ::Nothing) = transform!(data, t, v⃗, R, α, dc)
apply!(data, t, v⃗, R, α, dc, plan) = transform!(data, t, v⃗, R, α, dc; plan)

function run_transform(io, c, dc; with_plan)
    v⃗, R, α = bms(parameters(c.ℓα))
    short = modes(1, c.ℓₘₐₓ, Nᵗ_short, length(c.components))
    plan = nothing
    first_call = @elapsed begin
        with_plan && (plan = Scri.TransformPlan(short, dc))
        apply!(short, times(Nᵗ_short), v⃗, R, α, dc, plan)
    end
    info(io, "first call", first_call)
    return time_transform(io, c, dc, v⃗, R, α, plan)
end

function time_transform(io, c, dc, v⃗, R, α, plan)
    Nᵈ = length(c.components)
    t = times(c.Nᵗ)
    # The first call transforms the data in place, as a user would, so that the peak memory
    # is that of holding one copy of the data.
    data₁ = modes(1, c.ℓₘₐₓ, c.Nᵗ, Nᵈ)
    GC.gc()
    memory(io, "before")
    samples = [sample(() -> apply!(data₁, t, v⃗, R, α, dc, plan))]
    memory(io, "after one call")
    # The later calls each transform a fresh copy, which is made outside the timing.
    data = modes(1, c.ℓₘₐₓ, c.Nᵗ, Nᵈ)
    work = data₁
    repeat_samples!(
        samples, () -> copyto!(work, data), w -> apply!(w, t, v⃗, R, α, dc, plan); c.budget
    )
    record(io, "transform", samples)
    if !isnothing(plan)
        build = repeat_samples!(
            Sample[],
            nothing,
            () -> Scri.TransformPlan(c.ℓₘₐₓ, dc);
            budget=min(c.budget, 5.0),
        )
        record(io, "TransformPlan", build)
    end
    return nothing
end

################################################################################
# The objective with its gradient, on the branch.

function run_gradient(io, c, dc)
    θ = parameters(c.ℓα)
    Nᵈ = length(c.components)
    first_call = @elapsed begin
        t = times(Nᵗ_short)
        plan = ObjectivePlan(
            modes(2, c.ℓₘₐₓ, Nᵗ_short, Nᵈ), t, dc, compute_t′(t, 0.1, 2.0); c.ℓ_fit
        )
        data = modes(1, c.ℓₘₐₓ, Nᵗ_short, Nᵈ)
        objective_and_gradient!(similar(θ), θ, data, t, dc, plan)
        transform_objective(data, t, bms_from_parameters(θ)..., dc, plan)
    end
    info(io, "first call", first_call)
    return time_gradient(io, c, dc, θ)
end

function time_gradient(io, c, dc, θ)
    Nᵈ = length(c.components)
    t = times(c.Nᵗ)
    t′ = compute_t′(t, 0.1, 2.0)
    data = modes(1, c.ℓₘₐₓ, c.Nᵗ, Nᵈ)
    target = modes(2, c.ℓₘₐₓ, c.Nᵗ, Nᵈ)
    GC.gc()
    memory(io, "before the plan")
    plan_seconds = @elapsed plan = ObjectivePlan(target, t, dc, t′; c.ℓ_fit)
    ∇ = similar(θ)
    memory(io, "before")
    samples = [sample(() -> objective_and_gradient!(∇, θ, data, t, dc, plan))]
    memory(io, "after one call")
    repeat_samples!(
        samples, nothing, () -> objective_and_gradient!(∇, θ, data, t, dc, plan); c.budget
    )
    record(io, "gradient", samples)
    objective = repeat_samples!(
        Sample[],
        nothing,
        () -> transform_objective(data, t, bms_from_parameters(θ)..., dc, plan);
        c.budget,
    )
    record(io, "objective", objective)
    # The plan is built twice more, now that the peak memory has been read.
    build = [Sample(; time=plan_seconds, warmup=false)]
    for _ ∈ 1:2
        push!(build, sample(() -> ObjectivePlan(target, t, dc, t′; c.ℓ_fit)))
    end
    record(io, "ObjectivePlan", build)
    return nothing
end

################################################################################
# The gradient by finite differences through `transform!`.
#
# `transform!` on `main` chooses its own output times, which move with the transformation.
# So an objective that compares the transformed waveform with a fixed target has to resample
# the transformed modes onto the target's times at every evaluation.  Only the modes up to
# `ℓ_fit` are compared, so only they are resampled.

# A natural cubic spline through the samples `y` at the times `tk`, evaluated at the times
# `tq`, for one mode at a time: a tridiagonal solve from LinearAlgebra and a binary search
# for every query time, as one might write it in ten minutes.
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
            a = (tk[j + 1] - τq) / h[j]
            b = (τq - tk[j]) / h[j]
            out[m, q] =
                a * y[j] +
                b * y[j + 1] +
                ((a^3 - a) * ÿ[j] + (b^3 - b) * ÿ[j + 1]) * h[j]^2 / 6
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

# One evaluation of the objective: transform a copy of the data, resample the modes up to
# `ℓ_fit` onto the target's times, and sum the squares of the differences.
function evaluate!(work, out, data, target, t, tq, dc, θ)
    copyto!(work, data)
    d′, t′ = transform!(work, t, bms(θ)..., dc)
    Nᶠ = size(target, 1)
    for k ∈ axes(target, 3)
        @views resample_plain!(out[:, :, k], d′[1:Nᶠ, :, k], t′, tq)
    end
    return norm²(out, target)
end

function run_differences(io, c, dc)
    θ = parameters(c.ℓα)
    first_call = @elapsed evaluation(c, dc, θ, Nᵗ_short)(θ)
    info(io, "first call", first_call)
    return time_differences(io, c, θ, evaluation(c, dc, θ, c.Nᵗ))
end

# The objective for `Nᵗ` time samples, as a function of the parameters.  The target's times
# are fixed, and inside the window of output times for every transformation near this one.
function evaluation(c, dc, θ, Nᵗ)
    Nᵈ = length(c.components)
    t = times(Nᵗ)
    target = modes(2, c.ℓ_fit, Nᵗ, Nᵈ)
    data = modes(1, c.ℓₘₐₓ, Nᵗ, Nᵈ)
    work, out = similar(data), similar(target)
    copyto!(work, data)
    _, t′₀ = transform!(work, t, bms(θ)..., dc)
    tq = collect(LinRange(t′₀[begin] + 1, t′₀[end] - 1, Nᵗ))
    return θ -> evaluate!(work, out, data, target, t, tq, dc, θ)
end

function time_differences(io, c, θ, f::F) where {F}
    P = length(θ)
    ∇ = similar(θ)
    GC.gc()
    memory(io, "before")
    samples = [sample(() -> f(θ))]
    memory(io, "after one call")
    repeat_samples!(samples, nothing, () -> f(θ); c.budget)
    record(io, "evaluation", samples)
    predicted = (P + 1) * minimum(s.time for s ∈ samples)
    if predicted > c.limit
        record_estimate(io, "differences", predicted)
    else
        samples = [sample(() -> forward_differences!(f, ∇, θ))]
        memory(io, "after one gradient")
        repeat_samples!(samples, nothing, () -> forward_differences!(f, ∇, θ); c.budget)
        record(io, "differences", samples)
    end
    return nothing
end

################################################################################

function main(args)
    kind, ℓₘₐₓ, components, ℓα, ℓ_fit, Nᵗ, budget, limit, file = args
    c = (;
        ℓₘₐₓ=parse(Int, ℓₘₐₓ),
        components=Tuple(Symbol.(split(components, ','))),
        ℓα=parse(Int, ℓα),
        ℓ_fit=parse(Int, ℓ_fit),
        Nᵗ=parse(Int, Nᵗ),
        budget=parse(Float64, budget),
        limit=parse(Float64, limit),
    )
    dc = Scri.DataComponents(c.components...)
    open(file, "w") do io
        info(io, "load average", Sys.loadavg()[1])
        info(io, "threads", Threads.nthreads())
        info(io, "BLAS threads", BLAS.get_num_threads())
        info(io, "load", load_seconds)
        if kind == "transform"
            run_transform(io, c, dc; with_plan=false)
        elseif kind == "transform_plan"
            run_transform(io, c, dc; with_plan=true)
        elseif kind == "gradient"
            run_gradient(io, c, dc)
        elseif kind == "differences"
            run_differences(io, c, dc)
        else
            error("unknown kind `$kind`")
        end
        return memory(io, "at the end")
    end
    return nothing
end
main(ARGS)
