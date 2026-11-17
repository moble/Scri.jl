# The part of the comparison that runs on the branch that is checked out: one objective,
# the objective with its gradient from `objective_and_gradient!`, and the gradient from
# ForwardDiff through `transform_objective`.  Run by `run.jl`; see there.

using Scri
import ForwardDiff
include("problem.jl")

function time_branch(io, ℓₘₐₓ, Nᵗ, components, ℓα; budget, limit)
    (; t, data, target, θ) = problem(ℓₘₐₓ, Nᵗ, components, ℓα)
    P = length(θ)
    dc = DataComponents(components...)
    plan = ObjectivePlan(target, t, dc, compute_t′(t, 0.1, 2.0))
    J(θ) = transform_objective(data, t, bms_from_parameters(θ)..., dc, plan)
    ∇ = similar(θ)
    tJ, n = least_time(() -> J(θ); budget)
    record(io, ℓₘₐₓ, components, P, "objective", tJ, n, "measured")
    tG, n = least_time(() -> objective_and_gradient!(∇, θ, data, t, dc, plan); budget)
    record(io, ℓₘₐₓ, components, P, "gradient", tG, n, "measured")
    # ForwardDiff takes between one and a little over two objectives for each parameter.
    # It is skipped if even one objective for each would take longer than the limit.
    if P * tJ > limit
        record(io, ℓₘₐₓ, components, P, "ForwardDiff", 1.6 * P * tJ, 0, "skipped")
    else
        tF, n = least_time(() -> ForwardDiff.gradient(J, θ); budget)
        record(io, ℓₘₐₓ, components, P, "ForwardDiff", tF, n, "measured")
    end
    # The objective again, in case the machine's load changed during the measurements.
    tJ′, n = least_time(() -> J(θ); budget)
    record(io, ℓₘₐₓ, components, P, "objective, repeated", tJ′, n, "measured")
    return nothing
end

function main(args)
    (; file, Nᵗ, ℓs, ℓαs, budget, limit) = parse_options(args)
    println("this branch, threads = $(Threads.nthreads()), Nᵗ = $Nᵗ")
    open(file, "w") do io
        for ℓₘₐₓ ∈ ℓs, components ∈ component_sets, ℓα ∈ ℓαs
            ℓα ≤ ℓₘₐₓ || continue
            time_branch(io, ℓₘₐₓ, Nᵗ, components, ℓα; budget, limit)
        end
    end
end
main(ARGS)
