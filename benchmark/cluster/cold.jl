# A cold start: a new process that loads Scri, transforms the data once at full size, and
# exits, as a job on a cluster does.  Run by `run.jl`, which times the whole process from
# launch to exit; see there.
#
#     julia --project=<environment> cold.jl <ℓₘₐₓ> <components> <ℓα> <Nᵗ> <file>
#
# The results file has lines of the kind that `benchmark/headline/case.jl` writes: the
# seconds to load the package, the seconds of the transformation, which include compiling
# it, and the peak memory.  The transformation has a fixed boost and rotation, like those
# of `case.jl`, and a supertranslation with random modes up to `ℓα`.  Nothing else is
# loaded, so that the process does only what a user's job would do.

const process_start = time_ns()
using Scri
const load_seconds = (time_ns() - process_start) / 1e9

# Random modes up to `ℓα` of a real supertranslation, laid out as in `case.jl`.
function supertranslation(ℓα)
    a = 0.02 .* randn((ℓα + 1)^2)
    α = zeros(ComplexF64, (ℓα + 1)^2)
    for ℓ ∈ 0:ℓα
        i₀ = ℓ^2 + ℓ + 1
        α[i₀] = a[i₀]
        for m ∈ 1:ℓ
            αₘ = complex(a[i₀ + m], a[i₀ - m])
            α[i₀ + m] = αₘ
            α[i₀ - m] = (-1)^m * conj(αₘ)
        end
    end
    return α
end

function main(args)
    ℓₘₐₓ, ℓα, Nᵗ = parse(Int, args[1]), parse(Int, args[3]), parse(Int, args[4])
    components = Symbol.(split(args[2], ','))
    dc = Scri.DataComponents(components...)
    data = randn(ComplexF64, (ℓₘₐₓ + 1)^2, Nᵗ, length(components))
    t = collect(LinRange(0.0, (Nᵗ - 1) / 2, Nᵗ))
    v⃗ = Scri.QuatVec(0.02, 0.01, -0.03)
    R = Scri.Rotor(exp(Scri.QuatVec(0.1, -0.2, 0.3) / 2))
    α = supertranslation(ℓα)
    seconds = @elapsed transform!(data, t, v⃗, R, α, dc)
    open(args[5], "w") do io
        println(io, "info\tload\t", load_seconds)
        println(io, "info\tfirst call at full size\t", seconds)
        return println(io, "info\tmemory after one call\t", Sys.maxrss())
    end
    return nothing
end
main(ARGS)
