# A timing comparison of the ways to get the gradient of the objective with respect to the
# parameters of a BMS transformation.  Run it from anywhere, on a quiet machine:
#
#     julia benchmark/comparison/run.jl
#
# It compares, for the same data and the same transformation,
#
#   - `objective_and_gradient!` on the branch that is checked out;
#   - `ForwardDiff.gradient` through `transform_objective` on that branch; and
#   - finite differences through `transform!` on `main`, forward and central, for three
#     versions of the objective that `main.jl` describes; and
#   - if `--python` is given, finite differences through the transformation of the original
#     Python `scri` package, as `python_scri.py` describes.
#
# Every problem has 20,000 time samples.  The components are the strain alone; the strain
# with ψ₄, ψ₃, and ψ₂; and the strain with all five Weyl components.  The band limits and
# the numbers of parameters are options, with defaults below.
#
# The script changes nothing in the repository except to add its report, as
# `benchmark/comparison_<date>.txt`.  It reads `main` with `git archive`, into a temporary
# directory, and runs the two halves one after the other in processes of their own:
# `branch.jl` in the test environment of this checkout, and `main.jl` in a temporary
# environment that holds the copy of `main`.  Nothing runs concurrently.
#
# Options, each as `--name=value`:
#
#     --threads=8       threads for each half (Julia's; BLAS keeps its own default)
#     --nt=20000        time samples
#     --lmax=8,16       band limits of the data
#     --lalpha=0,4,8    band limits of the supertranslation, so that there are 7, 31, and
#                       87 parameters; those above a data band limit are skipped for it
#     --budget=10       seconds to spend repeating each measurement; each is the least of
#                       at least three timings, unless one timing alone exceeds the budget
#     --limit=900       the longest, in seconds, that one finite-difference gradient may be
#                       expected to take; a longer one is estimated from a single evaluation
#                       instead, and ForwardDiff is skipped if it would take longer
#     --difflimit=900   the same limit for the finite differences on `main` alone, which
#                       defaults to `--limit`; a small value, such as 60, measures only the
#                       cheap gradients, which shows how good the estimates are, and
#                       shortens the run from hours to well under one
#     --python=PATH     the Python executable of an environment that holds the original
#                       `scri` package (`python_scri.py` says how to make one); without
#                       this option, the Python part is left out
#
# `--quick` runs a small version in about a minute, to check that the script works, and
# prints its report without saving it.
#
# The default run takes a few hours, most of it in the finite differences on `main` with
# the plain resampling.  Raising `--limit` measures more of them and takes longer.  The
# Python part never runs a gradient: it times one evaluation, and the report multiplies.

using Dates: Dates
using Printf: @sprintf

const scri_uuid = "ccc3add2-6548-56f5-a99b-c07b0a26ef59"

function parse_arguments(args)
    options = Dict(
        "threads" => "8",
        "nt" => "20000",
        "lmax" => "8,16",
        "lalpha" => "0,4,8",
        "budget" => "10",
        "limit" => "900",
        "python" => "",
        "difflimit" => "",
    )
    quick = "--quick" ∈ args
    if quick
        merge!(
            options,
            Dict(
                "nt" => "400",
                "lmax" => "8",
                "lalpha" => "0,2",
                "budget" => "1",
                "limit" => "5",
            ),
        )
    end
    for a ∈ args
        a == "--quick" && continue
        m = match(r"^--(\w+)=(.+)$", a)
        (isnothing(m) || !haskey(options, m[1])) && error("unknown option `$a`")
        options[m[1]] = m[2]
    end
    return options, quick
end

# The versions of the packages that matter, in the environment `env`.
function versions(julia, env)
    code = """
    using Pkg
    for p in sort(collect(values(Pkg.dependencies())); by=p -> p.name)
        p.name in ("ForwardDiff", "OhMyThreads", "Quaternionic", "SphericalFunctions") &&
            print(p.name, " ", p.version, "; ")
    end
    """
    return strip(
        read(`$julia --startup-file=no --project=$env -e $code`, String), [' ', ';']
    )
end

# Read the tab-separated results of one half into a dictionary from the problem to its
# measurements, in the order they were made.
function read_results(file)
    problems = Tuple{Int,String,Int}[]
    results = Dict{Tuple{Int,String,Int},Vector{NamedTuple}}()
    for line ∈ eachline(file)
        startswith(line, '#') && continue
        ℓₘₐₓ, components, P, name, seconds, _, kind = split(line, '\t')
        key = (parse(Int, ℓₘₐₓ), String(components), parse(Int, P))
        key ∈ problems || push!(problems, key)
        push!(
            get!(results, key, NamedTuple[]),
            (; name=String(name), seconds=parse(Float64, seconds), kind=String(kind)),
        )
    end
    return problems, results
end

function report(io, options, header, branch_file, main_file, python_file)
    println(io, header)
    problems, branch = read_results(branch_file)
    _, main = read_results(main_file)
    python = isnothing(python_file) ? nothing : read_results(python_file)[2]
    seconds(rows, name) = only(r for r ∈ rows if r.name == name)
    mark(r) =
        if r.kind == "measured"
            " "
        elseif r.kind == "estimated"
            "~"
        else
            "?"
        end
    for key ∈ problems
        ℓₘₐₓ, components, P = key
        b, m = branch[key], main[key]
        J = min(seconds(b, "objective").seconds, seconds(b, "objective, repeated").seconds)
        G = seconds(b, "gradient").seconds
        F = seconds(b, "ForwardDiff")
        println(io)
        println(io, "ℓₘₐₓ = $ℓₘₐₓ │ $(replace(components, "," => " ")) │ P = $P")
        println(io, "  this branch", " "^41, "seconds   ÷ objective   ÷ gradient")
        println(io, @sprintf("    %-46s %10.4f %12.1f", "one objective", J, 1.0))
        println(
            io,
            @sprintf(
                "    %-46s %10.4f %12.1f %12.1f", "objective_and_gradient!", G, G / J, 1.0
            ),
        )
        if F.kind == "skipped"
            println(
                io,
                @sprintf(
                    "    %-46s    skipped: it would take about %.0f s",
                    "ForwardDiff",
                    F.seconds
                ),
            )
        else
            println(
                io,
                @sprintf(
                    "    %-46s %10.4f %12.1f %12.1f",
                    "ForwardDiff through transform_objective",
                    F.seconds,
                    F.seconds / J,
                    F.seconds / G
                ),
            )
        end
        println(io, "  main: finite differences through transform!")
        for variant ∈ ("no resampling", "careful resampling", "plain resampling")
            one = seconds(m, "$variant: one evaluation")
            println(
                io, @sprintf("    %-46s %10.4f", "$variant: one evaluation", one.seconds)
            )
            for (kind, N) ∈ (("forward", P + 1), ("central", 2P + 1))
                r = seconds(m, "$variant: $kind")
                println(
                    io,
                    @sprintf(
                        "    %-46s%s%10.4f %25.1f",
                        "$variant: $kind differences ($N)",
                        mark(r),
                        r.seconds,
                        r.seconds / G
                    ),
                )
            end
        end
        if !isnothing(python)
            println(io, "  Python scri: finite differences through its transformation")
            rows = python[key]
            one = seconds(rows, "python: one evaluation")
            println(
                io, @sprintf("    %-46s%s%10.4f", "one evaluation", mark(one), one.seconds)
            )
            for (kind, N) ∈ (("forward", P + 1), ("central", 2P + 1))
                r = seconds(rows, "python: $kind")
                println(
                    io,
                    @sprintf(
                        "    %-46s%s%10.4f %25.1f",
                        "$kind differences ($N)",
                        mark(r),
                        r.seconds,
                        r.seconds / G
                    ),
                )
            end
        end
        # The speedup that a count of evaluations predicts, as a check on the measurements.
        one = seconds(m, "no resampling: one evaluation").seconds
        println(
            io,
            @sprintf(
                "  check: (P + 1) × one evaluation on main ÷ the gradient = %.1f",
                (P + 1) * one / G
            ),
        )
    end
    println(io)
    println(io, "A time marked ~ was not measured.  For a gradient, it is the number of")
    println(io, "evaluations times the time of one, because the gradient would have taken")
    println(io, "longer than the limit of $(options["difflimit"]) s.")
    if isnothing(python)
        println(io, "Not included: the original Python `scri` package.")
    else
        println(
            io, "For Python `scri`, no gradient is ever run, and an evaluation marked ~ was"
        )
        println(
            io, "timed with fewer time samples and scaled up.  It transforms the shear and"
        )
        println(
            io, "all five Weyl components for either of the sets that has Weyl components,"
        )
        println(
            io, "so that transformation was timed once, for the set with three of them, and"
        )
        println(io, "the set with all five repeats those times.")
    end
    return nothing
end

function main(args)
    options, quick = parse_arguments(args)
    here = @__DIR__
    root = normpath(joinpath(here, "..", ".."))
    julia = Base.julia_cmd()
    threads = options["threads"]

    # A copy of `main`, and an environment that uses it.  The directory is kept afterwards,
    # because it also holds each part's results as they are measured, and those are worth
    # having if a later part fails.
    work = mktempdir(; cleanup=false)
    println("Working in $work")
    main_copy, env = joinpath(work, "main"), joinpath(work, "env")
    mkpath(main_copy)
    mkpath(env)
    run(pipeline(`git -C $root archive main`, `tar -x -C $main_copy`))
    write(
        joinpath(env, "Project.toml"),
        "[deps]\nScri = \"$scri_uuid\"\n\n[sources]\nScri = {path = \"$main_copy\"}\n",
    )
    run(
        `$julia --startup-file=no --project=$env -e 'using Pkg; Pkg.instantiate(); using Scri'`,
    )
    test_env = joinpath(root, "test")

    load = Sys.loadavg()
    load[1] > 2 && @warn "The machine is not quiet: its load average is $(load[1])."
    shared = [options[k] for k ∈ ("nt", "lmax", "lalpha", "budget", "limit")]
    isempty(options["difflimit"]) && (options["difflimit"] = options["limit"])
    shared_main = [options[k] for k ∈ ("nt", "lmax", "lalpha", "budget", "difflimit")]
    branch_file, main_file = joinpath(work, "branch.tsv"), joinpath(work, "main.tsv")
    run(
        `$julia --startup-file=no --project=$test_env -t $threads $(joinpath(here, "branch.jl")) $branch_file $shared`,
    )
    run(
        `$julia --startup-file=no --project=$env -t $threads $(joinpath(here, "main.jl")) $main_file $shared_main`,
    )

    # The original Python package, if an environment for it was given.
    python_file = nothing
    python_versions = ""
    if !isempty(options["python"])
        file = joinpath(work, "python.tsv")
        script = joinpath(here, "python_scri.py")
        try
            run(`$(options["python"]) $script $file $shared`)
            python_file = file
            python_versions = "\n# Python scri: " * lstrip(readline(file), '#')
        catch e
            @warn "The Python part failed, and is left out of the report." exception = e
        end
    end

    git(args...) = strip(read(`git -C $root $args`, String))
    header = """
    # A comparison of ways to get the gradient of the objective ($(Dates.today()))
    # $(Sys.cpu_info()[1].model), $threads threads, Julia $VERSION, $(options["nt"]) time samples
    # this branch: $(git("rev-parse", "--abbrev-ref", "HEAD")) at $(git("rev-parse", "--short", "HEAD")); $(versions(julia, test_env))
    # main: at $(git("rev-parse", "--short", "main")); $(versions(julia, env))$python_versions
    # load average at the start: $(@sprintf("%.1f", load[1])); at the end: $(@sprintf("%.1f", Sys.loadavg()[1]))
    # Each time is the least of several, in seconds.  The last column is the time as a
    # multiple of the time for `objective_and_gradient!`, which is the speedup it gives."""
    report(stdout, options, header, branch_file, main_file, python_file)
    if !quick
        file = joinpath(root, "benchmark", "comparison_$(Dates.today()).txt")
        isfile(file) &&
            (file = replace(file, ".txt" => Dates.format(Dates.now(), "_HHMM") * ".txt"))
        open(file, "w") do io
            return report(io, options, header, branch_file, main_file, python_file)
        end
        println("\nSaved as $file")
    end
    return nothing
end
main(ARGS)
