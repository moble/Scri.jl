# The headline benchmark: what someone who uses `main`, the branch that is checked out, or
# the original Python `scri` package would experience, in time and in peak memory.  Run it
# from anywhere, on a quiet machine:
#
#     julia benchmark/headline/run.jl --python=PATH
#
# It measures three things, always for waveforms with 20,000 time samples:
#
#   1. a plain transformation, by `transform!` on `main` and on the branch, and by the
#      Python package;
#   2. the same on the branch with a `TransformPlan`; and
#   3. the gradient of an objective that compares the modes up to `ℓ_fit`: by
#      `objective_and_gradient!` on the branch, and by finite differences through the
#      transformation, with a plainly written interpolation onto the target's times, on
#      `main` and in Python.  The finite differences are forward differences, which take
#      one evaluation of the objective for each parameter and one more.
#
# The components are the strain alone, and the strain with all five Weyl components.
#
# Every case runs in a process of its own, one after another, so that each case's peak
# memory is its own and nothing runs concurrently.  `case.jl` and `case.py` describe what
# one case does.  The Julia cases run in two temporary environments, one that holds the
# checkout and one that holds a copy of `main` made with `git archive`, each with
# Chairmarks, which takes the timings.
#
# The script changes nothing in the repository except to add its report, as
# `benchmark/headline_<date>.txt`.  The results of each case are kept in the working
# directory as they are measured, and `--work=DIR` with that directory resumes a run that
# was stopped, repeating only the cases that did not finish.
#
# Options, each as `--name=value`:
#
#     --threads=8,1     Julia's threads; every Julia case is run once for each number
#                       (BLAS keeps its own default)
#     --nt=20000        time samples
#     --lmax=8,12,16    band limits of the data
#     --lalpha=4,8      band limits of the supertranslation, so that there are 31 and 87
#                       parameters
#     --lfit=8          the band limit of the modes that the objective compares
#     --budget=10       seconds to spend repeating each measurement; each is the least of
#                       at least three timings (two, if one timing takes over 100 seconds)
#     --limit=600       the longest, in seconds, that one gradient by finite differences may
#                       be expected to take; a longer one is not run, and its time is the
#                       number of evaluations times the time of one
#     --python=PATH     the Python executable of an environment that holds the original
#                       `scri` package (`case.py` says how to make one); without this
#                       option, the Python cases are left out
#     --work=DIR        the working directory of an earlier run, to resume it
#
# `--quick` runs a small version in a few minutes, to check that the script works, and
# prints its report without saving it.  `--report`, with `--work=DIR`, runs nothing: it
# prints the report of the cases that have finished so far in a run that is under way.

using Dates: Dates
using Printf: @sprintf

const component_sets = (("h",), ("h", "ψ₄", "ψ₃", "ψ₂", "ψ₁", "ψ₀"))

function parse_arguments(args)
    options = Dict(
        "threads" => "8,1",
        "nt" => "20000",
        "lmax" => "8,12,16",
        "lalpha" => "4,8",
        "lfit" => "8",
        "budget" => "10",
        "limit" => "600",
        "python" => "",
        "work" => "",
    )
    quick = "--quick" ∈ args
    report_only = "--report" ∈ args
    if quick
        merge!(
            options,
            Dict(
                "threads" => "4,1",
                "nt" => "400",
                "lmax" => "8,10",
                "lalpha" => "2",
                "budget" => "1",
                "limit" => "5",
            ),
        )
    end
    for a ∈ args
        a ∈ ("--quick", "--report") && continue
        m = match(r"^--(\w+)=(.+)$", a)
        (isnothing(m) || !haskey(options, m[1])) && error("unknown option `$a`")
        options[m[1]] = m[2]
    end
    return options, quick, report_only
end

integers(s) = parse.(Int, split(s, ','))

# The cases, in the order they run: every Julia case with the first number of threads, then
# the Python cases, which take the longest, and then the Julia cases with the other numbers
# of threads.
function all_cases(options)
    problems = [
        (; ℓₘₐₓ, components, ℓα) for ℓₘₐₓ ∈ integers(options["lmax"]) for
        components ∈ component_sets for ℓα ∈ integers(options["lalpha"]) if ℓα ≤ ℓₘₐₓ
    ]
    julia_cases(threads) = [
        (; version, threads, kind, p...) for p ∈ problems for (version, kind) ∈ (
            ("branch", "transform"),
            ("branch", "transform_plan"),
            ("branch", "gradient"),
            ("main", "transform"),
            ("main", "differences"),
        )
    ]
    threads = integers(options["threads"])
    cases = julia_cases(threads[1])
    if !isempty(options["python"])
        for p ∈ problems, kind ∈ ("transform", "differences")
            push!(cases, (; version="python", threads=0, kind, p...))
        end
    end
    for n ∈ threads[2:end]
        append!(cases, julia_cases(n))
    end
    return problems, cases
end

function name(c)
    return "$(c.version)-t$(c.threads)-$(c.kind)-l$(c.ℓₘₐₓ)-c$(length(c.components))-a$(c.ℓα)"
end

# A temporary environment that holds the copy of Scri at `path`, and Chairmarks.
function make_environment(julia, env, path)
    isfile(joinpath(env, "Manifest.toml")) && return nothing
    mkpath(env)
    code = """
    using Pkg
    Pkg.develop(; path=ARGS[1])
    Pkg.add("Chairmarks")
    using Scri, Chairmarks
    """
    run(`$julia --startup-file=no --project=$env -e $code $path`)
    return nothing
end

# The versions of the packages that matter, in the environment `env`.
function versions(julia, env)
    code = """
    using Pkg
    for p in sort(collect(values(Pkg.dependencies())); by=p -> p.name)
        p.name in ("Chairmarks", "OhMyThreads", "Quaternionic", "SphericalFunctions") &&
            print(p.name, " ", p.version, "; ")
    end
    """
    return strip(
        read(`$julia --startup-file=no --project=$env -e $code`, String), [' ', ';']
    )
end

# Run one case in a process of its own, unless its results are already there.
function run_case(c, options, julia, work)
    here = @__DIR__
    file = joinpath(work, "results", name(c) * ".tsv")
    isfile(file) && return nothing
    part = file * ".part"
    args = [
        c.kind,
        string(c.ℓₘₐₓ),
        join(c.components, ','),
        string(c.ℓα),
        options["lfit"],
        options["nt"],
        options["budget"],
        options["limit"],
        part,
    ]
    cmd = if c.version == "python"
        `$(options["python"]) $(joinpath(here, "case.py")) $args`
    else
        env = joinpath(work, "env-" * c.version)
        `$julia --startup-file=no --project=$env -t $(c.threads) $(joinpath(here, "case.jl")) $args`
    end
    println("$(Dates.format(Dates.now(), "HH:MM:SS"))  $(name(c))")
    flush(stdout)
    seconds = @elapsed ok = success(pipeline(cmd; stdout, stderr))
    if ok
        mv(part, file)
    else
        @warn "The case $(name(c)) failed, and is left out of the report."
        isfile(part) && mv(part, file * ".failed"; force=true)
    end
    println("    $(round(seconds; digits=1)) s in all")
    return flush(stdout)
end

# The results of one case: its information lines, and its timings by name.
function read_case(file)
    information = Dict{String,String}()
    timings = Dict{String,NamedTuple}()
    for line ∈ eachline(file)
        fields = split(line, '\t')
        if fields[1] == "info"
            information[fields[2]] = fields[3]
        else
            _, quantity, least, median, n, bytes, gc, how = fields
            timings[quantity] = (;
                seconds=parse(Float64, least),
                median=parse(Float64, median),
                n=parse(Int, n),
                bytes=parse(Float64, bytes),
                gc=parse(Float64, gc),
                how=String(how),
            )
        end
    end
    return (; information, timings)
end

function read_results(cases, work)
    results = Dict{String,NamedTuple}()
    for c ∈ cases
        file = joinpath(work, "results", name(c) * ".tsv")
        isfile(file) && (results[name(c)] = read_case(file))
    end
    return results
end

# Seconds to three significant figures, and bytes as gigabytes.
function digits3(x)
    return if x < 0.9995
        @sprintf("%.3f", x)
    elseif x < 9.995
        @sprintf("%.2f", x)
    elseif x < 99.95
        @sprintf("%.1f", x)
    else
        @sprintf("%.0f", x)
    end
end
ratio(x) = x < 9.95 ? @sprintf("%.1f", x) : @sprintf("%.0f", x)

# A table with the problem in the first columns.  `columns` holds each remaining column's
# heading and a function from the problem to the text of its entry.
function table(io, title, problems, columns)
    headings = ["ℓₘₐₓ"; "components"; "P"; first.(columns)]
    rows = [
        [
            string(p.ℓₘₐₓ)
            (length(p.components) == 1 ? "h" : "h ψ₄ … ψ₀")
            string(6 + (p.ℓα + 1)^2)
            [entry(p) for (_, entry) ∈ columns]
        ] for p ∈ problems
    ]
    widths = [
        maximum(textwidth, [headings[i]; [r[i] for r ∈ rows]]) for i ∈ eachindex(headings)
    ]
    line(cells) = join(
        [
            i == 2 ? rpad(cells[i], widths[i]) : lpad(cells[i], widths[i]) for
            i ∈ eachindex(cells)
        ],
        " │ ",
    )
    println(io)
    println(io, title)
    println(io, "  ", line(headings))
    println(io, "  ", join(["─"^w for w ∈ widths], "─┼─"))
    for r ∈ rows
        println(io, "  ", line(r))
    end
    return nothing
end

function report(io, options, header, problems, cases, work)
    results = read_results(cases, work)
    threads = integers(options["threads"])
    case(version, n, kind, p) =
        get(results, name((; version, threads=n, kind, p...)), nothing)
    # The timing of `quantity` in a case, or nothing if the case is missing.
    function timing(version, n, kind, p, quantity)
        r = case(version, n, kind, p)
        return isnothing(r) ? nothing : get(r.timings, quantity, nothing)
    end
    function seconds(version, n, kind, p, quantity)
        x = timing(version, n, kind, p, quantity)
        return isnothing(x) ? "—" : (x.how == "estimated" ? "~" : "") * digits3(x.seconds)
    end
    function speedup(a, b)
        (isnothing(a) || isnothing(b)) && return "—"
        return (a.how == "estimated" ? "~" : "") * ratio(a.seconds / b.seconds)
    end
    # The peak memory of the process after its first call at full size, in gigabytes.
    function peak(version, n, kind, p)
        r = case(version, n, kind, p)
        isnothing(r) && return "—"
        i = r.information
        bytes = get(i, "memory after one gradient", get(i, "memory after one call", ""))
        return isempty(bytes) ? "—" : @sprintf("%.2f", parse(Float64, bytes) / 2^30)
    end
    function information(version, n, kind, p, key)
        r = case(version, n, kind, p)
        (isnothing(r) || !haskey(r.information, key)) && return "—"
        return digits3(parse(Float64, r.information[key]))
    end
    python = !isempty(options["python"])
    python_column(heading, entry) = python ? [(heading, entry)] : []

    println(io, header)
    println(io)
    println(
        io, "In the headings, `main/8` is `main` with 8 Julia threads, `branch/1` is the"
    )
    println(
        io, "branch with 1 thread, which is what Julia starts with by default, and so on."
    )
    println(
        io, "`plan` is the branch with a `TransformPlan`.  P is the number of parameters."
    )
    println(
        io, "A time marked ~ was not measured: it is the number of evaluations times the"
    )
    println(
        io, "time of one, because the gradient would have taken longer than the limit of"
    )
    println(
        io,
        "$(options["limit"]) s.  Memory is the peak resident set size of the process just after",
    )
    println(
        io, "its first call at full size, in gigabytes (2³⁰ bytes); it includes Julia or"
    )
    println(io, "Python itself, the compiled code, and the data.")
    println(io)
    println(
        io, "The objective's plan is made with `ℓ_fit`, which costs one more pass over the"
    )
    println(
        io,
        "pixels than a plan that compares every mode, even when `ℓ_fit` is the data's band",
    )
    println(
        io, "limit.  So these gradients take longer than those of `benchmark/comparison`."
    )

    println(io)
    println(io, "="^78)
    println(io, "1 and 2.  A plain transformation, with and without a TransformPlan")
    println(io, "="^78)
    T(version, kind) = n -> p -> seconds(version, n, kind, p, "transform")
    table(
        io,
        "Seconds for one call",
        problems,
        [
            python_column("Python", p -> seconds("python", 0, "transform", p, "transform"))
            [
                (heading * "/$n", entry(n)) for n ∈ threads for (heading, entry) ∈ (
                    ("main", T("main", "transform")),
                    ("branch", T("branch", "transform")),
                    ("plan", T("branch", "transform_plan")),
                )
            ]
        ],
    )
    table(
        io,
        "Ratios of those times",
        problems,
        [
            (heading * "/$n", entry(n)) for n ∈ threads for (heading, entry) ∈ (
                (
                    "main ÷ branch",
                    n ->
                        p -> speedup(
                            timing("main", n, "transform", p, "transform"),
                            timing("branch", n, "transform", p, "transform"),
                        ),
                ),
                (
                    "Python ÷ branch",
                    n ->
                        p -> speedup(
                            timing("python", 0, "transform", p, "transform"),
                            timing("branch", n, "transform", p, "transform"),
                        ),
                ),
                (
                    "branch ÷ plan",
                    n ->
                        p -> speedup(
                            timing("branch", n, "transform", p, "transform"),
                            timing("branch", n, "transform_plan", p, "transform"),
                        ),
                ),
            ) if python || !startswith(heading, "Python")
        ],
    )
    M(version, kind) = n -> p -> peak(version, n, kind, p)
    table(
        io,
        "Peak memory, GB",
        problems,
        [
            python_column("Python", p -> peak("python", 0, "transform", p))
            [
                (heading * "/$n", entry(n)) for n ∈ threads for (heading, entry) ∈ (
                    ("main", M("main", "transform")),
                    ("branch", M("branch", "transform")),
                    ("plan", M("branch", "transform_plan")),
                )
            ]
        ],
    )
    table(
        io,
        "Costs paid once: seconds to build the TransformPlan, and seconds for the first\n" *
        "call in a new process, which compiles (on 400 time samples)",
        problems,
        [
            [
                (
                    "build plan/$n",
                    p -> seconds("branch", n, "transform_plan", p, "TransformPlan"),
                ) for n ∈ threads
            ]
            python_column(
                "first, Python", p -> information("python", 0, "transform", p, "first call")
            )
            [
                (
                    "first, $version/$n",
                    p -> information(version, n, "transform", p, "first call"),
                ) for n ∈ threads for version ∈ ("main", "branch")
            ]
        ],
    )

    println(io)
    println(io, "="^78)
    println(
        io,
        "3.  The gradient of the objective, comparing modes up to ℓ = $(options["lfit"])",
    )
    println(io, "="^78)
    G(n) = p -> timing("branch", n, "gradient", p, "gradient")
    table(
        io,
        "Seconds for one gradient: objective_and_gradient! on the branch, and finite\n" *
        "differences on main and in Python (forward differences, which take P + 1\n" *
        "evaluations of the objective)",
        problems,
        [
            [
                ("branch/$n", p -> seconds("branch", n, "gradient", p, "gradient")) for
                n ∈ threads
            ]
            [
                ("main/$n", p -> seconds("main", n, "differences", p, "differences")) for
                n ∈ threads
            ]
            python_column(
                "Python", p -> seconds("python", 0, "differences", p, "differences")
            )
        ],
    )
    table(
        io,
        "Ratios: each way's time for a gradient ÷ the branch's",
        problems,
        [
            (heading * "/$n", entry(n)) for n ∈ threads for (heading, entry) ∈ (
                (
                    "main ÷ branch",
                    n ->
                        p -> speedup(
                            timing("main", n, "differences", p, "differences"), G(n)(p)
                        ),
                ),
                (
                    "Python ÷ branch",
                    n ->
                        p -> speedup(
                            timing("python", 0, "differences", p, "differences"),
                            G(n)(p),
                        ),
                ),
            ) if python || !startswith(heading, "Python")
        ],
    )
    table(
        io,
        "Seconds for one evaluation of the objective: on the branch with its plan, and on\n" *
        "main and in Python as a transformation, an interpolation, and a sum; and the\n" *
        "branch's gradient as a multiple of its objective",
        problems,
        [
            [
                ("branch/$n", p -> seconds("branch", n, "gradient", p, "objective")) for
                n ∈ threads
            ]
            [
                ("main/$n", p -> seconds("main", n, "differences", p, "evaluation")) for
                n ∈ threads
            ]
            python_column(
                "Python", p -> seconds("python", 0, "differences", p, "evaluation")
            )
            [
                (
                    "gradient ÷ objective/$n",
                    p -> speedup(G(n)(p), timing("branch", n, "gradient", p, "objective")),
                ) for n ∈ threads
            ]
        ],
    )
    table(
        io,
        "Peak memory, GB",
        problems,
        [
            [("branch/$n", M("branch", "gradient")(n)) for n ∈ threads]
            [("main/$n", M("main", "differences")(n)) for n ∈ threads]
            python_column("Python", p -> peak("python", 0, "differences", p))
        ],
    )
    table(
        io,
        "Costs paid once on the branch: seconds to build the ObjectivePlan, and seconds for\n" *
        "the first call in a new process, which compiles (on 400 time samples)",
        problems,
        [
            [
                (
                    "build plan/$n",
                    p -> seconds("branch", n, "gradient", p, "ObjectivePlan"),
                ) for n ∈ threads
            ]
            [
                ("first/$n", p -> information("branch", n, "gradient", p, "first call")) for
                n ∈ threads
            ]
        ],
    )

    missing_cases = [name(c) for c ∈ cases if !haskey(results, name(c))]
    println(io)
    isempty(missing_cases) ||
        println(io, "Cases that failed or did not run: ", join(missing_cases, ", "))
    return nothing
end

function main(args)
    options, quick, report_only = parse_arguments(args)
    here = @__DIR__
    root = normpath(joinpath(here, "..", ".."))
    julia = Base.julia_cmd()

    # The working directory holds a copy of `main`, the two environments, and the results
    # of each case.  It is kept afterwards.
    work = isempty(options["work"]) ? mktempdir(; cleanup=false) : options["work"]
    println("Working in $work")
    mkpath(joinpath(work, "results"))
    main_copy = joinpath(work, "main")
    if !isdir(main_copy)
        mkpath(main_copy)
        run(pipeline(`git -C $root archive main`, `tar -x -C $main_copy`))
    end
    make_environment(julia, joinpath(work, "env-main"), main_copy)
    make_environment(julia, joinpath(work, "env-branch"), root)

    load = round(Sys.loadavg()[1]; digits=1)
    load > 2 && @warn "The machine is not quiet: its load average is $load."
    problems, cases = all_cases(options)
    began = time()
    if !report_only
        for (i, c) ∈ enumerate(cases)
            print("[$i/$(length(cases))] ")
            run_case(c, options, julia, work)
        end
    end
    minutes = round((time() - began) / 60; digits=1)

    git(args...) = strip(read(`git -C $root $args`, String))
    python_versions = ""
    blas = Dict{String,String}()
    for r ∈ values(read_results(cases, work))
        i = r.information
        if haskey(i, "versions")
            python_versions = "\n# Python scri: $(i["versions"]), with $(i["threads"]) threads for numba"
        else
            blas[i["threads"]] = i["BLAS threads"]
        end
    end
    blas_threads = join(
        ["$(blas[n]) with $n" for n ∈ split(options["threads"], ',') if haskey(blas, n)],
        ", ",
    )
    header = """
    # The headline benchmark ($(Dates.today())): time and peak memory for $(options["nt"]) time samples
    # $(Sys.cpu_info()[1].model), $(Sys.total_memory() ÷ 2^30) GB, Julia $VERSION
    # branch: $(git("rev-parse", "--abbrev-ref", "HEAD")) at $(git("rev-parse", "--short", "HEAD")); $(versions(julia, joinpath(work, "env-branch")))
    # main: at $(git("rev-parse", "--short", "main")); $(versions(julia, joinpath(work, "env-main")))$python_versions
    # Julia's threads: $(options["threads"]); BLAS threads, which keep their default: $blas_threads
    # Each case ran in a process of its own, one after another, in $minutes minutes.
    # Load average before the first case: $load
    # Each time is the least of several, in seconds, measured by Chairmarks."""
    report(stdout, options, header, problems, cases, work)
    if !quick && !report_only
        file = joinpath(root, "benchmark", "headline_$(Dates.today()).txt")
        isfile(file) &&
            (file = replace(file, ".txt" => Dates.format(Dates.now(), "_HHMM") * ".txt"))
        open(file, "w") do io
            return report(io, options, header, problems, cases, work)
        end
        println("\nSaved as $file")
    end
    return nothing
end
main(ARGS)
