# The cluster benchmark: Scri.jl 1.0 against 0.1 and the original Python `scri` package,
# and how 1.0 scales with cores, each run as it would be run under SLURM on a cluster such
# as mbot.  `submit.sh` submits the jobs that run this script; see there.
#
# Every case has 20,000 time samples (by default), ℓₘₐₓ = 8, and a supertranslation with modes up to
# ℓα = 4, so that there are 31 parameters; the objective compares the modes up to ℓ = 8.
# The cases are those of `benchmark/headline/case.jl` and `case.py`, which describe what
# each one measures, and a cold start, by `cold.jl` and `cold.py`: a new process that loads
# the package, transforms the strain once, and exits.
#
# Both Julia versions are installed from the General registry, each into an environment
# of its own, as a user would install them.  Each case runs in a process of its own, inside
# a SLURM job that allocates the case's number of cores.  Julia is started with one thread
# for each of those cores, and BLAS and numba keep their own defaults, which follow the
# allocation: the variables that would override them are removed from each case's
# environment.
#
# BLAS's default is half the cores, so the matrix products of 1.0 run on half the threads
# that its loops do.  To see what that costs, `--blas` runs the cases of 1.0 for the
# scaling test again with BLAS's threads set by `OPENBLAS_NUM_THREADS`, and the report
# compares them with the default in a section of its own.
#
# Run it with `--work=DIR` and one of
#
#     --setup        build the two environments, and run every kind of case on a short
#                    problem, to check that each one works
#     --cores=N      run the cases for a job with N cores: with one core, the comparison
#                    of the three versions, the cold starts, and the first point of the
#                    scaling test; with more, the cases of 1.0 for the scaling test
#     --report       print the report of the cases that have finished, and save it as
#                    `benchmark/cluster_<date>.txt`
#
# Other options, each as `--name=value`:
#
#     --python=PATH  the Python executable of an environment that holds `scri`, which
#                    `--setup` and `--cores=1` need
#     --budget=10    seconds to spend repeating each measurement, as for `case.jl`
#     --limit=600    the longest that one gradient by finite differences may be expected
#                    to take before it is estimated instead, as for `case.jl`
#     --cold=3       cold starts of each version
#     --nt=20000     time samples; a small number, such as 400, checks the whole run in a
#                    few minutes
#     --blas=LIST    with `--cores=N`, run only the cases of 1.0 for the scaling test, once
#                    for each number of BLAS threads in the comma-separated list, where
#                    `all` means N; for example, `--blas=1,all`
#
# The results of each case are kept in `DIR/results`, in the format of the headline
# benchmark, and a case whose results are there already is not run again.

using Dates: Dates
using Printf: @sprintf

module Headline
include(joinpath(@__DIR__, "..", "headline", "run.jl"))
end

const ℓₘₐₓ, ℓα, ℓ_fit = 8, 4, 8
const strain = ("h",)
const all_components = ("h", "ψ₄", "ψ₃", "ψ₂", "ψ₁", "ψ₀")
const new, old = "v1.0.0", "v0.1.0"
const here = @__DIR__
const headline = normpath(joinpath(here, "..", "headline"))

# Variables that would override the numbers of threads that Julia, BLAS, and numba choose
# for themselves.  They are removed from the environment of every case.
const thread_variables = (
    "JULIA_NUM_THREADS",
    "JULIA_CPU_THREADS",
    "OPENBLAS_NUM_THREADS",
    "GOTO_NUM_THREADS",
    "OMP_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMBA_NUM_THREADS",
)
clean_environment() = Dict(k => v for (k, v) ∈ ENV if k ∉ thread_variables)

function parse_arguments(args)
    options = Dict(
        "work" => "",
        "python" => "",
        "cores" => "",
        "nt" => "20000",
        "budget" => "10",
        "limit" => "600",
        "cold" => "3",
        "blas" => "",
    )
    mode = ""
    for a ∈ args
        if a ∈ ("--setup", "--report")
            mode = a[3:end]
            continue
        end
        m = match(r"^--(\w+)=(.+)$", a)
        (isnothing(m) || !haskey(options, m[1])) && error("unknown option `$a`")
        options[m[1]] = m[2]
    end
    isempty(options["cores"]) || (mode = "cores")
    isempty(mode) && error("give one of `--setup`, `--cores=N`, or `--report`")
    isempty(options["work"]) && error("give the working directory as `--work=DIR`")
    return mode, options
end

# A case, with the number of BLAS threads it sets, or `nothing` for BLAS's default.
function case(version, threads, kind, components; blas=nothing)
    return (; version, threads, kind, ℓₘₐₓ, components, ℓα, blas)
end

# The name of a case's results, which is that of the headline benchmark, and ends in
# `-b<threads>` for a case that sets the number of BLAS threads.
name(c) = Headline.name(c) * (isnothing(c.blas) ? "" : "-b$(c.blas)")

# The cases for a job with `cores` cores, in the order they run.
function cases(cores)
    scaling = [
        case(new, cores, kind, components) for kind ∈ ("transform", "gradient") for
        components ∈ (strain, all_components)
    ]
    cores > 1 && return scaling
    return [
        scaling
        case(old, 1, "transform", strain)
        case(old, 1, "transform", all_components)
        case(old, 1, "differences", strain)
        case("python", 1, "transform", strain)
        case("python", 1, "transform", all_components)
        case("python", 1, "differences", strain)
        [case(version, 1, "cold", strain) for version ∈ (new, old, "python")]
    ]
end

# The cases of 1.0 for the scaling test, with each of the numbers of BLAS threads in
# `blas`, a comma-separated list in which `all` means `cores`.  A number equal to BLAS's
# default is run too, as a repeat of the default case.
function cases(cores, blas)
    settings = unique(
        s == "all" ? cores : parse(Int, s) for s ∈ split(blas, ',') if !isempty(s)
    )
    all(≥(1), settings) || error("the numbers of BLAS threads must be positive: $blas")
    return [
        case(new, cores, kind, components; blas=n) for n ∈ settings for
        kind ∈ ("transform", "gradient") for components ∈ (strain, all_components)
    ]
end

################################################################################
# Running the cases

# A new environment that holds the registered `version` of Scri, and Chairmarks.  The
# registry is updated first, since a copy that is out of date may not list 1.0.0.  A file
# named `ready` marks an environment that was finished.
function make_environment(julia, env, version)
    isfile(joinpath(env, "ready")) && return nothing
    mkpath(env)
    code = """
    using Pkg
    Pkg.Registry.update()
    Pkg.add(; name="Scri", version=ARGS[1])
    Pkg.add("Chairmarks")
    Pkg.precompile()
    using Scri, Chairmarks
    """
    run(
        setenv(
            `$julia --startup-file=no --project=$env -e $code $version`, clean_environment()
        ),
    )
    touch(joinpath(env, "ready"))
    return nothing
end

function command(c, args, options, julia, work)
    directory, script = c.kind == "cold" ? (here, "cold") : (headline, "case")
    cmd = if c.version == "python"
        `$(options["python"]) $(joinpath(directory, script * ".py")) $args`
    else
        env = joinpath(work, "env-" * c.version)
        `$julia --startup-file=no --project=$env -t $(c.threads) $(joinpath(directory, script * ".jl")) $args`
    end
    environment = clean_environment()
    isnothing(c.blas) || (environment["OPENBLAS_NUM_THREADS"] = string(c.blas))
    return setenv(cmd, environment)
end

# The seconds that the command took from launch to exit, or nothing if it failed.
function run_process(cmd)
    seconds = @elapsed ok = success(pipeline(cmd; stdout, stderr))
    return ok ? seconds : nothing
end

# Lines of a results file, as `case.jl` writes them.
info_line(name, value) = join(("info", name, value), '\t')
function time_line(name, seconds)
    s = sort(seconds)
    n = length(s)
    median = isodd(n) ? s[(n + 1) ÷ 2] : (s[n ÷ 2] + s[n ÷ 2 + 1]) / 2
    return join(("time", name, s[1], median, n, 0, 0, "measured"), '\t')
end

# Where and how the case ran.
function context(c)
    return [
        info_line("host", gethostname())
        info_line("cpu", strip(Sys.cpu_info()[1].model))
        info_line("job", get(ENV, "SLURM_JOB_ID", "none"))
        info_line("cores", c.threads)
    ]
end

# Run one case, unless its results are in `directory` already, and say whether it worked.
function run_case(c, options, julia, work, directory, Nᵗ)
    file = joinpath(directory, name(c) * ".tsv")
    isfile(file) && return true
    println("$(Dates.format(Dates.now(), "HH:MM:SS"))  $(name(c))")
    flush(stdout)
    began = time()
    ok = if c.kind == "cold"
        run_cold(c, options, julia, work, file, Nᵗ)
    else
        run_timings(c, options, julia, work, file, Nᵗ)
    end
    ok || @warn "The case $(name(c)) failed, and is left out of the report."
    println("    $(round(time() - began; digits=1)) s in all")
    flush(stdout)
    return ok
end

function run_timings(c, options, julia, work, file, Nᵗ)
    part = file * ".part"
    args = [
        c.kind,
        string(c.ℓₘₐₓ),
        join(c.components, ','),
        string(c.ℓα),
        string(ℓ_fit),
        string(Nᵗ),
        options["budget"],
        options["limit"],
        part,
    ]
    seconds = run_process(command(c, args, options, julia, work))
    if isnothing(seconds)
        isfile(part) && mv(part, file * ".failed"; force=true)
        return false
    end
    open(part, "a") do io
        return foreach(l -> println(io, l), [context(c); info_line("process", seconds)])
    end
    mv(part, file)
    return true
end

# Start the process `--cold` times, each time timing it from launch to exit.  The first
# start in a job may be the first on its node, and so the slowest, so its time is kept as
# well as the least.
function run_cold(c, options, julia, work, file, Nᵗ)
    process, load, first_call, memory = Float64[], Float64[], Float64[], Float64[]
    for i ∈ 1:parse(Int, options["cold"])
        part = "$file.part$i"
        args = [string(c.ℓₘₐₓ), join(c.components, ','), string(c.ℓα), string(Nᵗ), part]
        seconds = run_process(command(c, args, options, julia, work))
        isnothing(seconds) && return false
        information = Headline.read_case(part).information
        push!(process, seconds)
        push!(load, parse(Float64, information["load"]))
        push!(first_call, parse(Float64, information["first call at full size"]))
        push!(memory, parse(Float64, information["memory after one call"]))
        rm(part)
    end
    lines = [
        context(c)
        info_line("process, first start", process[1])
        info_line("memory after one call", maximum(memory))
        time_line("process", process)
        time_line("load", load)
        time_line("first call at full size", first_call)
    ]
    open(io -> foreach(l -> println(io, l), lines), file, "w")
    return true
end

function run_setup(options, julia, work)
    isempty(options["python"]) && error("`--setup` needs `--python=PATH`")
    make_environment(julia, joinpath(work, "env-" * new), "1.0.0")
    make_environment(julia, joinpath(work, "env-" * old), "0.1.0")
    check = mkpath(joinpath(work, "check"))
    short = merge(options, Dict("budget" => "0", "limit" => "5", "cold" => "1"))
    threads = min(4, Sys.EFFECTIVE_CPU_THREADS)
    failed = [
        name(c) for c ∈ unique([cases(1); cases(threads)]) if
        !run_case(c, short, julia, work, check, 400)
    ]
    isempty(failed) || error("These cases failed on a short problem: $(join(failed, ", "))")
    println("Every kind of case works on a short problem.")
    return nothing
end

function run_cores(options, julia, work)
    cores = parse(Int, options["cores"])
    allocated = get(ENV, "SLURM_CPUS_PER_TASK", "")
    isempty(allocated) ||
        parse(Int, allocated) == cores ||
        error("SLURM allocated $allocated cores to this job, not $cores")
    Sys.EFFECTIVE_CPU_THREADS == cores ||
        @warn "Julia sees $(Sys.EFFECTIVE_CPU_THREADS) cores available, not $cores, so the " *
            "defaults of BLAS and numba may not follow the allocation."
    selected = isempty(options["blas"]) ? cases(cores) : cases(cores, options["blas"])
    if cores == 1 && isempty(options["python"])
        @warn "Without `--python`, the Python cases are left out."
        filter!(c -> c.version != "python", selected)
    end
    directory = mkpath(joinpath(work, "results"))
    failures = 0
    for (i, c) ∈ enumerate(selected)
        print("[$i/$(length(selected))] ")
        run_case(c, options, julia, work, directory, parse(Int, options["nt"])) ||
            (failures += 1)
    end
    failures == 0 || error("$failures of the $(length(selected)) cases failed")
    return nothing
end

################################################################################
# The report

# The words of `text` filled into lines of at most `width` characters, with two spaces after
# the end of each sentence.
function fill_lines(text; width=80)
    lines = [""]
    for word ∈ split(text)
        separator = if isempty(lines[end])
            ""
        elseif endswith(lines[end], '.')
            "  "
        else
            " "
        end
        if textwidth(lines[end]) + textwidth(separator) + textwidth(word) ≤ width
            lines[end] *= separator * word
        else
            push!(lines, word)
        end
    end
    return join(lines, '\n')
end

# A table whose first column is left-aligned.
function grid(io, title, headings, rows)
    widths = [
        maximum(textwidth, [headings[i]; [r[i] for r ∈ rows]]) for i ∈ eachindex(headings)
    ]
    line(cells) = join(
        [
            i == 1 ? rpad(cells[i], widths[i]) : lpad(cells[i], widths[i]) for
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

function report(io, options, julia, work)
    directory = joinpath(work, "results")
    results = Dict(
        f[1:(end - 4)] => Headline.read_case(joinpath(directory, f)) for
        f ∈ readdir(directory) if endswith(f, ".tsv")
    )
    result(c) = get(results, name(c), nothing)
    function timing(c, quantity)
        r = result(c)
        return isnothing(r) ? nothing : get(r.timings, quantity, nothing)
    end
    function information(c, key)
        r = result(c)
        return isnothing(r) ? nothing : get(r.information, key, nothing)
    end
    function seconds(x)
        isnothing(x) && return "—"
        return (x.how == "estimated" ? "~" : "") * Headline.digits3(x.seconds)
    end
    function quotient(a, b)
        (isnothing(a) || isnothing(b)) && return "—"
        return (a.how == "estimated" ? "~" : "") * Headline.ratio(a.seconds / b.seconds)
    end
    # The peak memory of the process after its first call at full size, in gibibytes.
    function peak(c)
        r = result(c)
        isnothing(r) && return "—"
        i = r.information
        bytes = get(i, "memory after one gradient", get(i, "memory after one call", ""))
        return isempty(bytes) ? "—" : @sprintf("%.2f", parse(Float64, bytes) / 2^30)
    end
    versions = ("python", old, new)

    hosts = sort(
        unique(
            r.information["host"] for r ∈ values(results) if haskey(r.information, "host")
        ),
    )
    cpus = unique(
        r.information["cpu"] for r ∈ values(results) if haskey(r.information, "cpu")
    )
    python = something(information(case("python", 1, "transform", strain), "versions"), "—")
    numba = something(information(case("python", 1, "transform", strain), "threads"), "—")
    blas = something(information(case(new, 1, "transform", strain), "BLAS threads"), "—")
    P = 6 + (ℓα + 1)^2
    package_versions(label) = Headline.versions(julia, joinpath(work, "env-" * label))
    println(
        io,
        """
        # The cluster benchmark ($(Dates.today())): Scri.jl 1.0 and 0.1, and Python scri, under SLURM
        # $(join(cpus, "; ")), on $(join(hosts, ", ")); Julia $VERSION
        # 1.0: Scri 1.0.0; $(package_versions(new))
        # 0.1: Scri 0.1.0; $(package_versions(old))
        # Python: $python
        # Every case: $(options["nt"]) time samples, ℓₘₐₓ = $ℓₘₐₓ, a supertranslation with modes up to ℓ = $ℓα
        # ($P parameters), and an objective that compares the modes up to ℓ = $ℓ_fit.
        """,
    )
    paragraph = """
        Each case ran in a process of its own, inside a SLURM job that allocated it the
        cores shown, one job after another.  Julia was started with one thread for each
        allocated core, and BLAS and numba chose their own numbers of threads: with one
        core, BLAS had $blas and numba had $numba.  Each time is the least of several, in
        seconds.  The gradient on 0.1 and in Python is by forward differences, which take
        $(P + 1) evaluations of the objective; a time marked ~ was not measured, but is
        $(P + 1) times one evaluation, because the gradient would have taken longer than
        $(options["limit"]) seconds.  Memory is the peak resident set size of the process just
        after its first call at full size, in gibibytes (GiB, 2³⁰ bytes)."""
    println(io, fill_lines(paragraph))

    rows = [
        (
            "transform!, h",
            [(case(v, 1, "transform", strain), "transform") for v ∈ versions],
        ),
        (
            "transform!, h ψ₄ … ψ₀",
            [(case(v, 1, "transform", all_components), "transform") for v ∈ versions],
        ),
        (
            "objective, h",
            [
                (case("python", 1, "differences", strain), "evaluation"),
                (case(old, 1, "differences", strain), "evaluation"),
                (case(new, 1, "gradient", strain), "objective"),
            ],
        ),
        (
            "gradient, h",
            [
                (case("python", 1, "differences", strain), "differences"),
                (case(old, 1, "differences", strain), "differences"),
                (case(new, 1, "gradient", strain), "gradient"),
            ],
        ),
    ]
    grid(
        io,
        "==============================================================================\n" *
        "1.  One core: seconds for one call, and their ratios\n" *
        "==============================================================================\n",
        ["", "Python", "0.1", "1.0", "0.1 ÷ 1.0", "Python ÷ 1.0"],
        [
            [
                label
                [seconds(timing(w...)) for w ∈ ways]
                quotient(timing(ways[2]...), timing(ways[3]...))
                quotient(timing(ways[1]...), timing(ways[3]...))
            ] for (label, ways) ∈ rows
        ],
    )
    grid(
        io,
        "Peak memory, GiB",
        ["", "Python", "0.1", "1.0"],
        [
            [label; [peak(first(w)) for w ∈ ways]] for
            (label, ways) ∈ rows if label != "objective, h"
        ],
    )

    cold = [case(v, 1, "cold", strain) for v ∈ versions]
    first_start(c) =
        let x = information(c, "process, first start")
            isnothing(x) ? "—" : Headline.digits3(parse(Float64, x))
        end
    starts = something(timing(cold[3], "process"), (; n=0)).n
    grid(
        io,
        "==============================================================================\n" *
        "2.  A cold start: seconds for a new process that loads the package, transforms\n" *
        "    the strain once, and exits, with the least of $starts starts\n" *
        "==============================================================================\n",
        ["", "Python", "0.1", "1.0"],
        [
            ["the first start in the job"; first_start.(cold)],
            ["the least"; [seconds(timing(c, "process")) for c ∈ cold]],
            ["of which, loading"; [seconds(timing(c, "load")) for c ∈ cold]],
            [
                "of which, the transformation";
                [seconds(timing(c, "first call at full size")) for c ∈ cold]
            ],
            ["peak memory, GiB"; peak.(cold)],
        ],
    )

    core_counts = sort(
        unique(
            parse(Int, match(r"-t(\d+)-", k)[1]) for
            k ∈ keys(results) if startswith(k, new * "-t") && !occursin(r"-b\d+$", k)
        ),
    )
    kinds = [
        ("transform!, h", "transform", strain, "transform"),
        ("transform!, h ψ₄ … ψ₀", "transform", all_components, "transform"),
        ("gradient, h", "gradient", strain, "gradient"),
        ("gradient, h ψ₄ … ψ₀", "gradient", all_components, "gradient"),
    ]
    blas_threads(n) =
        something(information(case(new, n, "transform", strain), "BLAS threads"), "—")
    grid(
        io,
        "==============================================================================\n" *
        "3.  Scaling of 1.0: seconds for one call with the cores of the job, Julia's\n" *
        "    threads, and BLAS's threads\n" *
        "==============================================================================\n",
        ["cores"; "BLAS"; first.(kinds)],
        [
            [
                string(n);
                blas_threads(n);
                [seconds(timing(case(new, n, k, cs), q)) for (_, k, cs, q) ∈ kinds]
            ] for n ∈ core_counts
        ],
    )
    grid(
        io,
        "Speedup over one core",
        ["cores"; first.(kinds)],
        [
            [
                string(n)
                [
                    quotient(
                        timing(case(new, 1, k, cs), q), timing(case(new, n, k, cs), q)
                    ) for (_, k, cs, q) ∈ kinds
                ]
            ] for n ∈ core_counts
        ],
    )

    # Section 4, if any case set the number of BLAS threads: at each number of cores, the
    # default case and each number of BLAS threads that was set.  The number of BLAS
    # threads shown is the one the case reported, not the one it asked for.
    variants = Tuple{Int,Int}[]
    for k ∈ keys(results)
        m = match(r"-t(\d+)-.*-b(\d+)$", k)
        startswith(k, new * "-t") &&
            !isnothing(m) &&
            push!(variants, (parse(Int, m[1]), parse(Int, m[2])))
    end
    sort!(unique!(variants))
    if !isempty(variants)
        blas_row(n, label, blas) = [
            string(n)
            something(
                information(case(new, n, "transform", strain; blas), "BLAS threads"),
                "—",
            ) * label
            [seconds(timing(case(new, n, k, cs; blas), q)) for (_, k, cs, q) ∈ kinds]
        ]
        rows = []
        for n ∈ unique(first.(variants))
            push!(rows, blas_row(n, " (default)", nothing))
            append!(rows, [blas_row(n, "", b) for (m, b) ∈ variants if m == n])
        end
        grid(
            io,
            "==============================================================================\n" *
            "4.  BLAS's threads and 1.0: seconds for one call with the cores of the job,\n" *
            "    Julia's threads, and the BLAS threads shown\n" *
            "==============================================================================\n",
            ["cores"; "BLAS"; first.(kinds)],
            rows,
        )
    end

    failed = [f for f ∈ readdir(directory) if endswith(f, ".failed")]
    isempty(failed) || println(io, "\nCases that failed: ", join(failed, ", "))
    return nothing
end

function main(args)
    mode, options = parse_arguments(args)
    julia = Base.julia_cmd()
    work = abspath(options["work"])
    println("Working in $work")
    mode == "setup" && return run_setup(options, julia, work)
    mode == "cores" && return run_cores(options, julia, work)
    report(stdout, options, julia, work)
    root = normpath(joinpath(here, "..", ".."))
    file = joinpath(root, "benchmark", "cluster_$(Dates.today()).txt")
    isfile(file) &&
        (file = replace(file, ".txt" => Dates.format(Dates.now(), "_HHMM") * ".txt"))
    open(io -> report(io, options, julia, work), file, "w")
    println("\nSaved as $file")
    return nothing
end
main(ARGS)
