#!/usr/bin/env bash
# Submit the cluster benchmark to SLURM.  Run it on the head node of mbot, from a clone of
# Scri.jl, in a shell where `julia` and `conda` both work:
#
#     benchmark/cluster/submit.sh [WORK]
#
# WORK is the working directory, which holds the two Julia environments, the logs, and the
# results; it defaults to a new directory under ~/scri-cluster-benchmark.  Giving the
# directory of an earlier run resumes it: the cases that finished are not run again.
#
# The jobs, each of which waits for the one before it, so that no two of them run at once:
#
#   1. `scri-setup`, with 4 cores, installs Scri 1.0.0 and 0.1.0 from the General registry
#      into the working directory, precompiles them on a compute node, and checks every
#      kind of case on a short problem;
#   2. `scri-c1`, with 1 core, compares 1.0 with 0.1 and Python, and times cold starts;
#   3. `scri-c2` through `scri-c192`, one for each number of cores, time 1.0 for the
#      scaling test; and
#   4. `scri-report` prints the report into its log, logs/report.out, and saves it in this
#      clone as benchmark/cluster_<date>.txt.
#
# With BLAS set (see below), the jobs of step 3 run instead the cases that set the number of
# BLAS threads, and step 2 is left out unless CORES includes 1.
#
# If the setup fails, SLURM cancels every job after it.  If a later job fails, the jobs
# after it still run, and the report lists what is missing.  `run.jl` says what each case
# measures.  To print the report again at any time, run, on the head node,
#
#     julia --startup-file=no benchmark/cluster/run.jl --report --work=WORK
#
# These environment variables change the defaults:
#
#     PARTITION   the SLURM partition (nd357_0001)
#     CORES       the numbers of cores of the jobs (1 2 4 8 16 32 64 96 192); 1 must be
#                 among them for the comparison
#     PYTHON      the Python with `scri` (the one in the conda environment `scri`)
#     JULIA       the Julia (the one on the PATH)
#     NT          the number of time samples (20000); NT=400 checks everything in minutes
#     BLAS        numbers of BLAS threads, such as "1 all", where `all` means the job's
#                 cores; if given, each job runs only the cases of 1.0 for the scaling test,
#                 once with each number, CORES defaults to 2 4 8 16 32 64 96 192, and the
#                 logs' names end in `-blas`.  Give the WORK of an earlier run, so that the
#                 report compares these cases with BLAS's default there.

set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
work=${1:-$HOME/scri-cluster-benchmark/$(date +%Y-%m-%d-%H%M)}
partition=${PARTITION:-nd357_0001}
if [[ -n ${BLAS:-} ]]; then
    cores=${CORES:-2 4 8 16 32 64 96 192}
    blas="--blas=${BLAS// /,}"
    tag=-blas
else
    cores=${CORES:-1 2 4 8 16 32 64 96 192}
    blas=
    tag=
fi
julia=${JULIA:-$(command -v julia)}
# Not `conda run`, which can find the active environment's Python first on the PATH.  A
# `conda` installed inside another environment may list `scri` only by its path, so a
# path ending in `/envs/scri` counts too.
if [[ -z ${PYTHON:-} ]]; then
    prefix=$(conda env list 2>/dev/null |
        awk '$1 == "scri" || $NF ~ /\/envs\/scri$/ { print $NF; exit }')
    PYTHON=${prefix:+$prefix/bin/python}
fi
python=${PYTHON:-}
nt=${NT:-20000}

if ! scri_version=$("$python" -c 'import scri; print(scri.__version__)' 2>/dev/null); then
    echo "No Python with \`scri\` was found (tried \"$python\"); give its path as PYTHON." >&2
    exit 1
fi

mkdir -p "$work/logs"
work=$(cd "$work" && pwd)
echo "Working in $work"
echo "Julia: $julia ($("$julia" --startup-file=no --version))"
echo "Python: $python (scri $scri_version)"

run="$julia --startup-file=no $here/run.jl --work=$work --python=$python --nt=$nt"

submit() {
    sbatch --parsable --partition="$partition" --nodes=1 --ntasks=1 \
        --kill-on-invalid-dep=yes --chdir="$work" "$@"
}

job=$(submit --job-name=scri-setup --cpus-per-task=4 --mem=16G --time=1:00:00 \
    --output="$work/logs/setup$tag.out" --wrap="$run --setup")
setup=$job
echo "setup: job $setup"

after="afterok:$setup"
for c in $cores; do
    if ((c == 1)); then
        # Python peaked at 11 GB on a laptop with the five Weyl components.
        memory=24G
        limit=3:00:00
    else
        memory=$((2 * c < 16 ? 16 : (2 * c > 64 ? 64 : 2 * c)))G
        limit=0:45:00
    fi
    job=$(submit --job-name="scri-c$c" --cpus-per-task="$c" --mem="$memory" \
        --time="$limit" --dependency="$after" --output="$work/logs/cores-$c$tag.out" \
        --wrap="$run --cores=$c $blas")
    echo "$c cores: job $job"
    after="afterok:$setup,afterany:$job"
done

job=$(submit --job-name=scri-report --cpus-per-task=1 --mem=4G --time=0:15:00 \
    --dependency="afterany:$job" --output="$work/logs/report$tag.out" --wrap="$run --report")
echo "report: job $job"
echo "Follow them with: squeue -u $USER"
