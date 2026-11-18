# A cold start of the original Python `scri` package: a new process that imports it,
# transforms the data once at full size, and exits, as a job on a cluster does.  Run by
# `run.jl`, which times the whole process from launch to exit; see there.
#
#     python cold.py <ℓmax> <components> <ℓα> <Nt> <file>
#
# The data and the transformation are those of `benchmark/headline/case.py`, whose
# functions this script uses.  Importing it imports NumPy, `quaternion`, and `scri`, and
# nothing else of note, so the time to import it is the time to import `scri`.  The
# results file has the lines that `cold.jl` writes.

import os
import sys
import time

process_start = time.perf_counter()
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "headline"))
import case  # noqa: E402

load_seconds = time.perf_counter() - process_start


def main(args):
    ℓmax, components, ℓα, Nt, file = args
    data = case.make_data(int(ℓmax), int(Nt), tuple(components.split(",")))
    θ = case.parameters(int(ℓα))
    seconds = case.sample(lambda: case.transform(data, θ))
    results = case.Results(file)
    results.info("load", load_seconds)
    results.info("first call at full size", seconds)
    results.memory("after one call")


if __name__ == "__main__":
    main(sys.argv[1:])
