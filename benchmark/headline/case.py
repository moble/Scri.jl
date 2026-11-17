# One case of the headline benchmark for the original Python `scri` package, run in a
# process of its own, so that the peak memory of the process belongs to this case alone.
# Run by `run.jl` when it is given `--python`; see there.  It needs an environment that
# holds `scri`, such as
#
#     mamba create -n scri-timing -c conda-forge --override-channels scri
#
# The arguments are those of `case.jl`.  The kinds are
#
#   - `transform`: one call of the package's BMS transformation; and
#   - `differences`: the gradient of an objective by finite differences (the forward kind,
#     with one evaluation for each parameter and one more), each evaluation being a
#     transformation, an interpolation onto the target's fixed times with the package's
#     own `interpolate`, and a sum of squares over the modes up to `ℓ_fit`.
#
# The strain alone is transformed as a `WaveformModes` object.  With Weyl components, the
# data are an `AsymptoticBondiData` object, which is what `scri` recommends, and which
# always transforms the shear and all five Weyl components together.  Every option of the
# transformation is left at its default, except that the `AsymptoticBondiData` truncates
# products at the larger band limit, as its documentation advises.
#
# The keywords of `scri`'s `transform` describe the same transformation as the same
# keywords of Scri.jl's `BMS`; a cross-check on real waveforms on 2026-11-17 found no
# difference of sign or of order.  The supertranslation's modes here are drawn with
# NumPy's generator rather than Julia's, however, so the transformation is not exactly the
# one in `case.jl`.  It is one of the same size, which is all that the timing needs.  The
# data are random, as there.
#
# A gradient that would take longer than the limit is not run: its time is the number of
# evaluations times the time of one, and it is marked as an estimate.

import gc
import os
import resource
import statistics
import sys
import time

process_start = time.perf_counter()
import numpy as np
import quaternion
import scri

load_seconds = time.perf_counter() - process_start

Nt_short = 400
spin_weights = {"h": -2, "ψ₄": -2, "ψ₃": -1, "ψ₂": 0, "ψ₁": 1, "ψ₀": 2}
abd_names = {"ψ₄": "psi4", "ψ₃": "psi3", "ψ₂": "psi2", "ψ₁": "psi1", "ψ₀": "psi0"}


def parameters(ℓα):
    """The parameters of a small but generic BMS transformation whose supertranslation has
    modes up to `ℓα`, laid out as for `bms_from_parameters` in Scri.jl."""
    rng = np.random.default_rng(20261116)
    return np.concatenate(
        ([0.1, -0.2, 0.3, 0.02, 0.01, -0.03], 0.02 * rng.standard_normal((ℓα + 1) ** 2))
    )


def bms(θ):
    """The rotor (as its four components, which is the form that `scri` reads), the boost
    velocity, and the supertranslation modes for the parameters θ."""
    θ = np.asarray(θ, dtype=float)
    R = quaternion.as_float_array(quaternion.from_rotation_vector(θ[0:3]))
    η = np.linalg.norm(θ[3:6])
    v = (np.tanh(η) / η if η > 0 else 1.0) * θ[3:6]
    a = θ[6:]
    L = int(round(np.sqrt(a.size)))
    α = np.zeros(a.size, dtype=complex)
    for ℓ in range(L):
        i0 = ℓ * ℓ + ℓ
        α[i0] = a[i0]
        for m in range(1, ℓ + 1):
            αm = a[i0 + m] + 1j * a[i0 - m]
            α[i0 + m] = αm
            α[i0 - m] = (-1) ** m * np.conj(αm)
    return R, v, α


def random_modes(rng, Nt, ℓmax, s):
    """Random modes of spin weight `s` from ℓ = 0, with those below |s| set to zero."""
    shape = (Nt, (ℓmax + 1) ** 2)
    d = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    d[:, : s * s] = 0
    return d


def make_data(ℓmax, Nt, components):
    """Random data with `Nt` time samples at a step of 0.5."""
    rng = np.random.default_rng(1)
    t = 0.5 * np.arange(Nt, dtype=float)
    if components == ("h",):
        return scri.WaveformModes(
            t=t,
            data=random_modes(rng, Nt, ℓmax, -2)[:, 4:],
            ell_min=2,
            ell_max=ℓmax,
            frameType=scri.Inertial,
            dataType=scri.h,
            r_is_scaled_out=True,
            m_is_scaled_out=True,
        )
    abd = scri.AsymptoticBondiData(t, ℓmax, multiplication_truncator=max)
    abd.sigma = random_modes(rng, Nt, ℓmax, 2)
    for name, s in (("psi4", -2), ("psi3", -1), ("psi2", 0), ("psi1", 1), ("psi0", 2)):
        setattr(abd, name, random_modes(rng, Nt, ℓmax, s))
    return abd


def transform(data, θ):
    R, v, α = bms(θ)
    return data.transform(supertranslation=α, frame_rotation=R, boost_velocity=v)


class Objective:
    """The sum of squares of the difference between the transformed data and a random
    target, over the modes up to `ℓ_fit`.  The target's times are fixed by the first
    evaluation, inside the output times of the transformation."""

    def __init__(self, ℓmax, Nt, components, ℓ_fit):
        rng = np.random.default_rng(2)
        self.components = components
        self.Nt = Nt
        self.columns = (ℓ_fit + 1) ** 2
        self.data = make_data(ℓmax, Nt, components)
        self.t_target = None
        self.target = [
            random_modes(rng, Nt, ℓ_fit, spin_weights[c])[:, spin_weights[c] ** 2 :]
            for c in components
        ]

    def __call__(self, θ):
        d = transform(self.data, θ)
        if self.t_target is None:
            self.t_target = np.linspace(d.t[0] + 1, d.t[-1] - 1, self.Nt)
        d = d.interpolate(self.t_target)
        if self.components == ("h",):
            # The modes of a `WaveformModes` object start at ℓ = 2.
            return np.sum(np.abs(d.data[:, : self.columns - 4] - self.target[0]) ** 2)
        total = 0.0
        for c, target in zip(self.components, self.target):
            s = spin_weights[c]
            # The strain is twice the conjugate of the shear, as in `AsymptoticBondiData.h`.
            modes = np.asarray(2.0 * d.sigma.bar if c == "h" else getattr(d, abd_names[c]))
            total += np.sum(np.abs(modes[:, s * s : self.columns] - target) ** 2)
        return total


def forward_differences(f, θ, h=1e-6):
    f0 = f(θ)
    gradient = np.empty_like(θ)
    for p in range(θ.size):
        θp = θ.copy()
        θp[p] += h
        gradient[p] = (f(θp) - f0) / h
    return gradient


def sample(f):
    """The time of one call of `f`."""
    start = time.perf_counter()
    f()
    return time.perf_counter() - start


def repeat_samples(samples, f, budget):
    """Add timed calls until there are at least three of them (two, if a call takes longer
    than 100 seconds), and then until `budget` seconds have gone by or there are 30."""
    start = time.perf_counter()
    while len(samples) < 30:
        n = len(samples)
        enough = n >= 3 or (n == 2 and min(samples) > 100)
        if enough and time.perf_counter() - start >= budget:
            break
        samples.append(sample(f))
    return samples


class Results:
    """The lines of the results file, which `run.jl` reads back, as `case.jl` writes them."""

    def __init__(self, file):
        self.io = open(file, "w", encoding="utf-8")

    def info(self, name, value):
        print("info", name, value, sep="\t", file=self.io, flush=True)

    def memory(self, stage):
        # The peak resident set size, which macOS reports in bytes and Linux in kilobytes.
        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        self.info(f"memory {stage}", peak * (1 if sys.platform == "darwin" else 1024))

    def record(self, name, samples):
        row = ("time", name, min(samples), statistics.median(samples), len(samples), 0, 0)
        print(*row, "measured", sep="\t", file=self.io, flush=True)
        print(f"    {name}: {min(samples):.4g} s, least of {len(samples)}", flush=True)

    def record_estimate(self, name, seconds):
        print("time", name, seconds, seconds, 0, 0, 0, "estimated", sep="\t", file=self.io, flush=True)
        print(f"    {name}: about {seconds:.4g} s, estimated", flush=True)


def run_transform(results, ℓmax, components, ℓα, Nt, budget):
    θ = parameters(ℓα)
    results.info("first call", sample(lambda: transform(make_data(ℓmax, Nt_short, components), θ)))
    data = make_data(ℓmax, Nt, components)
    gc.collect()
    results.memory("before")
    samples = [sample(lambda: transform(data, θ))]
    results.memory("after one call")
    repeat_samples(samples, lambda: transform(data, θ), budget)
    results.record("transform", samples)


def run_differences(results, ℓmax, components, ℓα, ℓ_fit, Nt, budget, limit):
    θ = parameters(ℓα)
    results.info("first call", sample(lambda: Objective(ℓmax, Nt_short, components, ℓ_fit)(θ)))
    f = Objective(ℓmax, Nt, components, ℓ_fit)
    gc.collect()
    results.memory("before")
    samples = [sample(lambda: f(θ))]
    results.memory("after one call")
    repeat_samples(samples, lambda: f(θ), budget)
    results.record("evaluation", samples)
    predicted = (θ.size + 1) * min(samples)
    if predicted > limit:
        results.record_estimate("differences", predicted)
    else:
        # One gradient is enough: it is itself many evaluations.
        samples = [sample(lambda: forward_differences(f, θ))]
        results.memory("after one gradient")
        results.record("differences", samples)


def main(args):
    kind, ℓmax, components, ℓα, ℓ_fit, Nt, budget, limit, file = args
    ℓmax, ℓα, ℓ_fit, Nt = int(ℓmax), int(ℓα), int(ℓ_fit), int(Nt)
    budget, limit = float(budget), float(limit)
    components = tuple(components.split(","))
    import numba

    results = Results(file)
    results.info("load average", os.getloadavg()[0])
    results.info("threads", numba.get_num_threads())
    results.info(
        "versions",
        f"Python {sys.version.split()[0]}; scri {scri.__version__}; numpy {np.__version__}; "
        f"numba {numba.__version__}",
    )
    results.info("load", load_seconds)
    if kind == "transform":
        run_transform(results, ℓmax, components, ℓα, Nt, budget)
    elif kind == "differences":
        run_differences(results, ℓmax, components, ℓα, ℓ_fit, Nt, budget, limit)
    else:
        raise ValueError(f"unknown kind `{kind}`")
    results.memory("at the end")


if __name__ == "__main__":
    main(sys.argv[1:])
