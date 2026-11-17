# The part of the comparison that runs the original Python `scri` package: the time of one
# evaluation of an objective built from its BMS transformation, from which the time of a
# gradient by finite differences follows as a number of evaluations.  Run by `run.jl` when
# it is given `--python`; see there.  It needs an environment that holds `scri`, such as
#
#     mamba create -n scri-timing -c conda-forge --override-channels scri
#
# As on `main`, the transformation chooses its own output times, so one evaluation is a
# transformation, an interpolation onto the target's fixed times, and a sum of squares.
#
# The strain alone is transformed as a `WaveformModes` object.  The sets with Weyl
# components are transformed as an `AsymptoticBondiData` object, which is what `scri`
# recommends, and which always transforms the shear and all five Weyl components together.
# So the two larger sets cost the same here, and the second is not timed again: its rows
# repeat the first's.  Every option of the transformation is left at
# its default, except that the `AsymptoticBondiData` truncates products at the larger band
# limit, as its documentation advises.
#
# The conventions of `scri` differ from those of Scri.jl in the sign of the
# supertranslation and in the order of the factors, so the parameters do not describe the
# same transformation as in the other scripts.  It is one of the same size, which is all
# that the timing needs.  The data are random, as in the other scripts.
#
# The gradients themselves are never run.  One evaluation is timed, and if even that would
# take longer than the limit, or more memory than the machine has, it is timed with fewer
# time samples and scaled up, and marked as an estimate.

import os
import sys
import time

import numpy as np
import quaternion
import scri

component_sets = (("h",), ("h", "ψ₄", "ψ₃", "ψ₂"), ("h", "ψ₄", "ψ₃", "ψ₂", "ψ₁", "ψ₀"))
spin_weights = {"h": -2, "ψ₄": -2, "ψ₃": -1, "ψ₂": 0, "ψ₁": 1, "ψ₀": 2}
abd_names = {"ψ₄": "psi4", "ψ₃": "psi3", "ψ₂": "psi2", "ψ₁": "psi1", "ψ₀": "psi0"}


def bms(θ):
    """The rotor (as its four components, which is the form that `scri` reads), the boost
    velocity, and the supertranslation modes for the parameters θ, which are laid out as
    for `bms_from_parameters` in Scri.jl."""
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


class Problem:
    """One objective: random data with `Nt` time samples, and a random target on fixed
    times inside the output times of the transformation."""

    def __init__(self, ℓmax, Nt, components, ℓα):
        rng = np.random.default_rng(20261116)
        self.components = components
        self.t = 0.5 * np.arange(Nt, dtype=float)
        self.θ = np.concatenate(
            ([0.1, -0.2, 0.3, 0.02, 0.01, -0.03], 0.02 * rng.standard_normal((ℓα + 1) ** 2))
        )
        if components == ("h",):
            self.data = scri.WaveformModes(
                t=self.t,
                data=random_modes(rng, Nt, ℓmax, -2)[:, 4:],
                ell_min=2,
                ell_max=ℓmax,
                frameType=scri.Inertial,
                dataType=scri.h,
                r_is_scaled_out=True,
                m_is_scaled_out=True,
            )
        else:
            abd = scri.AsymptoticBondiData(self.t, ℓmax, multiplication_truncator=max)
            abd.sigma = random_modes(rng, Nt, ℓmax, 2)
            for name, s in (("psi4", -2), ("psi3", -1), ("psi2", 0), ("psi1", 1), ("psi0", 2)):
                setattr(abd, name, random_modes(rng, Nt, ℓmax, s))
            self.data = abd
        t_out = self.transformed(self.θ).t
        self.t_target = np.linspace(t_out[0] + 1, t_out[-1] - 1, Nt)
        self.target = [
            random_modes(rng, Nt, ℓmax, spin_weights[c])[:, spin_weights[c] ** 2 :]
            for c in components
        ]

    def transformed(self, θ):
        R, v, α = bms(θ)
        return self.data.transform(supertranslation=α, frame_rotation=R, boost_velocity=v)

    def objective(self, θ):
        d = self.transformed(θ).interpolate(self.t_target)
        if self.components == ("h",):
            return np.sum(np.abs(d.data - self.target[0]) ** 2)
        total = 0.0
        for c, target in zip(self.components, self.target):
            s = spin_weights[c]
            # The strain is twice the conjugate of the shear, as in `AsymptoticBondiData.h`.
            modes = np.asarray(2.0 * d.sigma.bar if c == "h" else getattr(d, abd_names[c]))
            total += np.sum(np.abs(modes[:, s * s :] - target) ** 2)
        return total



def least_time(f, budget):
    """The least of several timings of `f()`: at least three, and more until `budget`
    seconds have passed, unless the first alone takes longer than the budget."""
    start = time.perf_counter()
    f()
    best = time.perf_counter() - start
    n = 1
    nmin = 1 if best > budget else 3
    began = time.perf_counter()
    while n < nmin or (n < 50 and time.perf_counter() - began < budget):
        start = time.perf_counter()
        f()
        best = min(best, time.perf_counter() - start)
        n += 1
    return best, n


def grid_bytes(ℓmax, Nt, components):
    """A rough count of the memory that the transformation needs: a few copies of every
    field on the grid of directions, at every time."""
    if components == ("h",):
        fields, n = 1, 2 * (ℓmax + 1) + 1
    else:
        fields, n = 6, 2 * (2 * ℓmax) + 1
    return 4 * fields * n * n * Nt * 16


# The time of one evaluation with Weyl components, for each band limit of the data and of
# the supertranslation, so that the second set with Weyl components can repeat the first's.
with_weyl_components = {}


def time_problem(io, ℓmax, Nt, components, ℓα, budget, limit):
    P = 6 + (ℓα + 1) ** 2
    if components != ("h",) and (ℓmax, ℓα) in with_weyl_components:
        t_one, samples, kind = with_weyl_components[ℓmax, ℓα]
    else:
        # A small problem first, which compiles everything and predicts the time of the
        # full one.
        Nt_small = min(Nt, 1000)
        small = Problem(ℓmax, Nt_small, components, ℓα)
        small.objective(small.θ)
        t_small, _ = least_time(lambda: small.objective(small.θ), min(budget, 2.0))
        predicted = t_small * Nt / Nt_small
        memory = os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        if predicted > limit or grid_bytes(ℓmax, Nt, components) > 0.6 * memory:
            t_one, samples, kind = predicted, 0, "estimated"
        else:
            full = Problem(ℓmax, Nt, components, ℓα)
            t_one, samples = least_time(lambda: full.objective(full.θ), budget)
            kind = "measured"
        if components != ("h",):
            with_weyl_components[ℓmax, ℓα] = (t_one, samples, kind)
    rows = (
        ("python: one evaluation", t_one, samples, kind),
        ("python: forward", (P + 1) * t_one, 0, "estimated"),
        ("python: central", (2 * P + 1) * t_one, 0, "estimated"),
    )
    for name, seconds, n, how in rows:
        row = (ℓmax, ",".join(components), P, name, seconds, n, how)
        print("\t".join(map(str, row)), file=io, flush=True)
        note = f"least of {n}" if how == "measured" else how
        print(
            f"  ℓₘₐₓ = {ℓmax:2d} │ {' '.join(components):<18s} │ P = {P:2d} │ {name:<28s} │ "
            f"{seconds:10.4f} s │ {note}",
            flush=True,
        )


def main(args):
    file, Nt, ℓs, ℓαs, budget, limit = args
    Nt, budget, limit = int(Nt), float(budget), float(limit)
    ℓs = [int(x) for x in ℓs.split(",")]
    ℓαs = [int(x) for x in ℓαs.split(",")]
    import numba

    versions = (
        f"Python {sys.version.split()[0]}; scri {scri.__version__}; numpy {np.__version__}; "
        f"numba {numba.__version__} with {numba.get_num_threads()} threads"
    )
    print(f"Python scri, Nᵗ = {Nt}; {versions}", flush=True)
    with open(file, "w", encoding="utf-8") as io:
        print(f"#{versions}", file=io, flush=True)
        for ℓmax in ℓs:
            for components in component_sets:
                for ℓα in ℓαs:
                    if ℓα <= ℓmax:
                        time_problem(io, ℓmax, Nt, components, ℓα, budget, limit)


if __name__ == "__main__":
    main(sys.argv[1:])
