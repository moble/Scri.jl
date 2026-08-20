```@meta
CurrentModule = Scri
```

# Automatic Differentiation

[`transform!`](@ref) can be differentiated end-to-end with
[ForwardDiff](https://juliadiff.org/ForwardDiff.jl/) with respect to
the BMS parameters — boost velocity, rotation, and supertranslation
(hence also time and space translations).  This makes it possible to
drop a transformation directly into an optimization loop — for
example, minimizing some functional of the transformed data over the
BMS group — and let the autodiff machinery supply exact gradients.

There are two things to know before differentiating, described in
detail below:

  1. you *should* **fix the output time grid**, and
  2. you *must* **load `ForwardDiff`**.

## Fix the output time grid

Boosts and supertranslations mix time and direction, so the default
output grid `t′` is itself a function of the BMS parameters (see [The
output time grid](@ref)).  If you let `transform!` construct that grid
for you, the derivative of the data at output index `j` will mix the
derivative of the *field* with the motion of the *grid* — the sample
point at index `j` slides along the worldline as the parameters
change, and the returned derivative reflects both effects at once.
That is almost never what you want.  Moreover, that would also most
likely require re-interpolating the output data to a fixed time grid
anyway, which would be wasteful.

Instead, pass the `t′` keyword to pin the output grid in place.  With
a fixed grid the sample points stay put, so you get clean partial
derivatives of the field at fixed output times:

```julia
transform!(data, t, v⃗, R, α, dc; t′=my_grid)
```

The cost of fixing the grid is that you must choose its range in
advance, and it must remain within the valid span for *every*
parameter value the optimizer will try — not just the starting point.
Choose the range small enough that you are confident it will always be
covered.

The [`compute_t′(t, β, δt)`](@ref compute_t′(::Any, ::Real, ::Real))
method does this for you.  Given the input time array `t` and upper
bounds on the boost speed `β` and supertranslation magnitude `δt` you
expect to encounter, it returns a time array guaranteed to be valid
for every transformation within those bounds:

```julia
t′ = compute_t′(t, β_max, δt_max)
```

This does require some knowledge of the parameter ranges in advance,
but this might be tied into your optimization bounds.  Note that both
parameters are magnitudes rather than directed quantities.  In
particular, `δt_max` would be the largest absolute value of the
supertranslation function over the sphere with a worst-case
supertranslation.

The output will contain as many samples as the input `t`, being just a
scaled-and-shifted version of it.

!!! note

    For `transform!` itself, `t′` is required to have as many samples
    as `t`, because the result is written back into `data` in place.
    The optimization-oriented [`transform_objective`](@ref) and
    [`pixel_waveform`](@ref) described below have no such restriction:
    their `t′` may have any length, so a coarse (and correspondingly
    cheap) output grid is a one-line choice there.

## Load `ForwardDiff`

You *must* load `ForwardDiff` in your session.  Loading it activates a
package extension that keeps dual numbers out of the
parameter-independent pixel grid and its matrix factorizations.  These
quantities depend only on the band limit, not on the BMS parameters,
so their derivatives are identically zero; computing them with dual
numbers is not only wasteful but numerically fragile.  Without the
extension, derivatives with respect to the rotation come out as `NaN`.

This means that you must add `ForwardDiff` to the same project you're
using for `Scri`, and you must run `using ForwardDiff` in the same
session as you want to use it.  That is, do something like this once
in whichever project you will want to use autodiff with `Scri`:

```julia
using Pkg
Pkg.add("ForwardDiff")
```

And do something like this in every session where you want to actually
do that differentiation:

```julia
using ForwardDiff
```

## Preparing the data array

!!! tip

    This section applies when you need the full transformed *data*
    with derivatives.  If you only need a scalar objective — as in
    optimization of the BMS parameters — skip ahead to [Optimizing BMS
    parameters without the memory cost](@ref): with
    [`transform_objective`](@ref) the `data` array stays primal and
    none of the preparation below is necessary.

Because `transform!` works in place, the `data` array must already
have the dual element type before you call it.  One allocation-heavy
approach is inside the differentiated function to convert the array to
the dual type of the parameter being differentiated — for a parameter
`θ`, that is `typeof(θ)`:

```julia
function objective(θ)
    data_dual = Complex{typeof(θ)}.(data)
    data′, _ = transform!(data_dual, t, v⃗(θ), R(θ), α(θ), dc; t′=my_grid)
    # ... reduce data′ to a scalar ...
end

g = ForwardDiff.gradient(objective, θ₀)
```

This would re-allocate on every iteration, which could be wasteful.  A
more efficient approach is to allocate a dual-typed array once and
reuse it.  This is a little more complicated, but can be done nicely
using
[`PreallocationTools.jl`](https://docs.sciml.ai/PreallocationTools/stable/).

## Optimizing BMS parameters without the memory cost

The most common use of differentiation here is optimization: finding
the BMS parameters `θ = (v⃗, 𝔯, α)` that minimize the L² difference
between the transformed waveform and a fixed target.  Differentiating
`transform!` for this purpose is needlessly expensive.  With `N`
parameters, the dual copy of `data` is `(1+N)×` the size of an already
huge array; and because the distorted grid depends on `θ`, the five
synthesis matrices and the entire analysis stage are dragged into dual
arithmetic as well — only for all of that structure to be collapsed to
a single scalar at the end.

[`transform_objective`](@ref) computes that scalar directly:

```math
\mathrm{obj}(θ) = ∑_i w_i ∑_j ∑_k \left| d'_{ijk}(θ) -
\mathrm{target}_{ijk} \right|^2,
```

with `i` running over pixels, `j` over the fixed output times `t′`,
and `k` over the data components.  It fuses the ``L²`` accumulation
into the per-pixel transformation loop, so the transformed waveform is
never materialized, `data` is never modified (it is reused across
optimizer iterations), and dual numbers are confined to per-pixel
scalars and small task-local buffers.  The analysis stage (converting
back to waveform modes) is skipped entirely.

An optimization routine typically needs the objective and possibly its
gradient.

The recipe has three steps — the first two run *once*, outside the
optimization loop:

```julia
using Scri, ForwardDiff

# 1. A fixed output grid, valid for every parameter value the
#    optimizer may try (β_max and δt_max bound the boost speed and
#    supertranslation magnitude).  Any length works: a coarse grid
#    makes each iteration correspondingly cheaper.
t′, _ = compute_t′(t, β_max, δt_max)

# 2. Compute the target waveform's pixel values on that grid.  This is
#    θ-independent, so we build it once and reuse it on every iteration.
#    (If the target has a smaller band limit than `data`, zero-pad its
#    modes up to the same Nᵐ first — the pad-first workflow.)
target = pixel_waveform(h_target_modes, t_target, dc, t′)

# 3. The differentiable objective.  Note that `v⃗` and `𝔯` accept
#    plain vectors, so slices of the parameter vector pass straight
#    through.
objective(θ) = transform_objective(
    data, t, θ[1:3], θ[4:6], α(θ[7:end]), dc, target; t′=t′
)

# Compute the gradient at the starting point θ₀.
g = ForwardDiff.gradient(objective, θ₀)
```

A few things to know:

  - **Norm equivalence.**  The objective is a *pixel-domain* ``L²``
    norm.  The pixel grid has exactly as many points as modes and the
    synthesis matrices are full rank, so it vanishes exactly when the
    mode-domain difference vanishes, and defines an equivalent norm —
    but its numerical value differs from `sum(abs2, ...)` over mode
    weights.  Objective values are comparable across iterations, which
    is all an optimizer needs.

  - **Memory.**  Per thread, the dual-typed storage is a few buffers
    of size `O(Nᵈ Nᵗ)` — with `ntasks` threads and `N` parameters,
    roughly `ntasks × (3NᵗNᵈ + 2Nᵐ) × (1+N) × 16` bytes — compared to
    `(NᵐNᵗNᵈ + 5Nᵐ²) × (1+N) × 16` bytes for dual `data` plus dual
    synthesis matrices when differentiating through `transform!`.  For
    production-sized arrays this can be a reduction by an order of
    magnitude or two.

  - **Speed.**  Nothing important is lost by leaving the GEMM
    formulation behind: `Complex{Dual}` is not a BLAS type, so a dual
    `transform!` already runs on generic matrix multiplication.  The
    per-pixel dot products do the same work, and the analysis stage is
    eliminated entirely.

  - **Stability.**  The spline knots and Thomas factors stay primal
    (duals enter only because the values being interpolated are dual,
    so the conditioning of the interpolation is unchanged), and the
    results accumulate per-pixel → per-task → across tasks.  The
    result is deterministic for a fixed number of threads.

## Accuracy

Derivatives agree with finite differences at the `1e-8` level; see the
testitem `"transform!: ForwardDiff derivatives with respect to BMS
parameters"`.  This is most likely a statement about the accuracy of
the finite differences, rather than about the forward-mode autodiff —
which is likely far more accurate.
