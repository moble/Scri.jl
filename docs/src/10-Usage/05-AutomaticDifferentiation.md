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

Before differentiating, you *should* **fix the output time grid**, as
described below.

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
    The optimization-oriented [`transform_objective`](@ref) and its
    [`ObjectivePlan`](@ref), described below, have no such
    restriction: the plan's `t′` may have any length, so a coarse (and
    correspondingly cheap) output grid can be chosen instead.

## Using `ForwardDiff`

`ForwardDiff` is a dependency of `Scri`, which keeps dual numbers out
of the parameter-independent pixel grid and its matrix factorizations.
These quantities depend only on the band limit, not on the BMS
parameters, so their derivatives are identically zero; computing them
with dual numbers would be not only wasteful but numerically fragile,
and derivatives with respect to the rotation would come out as `NaN`.

[`objective_and_gradient!`](@ref), described below, needs nothing more
than `using Scri`.  To call `ForwardDiff` yourself, as in the examples
that differentiate [`transform!`](@ref) and
[`transform_objective`](@ref), add it to your own project, as for any
package you use directly:

```julia
using Pkg
Pkg.add("ForwardDiff")
```

and load it in the session:

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
huge array; and because the distorted grid depends on `θ`, the
synthesis matrices and the entire analysis stage are dragged into dual
arithmetic as well — only for all of that structure to be collapsed to
a single scalar at the end.

[`transform_objective`](@ref) computes that scalar directly:

```math
J(θ) = \frac{1}{𝒩} ∑_i w_i ∑_j v_j ∑_k c_k
       \left| d'_{ijk}(θ) - d_{\mathrm{t};ijk} \right|^2,
```

with `i` running over the pixels of a grid on the sphere, `j` over the
fixed output times `t′`, and `k` over the data components.  It fuses
the ``L²`` accumulation into the per-pixel transformation loop, so the
transformed waveform is never materialized, `data` is never modified
(it is reused across optimizer iterations), and dual numbers are
confined to per-pixel scalars and small task-local buffers.  The
analysis stage (converting back to waveform modes) is skipped
entirely.

Everything in ``J`` that does not depend on ``θ`` is stored in an
[`ObjectivePlan`](@ref), which is built once: the grid and its weights
``w_i``, the output times and their weights ``v_j``, the component
weights ``c_k``, the target's values ``d_{\mathrm{t};ijk}`` on the
grid, and the normalization ``𝒩``.  The grid is an exact quadrature
on the sphere for band-limited fields, and the time weights are the
trapezoid rule on `t′`, so ``J`` approximates the integral of ``|d' -
d_\mathrm{t}|^2`` over the sphere and the time window — a true ``L²``
norm.

The recipe has three steps — the first two run *once*, outside the
optimization loop:

```julia
using Scri, ForwardDiff

# 1. A fixed output grid, valid for every parameter value the
#    optimizer may try (β_max and δt_max bound the boost speed and
#    supertranslation magnitude).  Any length works: a coarse grid
#    makes each iteration correspondingly cheaper.
t′ = compute_t′(t, β_max, δt_max)

# 2. Everything that does not depend on the parameters, from the
#    target's modes and its own time samples.  (If the target has a
#    smaller band limit than `data`, zero-pad its modes up to the same
#    Nᵐ first.)
plan = ObjectivePlan(target_modes, t_target, dc, t′)

# 3. The differentiable objective.  Note that `v⃗` and `𝔯` accept
#    plain vectors, so slices of the parameter vector pass straight
#    through.  Three rotation components are a rotation vector; see
#    below.
objective(θ) = transform_objective(data, t, θ[1:3], θ[4:6], α(θ[7:end]), dc, plan)

# Compute the gradient at the starting point θ₀.  The identity
# transformation is θ = 0.
g = ForwardDiff.gradient(objective, θ₀)
```

Any optimizer that takes a function and its gradient will do, such as
a quasi-Newton method like L-BFGS, with the gradient from ForwardDiff
as above.

A few things to know:

- **Rotation parameters.**  Three components are read as a rotation
  vector ``𝐮``, the rotation by the angle ``|𝐮|`` about the axis
  ``𝐮``, which is the rotor ``\exp(𝐮/2)``; so the identity is zero,
  and every parameter is used.  Four components are read as a rotor
  `(w, x, y, z)` and normalized, which also works but leaves one
  redundant parameter; see [`rotor_argument`](@ref Scri.rotor_argument).

- **Weights and noise.**  The plan's keywords set the weights:
  `time_weights` for the ``v_j``, `component_weights` for the
  ``c_k``, and `weights` for per-pixel factors that multiply the
  quadrature weights, for example to down-weight part of the sky.
  `normalize=true` divides by the target's norm.  `ℓ_fit=L` compares
  only the modes with ``ℓ ≤ L`` of the difference, projected exactly by
  the quadrature, which leaves out high-``ℓ`` noise in the data or the
  target.

- **Memory.**  Per thread, the dual-typed storage is a few buffers of
  size `O(Nᵈ Nᵗ)` — with `ntasks` threads and `N` parameters, roughly
  `ntasks × (3NᵗNᵈ + 2Nᵐ) × (1+N) × 16` bytes — compared to `(NᵐNᵗNᵈ +
  5Nᵐ²) × (1+N) × 16` bytes for dual `data` plus dual synthesis
  matrices when differentiating through `transform!`.  For
  production-sized arrays this can be a reduction by an order of
  magnitude or two.

- **Speed.**  Nothing important is lost by leaving the GEMM
  formulation behind: `Complex{Dual}` is not a BLAS type, so a dual
  `transform!` already runs on generic matrix multiplication.  The
  per-pixel dot products do the same work, and the analysis stage is
  eliminated entirely.  The default grid has about twice as many
  points as there are modes, so an evaluation costs about twice what a
  sum over `transform!`'s own pixels would, but this is essential to
  ensure that the *squared* difference is integrated correctly,
  without aliasing.

- **Stability.**  The spline knots and Thomas factors stay primal
  (duals enter only because the values being interpolated are dual, so
  the conditioning of the interpolation is unchanged), and the results
  accumulate per-pixel → per-task → across tasks.  The result is
  deterministic for a fixed number of threads.

## Gradients at about twice the cost of the objective

ForwardDiff's gradient costs about one evaluation of the objective per
parameter: with a supertranslation up to ``ℓ = 4``, that is 31
evaluations.  [`objective_and_gradient!`](@ref) computes the objective
and its exact gradient together, at about twice the cost of the
objective alone, whatever the number of parameters: 2.5 evaluations at
``ℓ_\mathrm{max} = 8``, and 1.7 at ``ℓ_\mathrm{max} = 16``.  It takes
the transformation as one parameter vector,

```math
θ = (𝐮,\; 𝐛,\; a_1, …, a_{(ℓ_α+1)^2}),
```

a rotation vector ``𝐮`` (the frame rotation is ``\exp(𝐮/2)``), a
rapidity vector ``𝐛`` (the boost velocity is ``\tanh|𝐛|`` along
``𝐛``), and the real degrees of freedom of the supertranslation's
modes up to ``ℓ_α``.  Every ``θ`` is a valid transformation, so an
optimizer needs no bounds, and the identity is ``θ = 0``.

```julia
using Scri

t′ = compute_t′(t, β_max, δt_max)
plan = ObjectivePlan(target_modes, t_target, dc, t′)

# The objective and its gradient, in the form optimizers usually take.
fg!(∇, θ) = objective_and_gradient!(∇, θ, data, t, dc, plan)

θ₀ = zeros(6 + 9)  # the identity, with a supertranslation up to ℓ = 2
# ... minimize with any method that uses the value and the gradient ...
v⃗, R, α = bms_from_parameters(θ_best)
```

[`bms_from_parameters`](@ref) turns the result back into the boost
velocity, rotation, and supertranslation that [`transform!`](@ref)
takes, and [`parameters_from_bms`](@ref) gives the ``θ`` of a known
transformation, to start from it.  Plans with `ℓ_fit` work as well, at
the cost of one more pass over the pixels to project the transformed
waveform before the gradient is computed.

The gradient is computed by running the expensive part of each pixel's
transformation backwards: the residual is propagated through the
mixing of the components, the spline interpolation in time, and the
synthesis from the modes, which costs one more pass as expensive as the
synthesis.  The geometry of each pixel, which depends on the
parameters but costs only a few operations per pixel, is differentiated
with dual numbers.  The derivation is on the page [The gradient of the
objective](@ref objective_gradient).

## Accuracy

Derivatives agree with finite differences at the `1e-8` level; see the
testitem `"transform!: ForwardDiff derivatives with respect to BMS
parameters"`.  This is most likely a statement about the accuracy of
the finite differences, rather than about the forward-mode autodiff —
which is likely far more accurate.
