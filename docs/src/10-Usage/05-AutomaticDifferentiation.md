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
vector `θ`, that is `eltype(θ)`.  Here [`bms_from_parameters`](@ref)
turns `θ` into the boost velocity, rotation, and supertranslation; its
parameter layout is described [below](@ref "Gradients at a few times the
cost of the objective"):

```julia
function objective(θ)
    v⃗, R, α = bms_from_parameters(θ)
    data_dual = Complex{eltype(θ)}.(data)
    data′, _ = transform!(data_dual, t, v⃗, R, α, dc; t′=my_grid)
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
the BMS parameters `θ` that minimize the L² difference
between the transformed waveform and a fixed target.  Differentiating
`transform!` for this purpose is needlessly expensive.  With `N`
parameters, the dual copy of `data` is `(1+N)×` the size of an already
huge array; and because the distorted grid depends on `θ`, the
synthesis matrices and the entire analysis stage are dragged into dual
arithmetic as well — only for all of that structure to be collapsed to
a single scalar at the end.

[`transform_objective`](@ref) computes that scalar directly:

```math
Φ(θ) = \frac{1}{𝒩} ∑_j w_j ∑_i v_i ∑_k u_k
       \left| {}^{k}d'^{i}_{j}(θ) - {}^{k}c^{i}_{j} \right|^2,
```

with `j` running over the pixels of a grid on the sphere, `i` over the
fixed output times `t′`, and `k` over the data components, and with
``c`` the target.  It accumulates the ``L²`` sum as the pixels are
transformed, a block of pixels at a time, so the transformed waveform
is never stored in full, `data` is never modified (it is reused across
optimizer iterations), and dual numbers are confined to the geometry
of the pixels and to storage for one block of pixels.  The analysis
stage (converting back to waveform modes) is skipped entirely.

Everything in ``Φ`` that does not depend on ``θ`` is stored in an
[`ObjectivePlan`](@ref), which is built once: the grid and its weights
``w_j``, the output times and their weights ``v_i``, the component
weights ``u_k``, the target's values ``{}^{k}c^{i}_{j}`` on the grid,
and the normalization ``𝒩``.  The grid is an exact quadrature on the
sphere for band-limited fields, and the time weights are for the
trapezoidal rule on `t′`, so ``Φ`` approximates the integral of
``|d'-c|^2`` over the sphere and the time window — an ``L²`` norm.  A
boosted field is not band-limited, so the sphere part is spectrally
accurate rather than exact; see [`ObjectivePlan`](@ref) for the
precise condition and for the option to pass a finer grid.

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

# 3. The differentiable objective, with θ read by `bms_from_parameters`
#    as a rotation vector, a rapidity vector, and the supertranslation's
#    modes; see below.
objective(θ) = transform_objective(data, t, bms_from_parameters(θ)..., dc, plan)

# Compute the gradient at a starting point.  The identity
# transformation is θ = 0; here with a supertranslation up to ℓ = 2.
θ₀ = zeros(6 + 9)
g = ForwardDiff.gradient(objective, θ₀)
```

Any optimizer that takes a function and its gradient will do, such as
a quasi-Newton method like L-BFGS, with the gradient from ForwardDiff
as above.

A few things to know:

- **Parameters.**  The parameter vector ``θ`` is laid out as for
  [`objective_and_gradient!`](@ref), described below: a rotation
  vector ``\vec{r}``, a rapidity vector ``\vec{η}``, and the
  supertranslation's modes.  So the same starting point and the same
  result serve both ways of computing the gradient.
  `transform_objective` also accepts the transformation directly;
  there, three rotation components are read as a rotation vector, and
  four as a rotor `(w, x, y, z)`, which is normalized; see
  [`rotor_argument`](@ref Scri.rotor_argument).

- **Weights and noise.**  The plan's keywords set the weights:
  `time_weights` for the ``v_i``, `component_weights` for the ``u_k``,
  and `weights` for per-pixel factors that multiply the quadrature
  weights, for example to down-weight part of the sky.
  `normalize=true` divides by the target's norm.  `ℓ_fit=L` compares
  only the modes with ``ℓ ≤ L`` of the difference, projected exactly
  by the quadrature, which leaves out high-``ℓ`` noise in the data or
  the target.

- **Memory.**  The dual-typed storage is for one block of pixels at a
  time: their harmonics and their synthesized time series, a copy of
  each in plain numbers for the matrix products, and a few buffers of
  size `O(NᵗNᵈ)` for each thread.  The block is sized so that its
  harmonics and series take no more memory than `data` itself, or 128
  MiB if that is larger, but it holds at least four pixels for each
  thread; with the copies, the total is at most about twice that.
  Differentiating through `transform!` instead takes `(NᵐNᵗNᵈ + 5Nᵐ²)
  × (1+N) × 16` bytes, with `N` parameters, for dual `data` plus dual
  synthesis matrices, which is more than `1+N` times the size of
  `data`.

- **Speed.**  The synthesis of each block of pixels is one matrix
  product for each component.  `Complex{Dual}` is not a BLAS type, so
  the values and partials of the dual harmonics are first laid side by
  side in one matrix of plain complex numbers; the product, which is
  linear in them, is then a single BLAS call.  A dual `transform!`
  runs on generic matrix multiplication instead, and its analysis
  stage is eliminated here entirely.  The default grid has about twice
  as many points as there are modes, so an evaluation costs about
  twice what a sum over `transform!`'s own pixels would, but this is
  essential to ensure that the *squared* difference is integrated
  correctly, without aliasing.

- **Stability.**  The spline knots and Thomas factors stay primal
  (duals enter only because the values being interpolated are dual, so
  the conditioning of the interpolation is unchanged).  Each pixel's
  share of the objective is computed on its own, and the shares are
  summed in pixel order.  The result is deterministic for a fixed
  number of threads.

## Gradients at a few times the cost of the objective

ForwardDiff's gradient costs about one evaluation of the objective per
parameter: with a supertranslation up to ``ℓ = 4``, that is 31
evaluations.  [`objective_and_gradient!`](@ref) computes the objective
and its exact gradient together, at a few times the cost of the
objective alone, whatever the number of parameters: about four
evaluations at ``ℓ_\mathrm{max} = 16``, or five with `ℓ_fit`.  It
takes the transformation as one parameter vector,

```math
θ = (\vec{r},\; \vec{η},\; a_1, …, a_{(ℓ_α+1)^2}),
```

a rotation vector ``\vec{r}`` (the frame rotation is
``\exp(\vec{r}/2)``), a rapidity vector ``\vec{η}`` (the boost
velocity is ``\tanh|\vec{η}|`` along ``\vec{η}``), and the real
degrees of freedom of the supertranslation's modes up to ``ℓ_α``.
Every ``θ`` describes a valid BMS transformation, and the identity is
``θ = 0``.  The plan's fixed output grid `t′` is valid only for the
transformations it was built for, though: those with a boost speed up
to `β_max` and a supertranslation up to `δt_max`, as [above](@ref "Fix
the output time grid").  Past those bounds, the call throws an
`ArgumentError`, which stops most optimizers.  So either build `t′`
with a generous margin, or bound the optimizer, for example by keeping
the rapidity below ``\operatorname{atanh}(β_\mathrm{max})``.

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
