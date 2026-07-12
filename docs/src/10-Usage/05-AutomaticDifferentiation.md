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
detail below: you should **fix the output grid**, and **you must load
ForwardDiff**.

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

    At the moment, `t′` is required to have as many samples as `t`.
    This requirement will be relaxed in the future, but for now you
    might want to downsample the input `t` if you want a coarser grid.

## Load ForwardDiff

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

## Accuracy

Derivatives agree with finite differences at the `1e-8` level; see the
testitem `"transform!: ForwardDiff derivatives with respect to BMS
parameters"`.  This is most likely a statement about the accuracy of
finite differences, rather than about the forward-mode autodiff —
which is likely much more accurate.
