```@meta
CurrentModule = Scri
```

# [The gradient of the objective](@id objective_gradient)

One of the biggest challenges when comparing waveforms — whether from
numerical simulations, or from surrogate, analytical (e.g., PN or
EMRI), or quasi-analytical (e.g., EOB) models — is understanding how
the gauge freedom affects the comparison.  *All* waveforms have this
intrinsic BMS freedom.  One way to deal with this freedom is to
optimize the BMS transformation that brings one waveform into the best
agreement with another.  The obvious measure to use for "best
agreement" is the weighted ``L^2`` difference between the two
waveforms — integrated over both the sphere and time.  We call this
the objective function ``Φ(θ)``, where ``θ`` parameterize the BMS
transformation, implemented as the function
[`transform_objective`](@ref).  One crucial element in making this
efficient is the fact that the objective function eliminates the need
to compute the full transformed waveform; we don't need to transform
back to modes if we just need to *integrate* the difference from the
other waveform.  But we can do better if we think specifically about
the goal here.

To actually minimize this quantity efficiently, we don't just need its
*value*, we also need its *gradient* with respect to the parameters of
the BMS transformation, ``∂Φ/∂θ^p``.  There are six parameters for the
Lorentz transformation and ``(ℓ_α+1)^2`` more for a supertranslation
with modes up to ``ℓ_α``, so there are easily dozens of parameters,
and the obvious ways of computing the gradient — finite differences,
or dual numbers in every parameter — cost about one evaluation of
``Φ`` for each parameter.

This page explains how [`objective_and_gradient!`](@ref) computes the
value *and* gradient of the objective function with respect to all of
the parameters efficiently, at the cost of only a few evaluations of
the objective function itself — almost independently of the number of
parameters.  In brief, ``Φ`` is computed by a chain of four simple
steps, repeated for each direction on the sphere.  The parameters
enter that chain only through a small amount of geometry for each
direction, which is cheap to compute but complicated.  The data enter
through sums over every mode and every time, which are expensive but
simple.  We differentiate the geometry *forward*, with dual numbers,
and the four steps *backward*, by hand, and the two meet in the
middle.

## The objective

The parameters are

```math
θ = \left(\vec{r},\ \vec{η},\ a_1, …, a_{(ℓ_α+1)^2}\right):
```

a rotation vector ``\vec{r}``, a rapidity vector ``\vec{η}``, and the
real degrees of freedom of the supertranslation's modes.
[`bms_from_parameters`](@ref) turns them into the boost velocity, the
frame rotation, and the supertranslation that [`transform!`](@ref)
takes.  We write ``P`` for the number of parameters.

The objective is

```math
Φ(θ) = \frac{1}{𝒩} \sum_j w_j \sum_i v_i \sum_k u_k
       \left| μ^{ki}_{j}(θ) - c^{ki}_{j} \right|^2,
```

where

- ``j`` labels the *pixels*, which are the points of a grid of
  directions on the sphere of the new frame, and the ``w_j`` are
  weights for integrating over that sphere;
- ``i`` labels the output times ``t'_i``, and the ``v_i`` are weights
  for integrating over them;
- ``k`` labels the components of the data, such as the strain and
  ``ψ_4``, and the ``u_k`` weight them against each other;
- ``μ^{ki}_{j}(θ)`` is the transformed waveform: its component ``k``
  at pixel ``j`` and time ``t'_i``;
- ``c^{ki}_{j}`` is the waveform that it is compared with, at the
  same pixels and times; and
- ``𝒩`` is a normalization, either 1 or the same sum with ``μ`` set
  to zero.

So ``Φ`` is a weighted and squared ``L^2`` distance between the two
waveforms, over the sphere and the output times.  Everything in it
except ``μ`` is fixed before the minimization starts, and is stored
in an [`ObjectivePlan`](@ref), whose docstring describes the grid, the
weights, and the choices that the caller has.

Each pixel contributes to ``Φ`` on its own:

```math
Φ = \frac{1}{𝒩} \sum_j ν_j,
\qquad
ν_j = \sum_{i,k} W_{jki} \left| μ^{ki}_{j} - c^{ki}_{j} \right|^2,
\qquad
W_{jki} = w_j\, v_i\, u_k.
```

The computation of ``ν_j`` involves no pixel but ``j``.  So the pixel
index passes untouched through everything that follows, until the sum
over pixels at the very end, and we can differentiate one pixel at a
time and add the results.

### Indices

| Index | Labels | Number |
| :--- | :--- | :--- |
| ``p`` | the parameters ``θ^p`` | ``P`` |
| ``j`` | the pixels | ``N_\mathrm{pix}`` |
| ``ℓm`` | the modes; always a pair, treated as one index | ``N_\mathrm{modes}`` |
| ``h`` | the input times ``t_h`` | ``N_t`` |
| ``i`` | the output times ``t'_i`` | ``N_{t'}`` |
| ``k``, ``k'`` | the components of the data; ``s_k`` is the spin weight of component ``k`` | ``N_\mathrm{comp}`` |
| ``n`` | eight real numbers that describe the geometry of a pixel | ``8`` |

Indices are placed as in general relativity.  An upper index is a
slot of a vector, and a lower index is a slot of a covector; an upper
index in the denominator of a derivative counts as a lower index.  An
index that appears once up and once down in a product is summed, and
every other sum is written out.

The pixel index ``j`` is only a label.  It is always written as the
first lower index of whatever belongs to a pixel, and it is never
summed except by the ``\sum_j`` above.  The weights ``W_{jki}`` and
the times ``t_h``, ``t'_i``, and ``τ_{ji}`` are fixed arrays of
numbers, and all of their indices are labels of the same kind.

### The four stages

The transformed waveform at pixel ``j`` is computed from the modes of
the original waveform in three stages, and a fourth compares the
result with ``c``.  We name the stages ``σ``, ``ι``, ``μ``, and ``ν``
— for synthesis, interpolation, mixing, and norm — and we use the same
letter, with indices, for the output of each stage.  (On this page,
``σ`` never denotes the shear.)  The first three answer, in turn,
*where*, *when*, and *in which basis* the original waveform is to be
evaluated.

The input is the array ``d^{kh}{}_{ℓm}``: mode ``(ℓ,m)`` of component
``k`` of the original waveform at the input time ``t_h``.

**Synthesis, ``σ``.**  Because of the boost and the rotation, a
direction in the new frame corresponds to a different direction in the
original frame.  The first stage evaluates the original waveform in
that direction, at each of its own time samples:

```math
σ^{kh}_{j} = d^{kh}{}_{ℓm}\ {}_{s_k}Y^{ℓm}_{j}.
```

Here ``{}_{s}Y^{ℓm}_{j}`` stands for the spin-weighted spherical
harmonic evaluated at the rotor of pixel ``j`` in the original frame,
which [`aberration`](@ref) computes, and the sum runs over the modes
with ``ℓ ≥ |s_k|``.  The harmonics are functions of the rotor, and not
only of the direction, so this stage also accounts for the rotation of
the dyad.

We usually think of such a sum as a linear map applied to the modes,
with the harmonics as its coefficients.  Here the roles are reversed.
The modes are fixed data, and the harmonics are what depends on ``θ``.
So we regard ``σ`` as a linear map applied to the *harmonics*, with
the *data* as its coefficients, and the indices are placed
accordingly.

**Interpolation, ``ι``.**  At pixel ``j``, the output time ``t'_i``
corresponds to the original time

```math
τ_{ji} = κ_j^{-1}\, t'_i + α_j,
```

where ``κ_j`` is the conformal factor and ``α_j`` is the
supertranslation, both evaluated at that pixel.  This is the law
``t' = κ(t - α)`` solved for ``t``, with any sign convention for ``α``
already applied.  The second stage interpolates each component to
those times, with the natural cubic spline through its samples:

```math
ι^{ki}_{j} = S_{j}{}^{i}{}_{h}\ σ^{kh}_{j}.
```

The value of a spline at any time is a linear combination of the
samples it passes through, and ``S_{j}{}^{i}{}_{h}`` is the
coefficient of the sample at ``t_h`` in the value at ``τ_{ji}``.
These coefficients are real.  They depend on the times, but not on
the samples.

**Mixing, ``μ``.**  The components are defined relative to a tetrad,
and the tetrad of the new frame differs from the original one by a
rescaling and a null rotation.  So each new component is a combination
of the original ones, and the shear and the strain are also shifted by
a multiple of ``ð^2α``:

```math
μ^{ki}_{j} = M_{j}{}^{k}{}_{k'}(τ_{ji})\ ι^{k'i}_{j} + q^{k}_{j}.
```

These are the transformation laws of the page [BMS Action on
Fields](@ref "BMS Action on Fields"), which [`mix_components!`](@ref)
applies.  For example, if the components are the strain, ``ψ_3``, and
``ψ_4`` at ``ℐ^+``, in that order and in the default conventions, then

```math
\begin{aligned}
μ^{1i}_{j} &= κ_j^{-1} \left( ι^{1i}_{j} - \bar{ð}^2 α_j \right), \\
μ^{2i}_{j} &= κ_j^{-3} \left( ι^{2i}_{j} + b_j(τ_{ji})\, ι^{3i}_{j} \right), \\
μ^{3i}_{j} &= κ_j^{-3}\, ι^{3i}_{j},
\end{aligned}
```

from which ``M_{j}{}^{k}{}_{k'}`` and ``q^{k}_{j}`` can be read off.
The null-rotation parameter ``b = ðt'/2κ`` is linear in the original
time,

```math
b_j(τ) = b^{(0)}_{j} + b^{(1)}_{j}\, τ,
```

which is why ``M_{j}{}^{k}{}_{k'}`` depends on ``τ_{ji}``; see
[Computing ``ðt'/2κ``](@ref computing_eth_tprime_over_2kappa).

**Norm, ``ν``.**  The last stage compares the result with ``c`` and
reduces everything at the pixel to one real number:

```math
ν_j = \sum_{i,k} W_{jki} \left| μ^{ki}_{j} - c^{ki}_{j} \right|^2.
```

Each stage acts on one kind of index and passes the others along
untouched, and the pixel index passes through all four:

| Stage | Output | Acts on | Passes along |
| :--- | :--- | :--- | :--- |
| synthesis | ``σ^{kh}_{j} = d^{kh}{}_{ℓm}\ {}_{s_k}Y^{ℓm}_{j}`` | the modes ``ℓm`` | ``k``, ``h``, ``j`` |
| interpolation | ``ι^{ki}_{j} = S_{j}{}^{i}{}_{h}\ σ^{kh}_{j}`` | the times, ``h → i`` | ``k``, ``j`` |
| mixing | ``{μ^{ki}_{j} = M_{j}{}^{k}{}_{k'}(τ_{ji})\ ι^{k'i}_{j} + q^{k}_{j}}`` | the components, ``k' → k`` | ``i``, ``j`` |
| norm | ``{ν_j = \sum_{i,k} W_{jki} \lvert μ^{ki}_{j} - c^{ki}_{j} \rvert^2}`` | the times ``i`` and the components ``k`` | ``j`` |

Only the final sum, ``Φ = 𝒩^{-1} \sum_j ν_j``, acts on the pixels.

### Where the parameters enter

The data, the times, the weights, and ``c`` do not depend on ``θ``.
The parameters enter only through the *geometry* of each pixel, which
consists of two sets of numbers:

- the harmonics ``{}_{s}Y^{ℓm}_{j}``, for each spin weight ``s`` among
  the components, which enter ``σ``; and
- eight real numbers,

  ```math
  \left(g^1_j, …, g^8_j\right)
  = \left( κ_j^{-1},\ α_j,\
    \operatorname{Re} b^{(0)}_{j},\ \operatorname{Im} b^{(0)}_{j},\
    \operatorname{Re} b^{(1)}_{j},\ \operatorname{Im} b^{(1)}_{j},\
    \operatorname{Re} ð^2α_j,\ \operatorname{Im} ð^2α_j \right),
  ```

  which enter ``ι`` through the times ``τ_{ji}``, and ``μ`` through
  ``M_{j}{}^{k}{}_{k'}`` and ``q^{k}_{j}``.

So, for each pixel,

```math
ν_j = ν\Bigl(μ\bigl(ι\left(σ(Y_j),\, g_j\right),\, g_j\bigr)\Bigr),
```

with ``Y_j`` and ``g_j`` functions of ``θ``.

### Where the work is

A waveform has hundreds of modes and thousands of time samples, or
more of each.  The synthesis is a sum over every mode, for every time
sample of every component, which is ``N_\mathrm{comp}\, N_t\,
N_\mathrm{modes}`` multiplications for each pixel.  For all but the
smallest band limits, that is much more work than everything else
together.  The other stages never touch the modes: the interpolation
costs a few operations for each input and output time, and the mixing
and the norm cost a few for each output time.  The geometry never
touches the times: it costs a few operations for each mode.  So the
picture to keep in mind is of one expensive sum, ``σ``, with cheap
computations on either side of it.

## The chain rule, multiplied in two orders

Everything on this page follows from the chain rule, together with
one observation about the order in which its factors are multiplied.
Suppose that a real number ``ν`` is computed from parameters ``θ^p``
in stages: first some quantities ``a^A`` are computed from the
``θ^p``, then some other quantities ``b^B`` from the ``a^A``, and
then ``ν`` from the ``b^B``.  The chain rule says that

```math
\frac{∂ν}{∂θ^p}
= \frac{∂ν}{∂b^B}\, \frac{∂b^B}{∂a^A}\, \frac{∂a^A}{∂θ^p}.
```

The right-hand side is a product of three arrays — the Jacobians of
the three stages — with sums over ``A`` and ``B``.  The product may be
multiplied out starting from either end.  The answer is the same, but
the work is not, and the two choices are the two "modes" of automatic
differentiation.

**Forward mode** starts at the right, with the parameters.  It
computes ``∂a^A/∂θ^p`` along with ``a^A``, then

```math
\frac{∂b^B}{∂θ^p} = \frac{∂b^B}{∂a^A}\, \frac{∂a^A}{∂θ^p}
```

along with ``b^B``, and so on to the end.  At every stage it holds the
derivatives of the current quantities with respect to all ``P``
parameters.  In geometric terms, it pushes the ``P`` coordinate
vectors ``∂/∂θ^p`` forward through each map.  This is what dual
numbers do: every number is stored together with its ``P``
derivatives, and every arithmetic operation updates both.  Nothing has
to be remembered from one stage to the next, but all the work of every
stage is done ``P + 1`` times.

**Reverse mode** starts at the left, with the result.  It computes the
covector ``∂ν/∂b^B``, then

```math
\frac{∂ν}{∂a^A} = \frac{∂ν}{∂b^B}\, \frac{∂b^B}{∂a^A},
```

and so on back to the beginning.  At every stage it holds the
derivatives of the *one* number ``ν`` with respect to the current
quantities.  In geometric terms, it pulls the gradient ``dν`` back
through each map.  There is one covector, however many parameters
there are, so the work of every stage is done roughly twice: once to
evaluate it, and once to pull the covector back through it.  The
price is paid in bookkeeping.  The Jacobians are needed in the
opposite order from the one in which the evaluation produces them, so
the whole evaluation must be finished, and enough of it remembered,
before the backward pass can start.  A general-purpose tool has to
record every operation as it goes.

| | Forward mode | Reverse mode |
| :--- | :--- | :--- |
| Starts from | the parameters | the result |
| Holds, at each stage | ``P`` vectors: the derivatives of the current quantities with respect to each ``θ^p`` | one covector: the derivatives of ``ν`` with respect to the current quantities |
| Work | ``P + 1`` evaluations | a few evaluations, whatever ``P`` is |
| Must remember | nothing | the evaluation |

Forward mode suits a computation with few inputs, then, and reverse
mode suits one with few outputs.  But nothing forces us to use one
mode throughout.  The product can be cut at any intermediate quantity,
say ``a^A``: we push the ``P`` vectors forward as far as
``∂a^A/∂θ^p``, pull the covector back as far as ``∂ν/∂a^A``, and
contract the two there.

!!! note "Other names for these things"

    In the literature on automatic differentiation, ``∂a^A/∂θ^p`` is
    called a *tangent*, and the step that advances it is called a
    *pushforward* or a *Jacobian–vector product*.  The covector
    ``∂ν/∂a^A`` is called a *cotangent* or an *adjoint*, and the step
    that takes it back is called a *pullback* or a *vector–Jacobian
    product*.

## Derivatives with respect to complex quantities

One wrinkle remains.  The parameters are real, and so is each
``ν_j``, but most of the quantities between them are complex.  We handle a complex
quantity ``z`` as we would a complex scalar field in a Lagrangian, by
treating ``z`` and ``\bar{z}`` as independent variables, with

```math
\frac{∂}{∂z}
= \frac{1}{2} \left( \frac{∂}{∂(\operatorname{Re} z)}
  - \mathrm{i}\, \frac{∂}{∂(\operatorname{Im} z)} \right),
\qquad
\frac{∂}{∂\bar{z}}
= \frac{1}{2} \left( \frac{∂}{∂(\operatorname{Re} z)}
  + \mathrm{i}\, \frac{∂}{∂(\operatorname{Im} z)} \right).
```

Three facts are all that we need, in which ``ν`` stands for any real
quantity.

1. Because ``ν`` is real, ``∂ν/∂\bar{z}`` is the conjugate of
   ``∂ν/∂z``.  So the one complex number ``∂ν/∂z`` holds everything
   about how ``ν`` responds to ``z``, and for complex quantities
   ``z^a``

   ```math
   dν = \frac{∂ν}{∂z^a}\, dz^a + \frac{∂ν}{∂\bar{z}^a}\, d\bar{z}^a
      = 2 \operatorname{Re} \left[ \frac{∂ν}{∂z^a}\, dz^a \right].
   ```

2. If quantities ``ζ^b`` are *holomorphic* functions of the ``z^a`` —
   depending on the ``z^a``, but not on their conjugates — then the
   chain rule has its ordinary form, with no conjugates anywhere:

   ```math
   \frac{∂ν}{∂z^a} = \frac{∂ν}{∂ζ^b}\, \frac{∂ζ^b}{∂z^a}.
   ```

   The stages ``σ``, ``ι``, and ``μ`` are of this kind.

3. If the ``z^a`` depend on a *real* quantity ``x``, then by the first
   fact

   ```math
   \frac{∂ν}{∂x}
   = 2 \operatorname{Re} \left[ \frac{∂ν}{∂z^a}\, \frac{∂z^a}{∂x} \right].
   ```

   This is how the real parameters will meet the complex harmonics.

A conjugate appears only where a function is not holomorphic, and
here the only such function is the last one, ``ν``.

!!! details "The same thing in ChainRules"

    Julia's
    [ChainRules](https://juliadiff.org/ChainRulesCore.jl/stable/maths/complex.html)
    represents the same information by the conjugate of twice our
    derivative, ``∂ν/∂(\operatorname{Re} z) + \mathrm{i}\,
    ∂ν/∂(\operatorname{Im} z)``.  With that choice, the backward step
    through multiplication by a complex number ``A`` is multiplication
    by ``\bar{A}``.  With ours it is multiplication by ``A`` itself,
    so that every backward step below uses the same array as its
    forward step, with no conjugation.

## Forward through the geometry, backward through the data

For each pixel, the chain from the parameters to ``ν_j`` is

```math
θ^p
\ \longmapsto\ \left( {}_{s}Y^{ℓm}_{j},\ g^n_j \right)
\ \longmapsto\ σ^{kh}_{j}
\ \longmapsto\ ι^{ki}_{j}
\ \longmapsto\ μ^{ki}_{j}
\ \longmapsto\ ν_j,
```

where the ``g^n_j`` skip the synthesis and enter ``ι`` and ``μ``
directly.  We have to decide where to cut it.

Forward mode all the way is what `ForwardDiff.gradient` does when it
is applied to [`transform_objective`](@ref).  Each harmonic then has
``P`` derivatives attached to it, and the synthesis has to be done for
the values and again for each derivative: ``P + 1`` syntheses for each
pixel.

Reverse mode all the way would need only two sums of that size for
each pixel: the synthesis, and one more to pull the covector back
through it.  But it would also need the backward pass of the geometry:
of the products of rotors, the recursion for the harmonics, the
supertranslation and its derivatives, and the null-rotation parameter.
Written by hand, that would be long and delicate, and it would have to
be rewritten whenever the code for the geometry changed.

Fortunately, it hardly matters how the geometry is differentiated.
Nothing in it involves the time samples: the geometry of a pixel is
computed once, however long the waveform is, and so are its
derivatives.  Forward mode costs ``P`` times very little, where
reverse mode would cost a few times very little, and dual numbers do
it with no effort on our part.  What does matter is the direction in
which we pass through the synthesis: one covector must go back through
it, rather than ``P`` vectors forward.  So we cut the chain at the
geometry.  The gradient of one pixel's ``ν_j`` is then

```math
\frac{∂ν_j}{∂θ^p}
= 2 \operatorname{Re} \left[ \sum_s
    \frac{∂ν_j}{∂\,{}_{s}Y^{ℓm}_{j}}\,
    \frac{∂\,{}_{s}Y^{ℓm}_{j}}{∂θ^p} \right]
+ \frac{∂ν_j}{∂g^n_j}\, \frac{∂g^n_j}{∂θ^p},
```

where the sum over ``s`` runs over the spin weights among the
components.  In each product, the first factor is a covector that has
been pulled back from ``ν_j``, and the second is one of the ``P``
vectors that have been pushed forward from the parameters.  So there
are three things to compute:

1. ``∂\,{}_{s}Y^{ℓm}_{j}/∂θ^p`` and ``∂g^n_j/∂θ^p``, by forward mode
   through the geometry, with dual numbers;
2. ``∂ν_j/∂\,{}_{s}Y^{ℓm}_{j}``, by a backward pass through ``ν``,
   ``μ``, ``ι``, and ``σ``, written by hand; and
3. ``∂ν_j/∂g^n_j``, which is only eight numbers, by a combination of
   the two.

For example, take ``ℓ_\mathrm{max} = 16``, so that there are 289
modes, with 4000 time samples, two components, and a supertranslation
up to ``ℓ_α = 4``, so that ``P = 31``.  One synthesis is then about
2.3 million multiplications for each pixel.  Forward mode all the way
does 32 of them.  The cut chain does two, plus work that involves no
sum of that size.

The backward pass is easy to write by hand because the stages are so
simple.  The stages ``σ`` and ``ι`` are linear in what flows through
them, and ``μ`` is linear up to an offset.  So their Jacobians are
just the arrays that the stages themselves apply:

```math
\frac{∂σ^{kh}_{j}}{∂\,{}_{s_k}Y^{ℓm}_{j}} = d^{kh}{}_{ℓm},
\qquad
\frac{∂ι^{ki}_{j}}{∂σ^{k'h}_{j}} = δ^{k}{}_{k'}\, S_{j}{}^{i}{}_{h},
\qquad
\frac{∂μ^{ki}_{j}}{∂ι^{k'i'}_{j}}
= M_{j}{}^{k}{}_{k'}(τ_{ji})\, δ^{i}{}_{i'}.
```

None of these depends on the values that flow through the stages;
they depend only on the data and on the ``g^n_j``.  So, unlike a
general computation in reverse mode, ours has nothing to record on the
way forward except its final result, the residual
``μ^{ki}_{j} - c^{ki}_{j}``.

## The backward pass

Each step takes the derivatives of ``ν_j`` with respect to the output
of a stage, and returns its derivatives with respect to the input of
that stage.

**Through ``ν``.**  Writing ``|z|^2 = z \bar{z}`` and differentiating
with ``\bar{z}`` held fixed,

```math
\frac{∂ν_j}{∂μ^{ki}_{j}}
= W_{jki} \left( \bar{μ}^{ki}_{j} - \bar{c}^{ki}_{j} \right),
```

with no sum, because the indices of the weights are labels.  This is
the only conjugate in the backward pass.  The weights act as a metric:
they turn the conjugated residual, which is a vector, into a covector.

**Through ``μ``.**

```math
\frac{∂ν_j}{∂ι^{ki}_{j}}
= \frac{∂ν_j}{∂μ^{k'i}_{j}}\, M_{j}{}^{k'}{}_{k}(τ_{ji}).
```

The mixing acts on the component index and passes the time index
along, and so does this step: for each output time, the same
coefficients are summed over their other index.  The offset
``q^{k}_{j}`` does not depend on ``ι``, so it drops out here.  Its
dependence on the geometry is not lost; it enters through
``∂ν_j/∂g^n_j``, below.

**Through ``ι``.**

```math
\frac{∂ν_j}{∂σ^{kh}_{j}}
= \frac{∂ν_j}{∂ι^{ki}_{j}}\, S_{j}{}^{i}{}_{h}.
```

The interpolation acts on the time index and passes the component
index along.  The array ``S_{j}{}^{i}{}_{h}`` is never formed.  It
would be full, since every sample influences the spline everywhere,
but the interpolation applies it in a number of operations
proportional to ``N_t + N_{t'}``, by solving a tridiagonal system and
then evaluating a cubic on each interval.  The sum over ``i`` here is
done by the same operations in the opposite order, at the same cost;
[`spline_transpose!`](@ref Scri.spline_transpose!) does it.

!!! details "How the spline is run backward"

    Consider one component, and drop its index ``k``.  Suppose that
    ``t_{h_0} ≤ τ_{ji} < t_{h_0+1}``.  The value of the spline at
    ``τ_{ji}`` is a combination of the two samples ``σ^{h_0}_{j}`` and
    ``σ^{h_0+1}_{j}``, and of the second derivatives
    ``\ddot{σ}^{h_0}_{j}`` and ``\ddot{σ}^{h_0+1}_{j}`` of the spline
    at those two times:

    ```math
    ι^{i}_{j} = E_{j}{}^{i}{}_{h}\ σ^{h}_{j}
              + F_{j}{}^{i}{}_{h}\ \ddot{σ}^{h}_{j},
    ```

    where ``E_{j}{}^{i}{}_{h}`` and ``F_{j}{}^{i}{}_{h}`` are nonzero
    only for ``h = h_0`` and ``h = h_0 + 1``.  The second derivatives
    vanish at the first and last times, which is what makes the spline
    "natural", and at the times between they solve a tridiagonal
    system,

    ```math
    Q^{h}{}_{h'}\, \ddot{σ}^{h'}_{j} = B^{h}{}_{h'}\, σ^{h'}_{j},
    ```

    where ``B^{h}{}_{h'}`` takes second differences of the samples.
    Neither ``Q`` nor ``B`` depends on the pixel.  So

    ```math
    S_{j}{}^{i}{}_{h}
    = E_{j}{}^{i}{}_{h}
    + F_{j}{}^{i}{}_{h'} \left(Q^{-1}\right)^{h'}{}_{h''} B^{h''}{}_{h},
    ```

    and the backward step is

    ```math
    \frac{∂ν_j}{∂σ^{h}_{j}}
    = \frac{∂ν_j}{∂ι^{i}_{j}}\, E_{j}{}^{i}{}_{h}
    + ζ_{jh''}\, B^{h''}{}_{h},
    \qquad \text{where} \qquad
    ζ_{jh''}\, Q^{h''}{}_{h'}
    = \frac{∂ν_j}{∂ι^{i}_{j}}\, F_{j}{}^{i}{}_{h'}.
    ```

    In words: each ``∂ν_j/∂ι^{i}_{j}`` is spread onto the two samples
    and the two second derivatives on either side of ``τ_{ji}``; the
    part on the second derivatives goes through a tridiagonal solve;
    and the result is spread onto the samples by the second
    differences.  The matrix of the coefficients ``Q^{h}{}_{h'}`` is
    symmetric.  So the system for ``ζ``, in which ``Q`` is contracted
    on its other index, has the same matrix as the forward solve, and
    reuses its factors.

    [`spline_transpose!`](@ref Scri.spline_transpose!) takes ``E`` and
    ``F`` from the same evaluation function as the forward
    interpolation, and chooses the interval of each ``τ_{ji}`` in the
    same way, including for a ``τ_{ji}`` that lands on a sample or, by
    roundoff, a hair outside the first or last one.  So the two agree
    to roundoff.

**Through ``σ``.**

```math
\frac{∂ν_j}{∂\,{}_{s}Y^{ℓm}_{j}}
= \sum_{k:\, s_k = s} \frac{∂ν_j}{∂σ^{kh}_{j}}\, d^{kh}{}_{ℓm}.
```

The sum over ``k`` collects the components that share the spin weight
``s``, and so share their harmonics.  The synthesis summed the data
over the modes, for each time and component; this step sums the same
data over the times and components, for each mode.  It is the same
number of multiplications, ``N_\mathrm{comp}\, N_t\,
N_\mathrm{modes}``, and it is the second of the two expensive
operations, and the last.

## The forward passes

**The geometry, with ``P`` derivatives.**  Each parameter is made into
a dual number that holds the value ``θ^p`` and the derivatives
``∂θ^p/∂θ^q = δ^p{}_q``.  The same functions that compute the geometry
for the objective then return it together with its derivatives:
[`bms_from_parameters`](@ref); [`aberration`](@ref), for each pixel's
rotor in the original frame; [`supertranslation_on_grid`](@ref
Scri.supertranslation_on_grid); [`compute_ðt′╱2κ`](@ref
Scri.compute_ðt′╱2κ); the conformal factor; and the calculator for the
harmonics, from SphericalFunctions.  This gives
``∂\,{}_{s}Y^{ℓm}_{j}/∂θ^p`` and ``∂g^n_j/∂θ^p``.

!!! details "Eight parameters at a time"

    A dual number that held all ``P`` derivatives would have a
    different type for each ``P``, and the code would have to be
    compiled again for every new number of parameters.  So the dual
    numbers hold eight derivatives, and the geometry is computed once
    for each group of eight parameters, with any places to spare in
    the last group left at zero.  The derivatives that come out are
    the same.  The extra work is that the values of the geometry are
    computed again for each group, which costs about as much as a few
    more derivatives.

Note that a quantity may depend on the parameters in more than one
way.  For example, ``ð^2α_j`` changes with the modes of the
supertranslation, and also with the direction at which it is
evaluated, which moves with the boost and the rotation.  Dual numbers
account for every such path without our having to list them, which is
exactly why we use them here.

**The eight numbers.**  We still need ``∂ν_j/∂g^n_j``.  By the third
fact about complex quantities,

```math
\frac{∂ν_j}{∂g^n_j}
= 2 \operatorname{Re} \left[ \frac{∂ν_j}{∂μ^{ki}_{j}}\,
  \left. \frac{∂μ^{ki}_{j}}{∂g^n_j} \right|_{σ} \right],
```

where the derivative of ``μ`` is taken through both ``ι`` and ``μ``,
with the synthesized values ``σ^{kh}_{j}`` held fixed.  The eight
numbers enter in three ways: ``κ_j^{-1}`` and ``α_j`` move the times
``τ_{ji}`` at which the spline is evaluated; ``κ_j^{-1}``,
``b^{(0)}_{j}``, and ``b^{(1)}_{j}`` set the coefficients
``M_{j}{}^{k}{}_{k'}``, which also depend on ``τ_{ji}``; and
``κ_j^{-1}`` and ``ð^2α_j`` set the offsets ``q^{k}_{j}``.  We could
differentiate each of these by hand.  But there are only eight
directions, so forward mode is cheap here too: ``ι`` and ``μ`` are
evaluated with dual numbers that hold derivatives with respect to the
eight ``g^n_j``, which gives ``μ^{ki}_{j}`` and
``∂μ^{ki}_{j}/∂g^n_j`` together.  This costs several times as much as
``ι`` and ``μ`` themselves, but it involves no sum over modes, and it
is the same whatever the number of parameters.

## Assembling the gradient

Summing over the pixels,

```math
\frac{∂Φ}{∂θ^p}
= \frac{1}{𝒩} \sum_j \frac{∂ν_j}{∂θ^p}
= \frac{1}{𝒩} \sum_j \left\{
    2 \operatorname{Re} \left[ \sum_s
      \frac{∂ν_j}{∂\,{}_{s}Y^{ℓm}_{j}}\,
      \frac{∂\,{}_{s}Y^{ℓm}_{j}}{∂θ^p} \right]
  + \frac{∂ν_j}{∂g^n_j}\, \frac{∂g^n_j}{∂θ^p}
\right\}.
```

Each pixel's term is complete in itself: its ``ν_j`` is evaluated, the
covector is pulled back to its harmonics, and the contraction with its
``P`` vectors is done, all without reference to any other pixel.  So
nothing has to be stored for the whole grid at once, and the pixels
can be shared among threads.

The diagram shows the whole computation for one pixel.  The
evaluation runs from left to right along the top, and the backward
pass returns from right to left below it.  The derivatives of the
geometry with respect to the parameters come down the left side, and
meet the covectors at the bottom.

```@raw html
<div class="gradient-flow-diagram" style="text-align: center; margin: 1.5em 0;">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1106 684" width="100%" style="max-width: 1000px;" role="img"
     aria-label="The flow of the objective and its gradient for one pixel">
  <style>
    .grad-flow text { fill: currentColor; text-anchor: middle; font-size: 17px; white-space: pre; }
    .grad-flow text.cap { font-size: 14.5px; }
    .grad-flow text.lab { font-size: 16px; }
    .grad-flow .mi { font-family: KaTeX_Math, "Times New Roman", serif; font-style: italic; }
    .grad-flow .mu { font-family: KaTeX_Main, "Times New Roman", serif; }
    .grad-flow .pr { font-family: "Times New Roman", "Liberation Serif", serif; }
    .grad-flow .box { stroke: currentColor; stroke-width: 1.2; fill: none; }
    .grad-flow .fwd { fill: rgba(56, 132, 255, 0.17); }
    .grad-flow .rev { fill: rgba(255, 140, 0, 0.22); }
    .grad-flow .plain { fill: rgba(128, 128, 128, 0.10); }
    .grad-flow .heavy { stroke-width: 3; }
    .grad-flow .edge { stroke: currentColor; stroke-width: 1.3; fill: none; marker-end: url(#grad-flow-arrow); }
    .grad-flow .edge.bare { marker-end: none; }
  </style>
  <defs>
    <marker id="grad-flow-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7"
            orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="currentColor"/>
    </marker>
  </defs>
  <g class="grad-flow">
    <path class="edge" d="M 95,8 L 95,50"/>
  <text x="104" y="26" class="lab" style="text-anchor: start"><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan></text>
  <path class="edge" d="M 180,104 L 232,104"/>
  <text x="205.0" y="93" class="lab" style="text-anchor: middle"><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan></text>
  <path class="edge" d="M 155,50 L 155,22 L 761,22 L 761,62"/>
  <text x="185" y="14" class="lab" style="text-anchor: start"><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan></text>
  <path class="edge" d="M 553,22 L 553,62"/>
  <path class="edge" d="M 450,104 L 466,104"/>
  <path class="edge" d="M 640,104 L 656,104"/>
  <path class="edge" d="M 866,104 L 882,104"/>
  <path class="edge" d="M 341,188 L 341,146"/>
  <path class="edge" d="M 341,238 L 341,280"/>
  <path class="edge" d="M 989,146 L 989,280"/>
  <text x="979" y="218" class="lab" style="text-anchor: end"><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> − </tspan><tspan
      class="mi">c</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan></text>
  <path class="edge" d="M 882,322 L 866,322"/>
  <path class="edge" d="M 656,322 L 640,322"/>
  <path class="edge" d="M 466,322 L 450,322"/>
  <path class="edge" d="M 989,364 L 989,448 L 907,448"/>
  <path class="edge" d="M 667,487 L 667,526"/>
  <path class="edge" d="M 341,364 L 341,526"/>
  <text x="331" y="453" class="lab" style="text-anchor: end"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan></text>
  <path class="edge" d="M 55,158 L 55,526"/>
  <text x="65" y="241" class="lab" style="text-anchor: start"><tspan
      class="mu">∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="10.9">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7">,</tspan></text>
  <text x="65" y="267" class="lab" style="text-anchor: start"><tspan
      class="mu">∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="2.3">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan></text>
    <g>
    <rect class="box fwd" x="10" y="50" width="170" height="108" rx="7"/>
    <text x="95" y="87" class="cap">geometry of the pixel,</text>
    <text x="95" y="109" class="cap">with its P derivatives</text>
    <text x="95" y="131"><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="10.9">,   </tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan></text>
  </g>
  <g>
    <rect class="box plain heavy" x="232" y="62" width="218" height="84" rx="7"/>
    <text x="341" y="98" class="cap">synthesis, <tspan class="mi">σ</tspan></text>
    <text x="341" y="120"><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="8.4"> = </tspan><tspan
      class="mi">d</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">ℓm</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="3" font-size="58%">k</tspan><tspan
      class="mi" dy="-8">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="466" y="62" width="174" height="84" rx="7"/>
    <text x="553" y="98" class="cap">interpolation, <tspan class="mi">ι</tspan></text>
    <text x="553" y="120"><tspan
      class="mi">ι</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> = </tspan><tspan
      class="mi">S</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"></tspan><tspan
      class="mi" dy="-7" font-size="72%">i</tspan><tspan
      class="mu" dy="7"></tspan><tspan
      class="mi" dy="5" font-size="72%">h</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">j</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="656" y="62" width="210" height="84" rx="7"/>
    <text x="761" y="98" class="cap">mixing, <tspan class="mi">μ</tspan></text>
    <text x="761" y="120"><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> = </tspan><tspan
      class="mi">M</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"></tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mu" dy="7"></tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">ι</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mi" dy="12" dx="-14.1" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="9.1"> + </tspan><tspan
      class="mi">q</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">k</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">j</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="882" y="62" width="214" height="84" rx="7"/>
    <text x="989" y="98" class="cap">norm, <tspan class="mi">ν</tspan></text>
    <text x="989" y="120"><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> = Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ik</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">W</tspan><tspan
      class="mi" dy="5" font-size="72%">jki</tspan><tspan
      class="mu" dy="-5"> |</tspan><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> − </tspan><tspan
      class="mi">c</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6">|</tspan><tspan
      class="mu" dy="-7" font-size="72%">2</tspan></text>
  </g>
  <g>
    <rect class="box io" x="266" y="188" width="150" height="50" rx="7"/>
    <text x="341" y="218"><tspan
      class="mu">data</tspan><tspan
      class="mu">  </tspan><tspan
      class="mi">d</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">ℓm</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="882" y="280" width="214" height="84" rx="7"/>
    <text x="989" y="305" class="cap">back through <tspan class="mi">ν</tspan></text>
    <text x="989" y="327"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> =</tspan></text>
    <text x="989" y="349"><tspan
      class="mi">W</tspan><tspan
      class="mi" dy="5" font-size="72%">jki</tspan><tspan
      class="mu" dy="-5"> (</tspan><tspan
      class="mi">μ</tspan><tspan
      class="mu" dy="-1.5" dx="-9.5">ˉ</tspan><tspan
      class="mi" dy="-5.5" dx="3.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> − </tspan><tspan
      class="mi">c</tspan><tspan
      class="mu" dy="-1.5" dx="-8.0">ˉ</tspan><tspan
      class="mi" dy="-5.5" dx="2.0" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6">)</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="656" y="280" width="210" height="84" rx="7"/>
    <text x="761" y="305" class="cap">back through <tspan class="mi">μ</tspan></text>
    <text x="761" y="327"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">ι</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6"> =</tspan></text>
    <text x="761" y="349"><tspan
      class="mu">(∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mi" dy="12" dx="-14.1" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="9.1">) </tspan><tspan
      class="mi">M</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"></tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="7"></tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="466" y="280" width="174" height="84" rx="7"/>
    <text x="553" y="305" class="cap">back through <tspan class="mi">ι</tspan></text>
    <text x="553" y="327"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="8.4"> =</tspan></text>
    <text x="553" y="349"><tspan
      class="mu">(∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">ι</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6">) </tspan><tspan
      class="mi">S</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"></tspan><tspan
      class="mi" dy="-7" font-size="72%">i</tspan><tspan
      class="mu" dy="7"></tspan><tspan
      class="mi" dy="5" font-size="72%">h</tspan></text>
  </g>
  <g>
    <rect class="box rev heavy" x="232" y="280" width="218" height="84" rx="7"/>
    <text x="341" y="305" class="cap">back through <tspan class="mi">σ</tspan></text>
    <text x="341" y="327"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="10.9"> =</tspan></text>
    <text x="341" y="349"><tspan
      class="mu">Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="mu" font-size="72%">:</tspan><tspan
      class="mi" font-size="72%">s</tspan><tspan
      class="mi" dy="3" font-size="58%">k</tspan><tspan
      class="mu" dy="-3" font-size="72%">=</tspan><tspan
      class="mi" font-size="72%">s</tspan><tspan
      class="mu" dy="-5"> (∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">σ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="8.4">) </tspan><tspan
      class="mi">d</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">kh</tspan><tspan
      class="mi" dy="12" dx="-13.5" font-size="72%">ℓm</tspan></text>
  </g>
  <g>
    <rect class="box fwd" x="427" y="409" width="480" height="78" rx="7"/>
    <text x="667" y="442" class="cap"><tspan class="mi">ι</tspan> and <tspan class="mi">μ</tspan>, with derivatives in the eight numbers, give ∂<tspan class="mi">μ</tspan>/∂g</text>
    <text x="667" y="464"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="2.3"> = 2 </tspan><tspan
      class="mu">Re</tspan><tspan
      class="mu"> [(∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6">) ∂</tspan><tspan
      class="mi">μ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ki</tspan><tspan
      class="mi" dy="12" dx="-10.6" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="5.6">/∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="2.3">]</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="40" y="526" width="1056" height="96" rx="7"/>
    <text x="568" y="557" class="cap">contract the covectors with the P vectors, and sum over the pixels</text>
    <text x="568" y="579"><tspan
      class="mu">∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7"> = 2 </tspan><tspan
      class="mu">Re</tspan><tspan
      class="mu"> [ Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mu" dy="-5"> (∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="10.9">) ∂</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">ℓm</tspan><tspan
      class="mi" dy="12" dx="-15.9" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="10.9">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7"> ] + (∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="2.3">) ∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mi" dy="12" dx="-7.3" font-size="72%">j</tspan><tspan
      class="mu" dy="-5" dx="2.3">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan></text>
    <text x="568" y="601"><tspan
      class="mu">∂Φ/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7"> = (1/𝒩) Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">ν</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan></text>
  </g>

  <g class="legend">
    <rect class="box fwd" x="80" y="648" width="22" height="16" rx="3"/>
    <text x="110" y="661" class="cap" style="text-anchor: start">forward mode: dual numbers</text>
    <rect class="box rev" x="356" y="648" width="22" height="16" rx="3"/>
    <text x="386" y="661" class="cap" style="text-anchor: start">backward pass, by hand</text>
    <rect class="box plain heavy" x="606" y="648" width="22" height="16" rx="3"/>
    <text x="636" y="661" class="cap" style="text-anchor: start">a sum over every mode, time, and component</text>
  </g>
  </g>
</svg>
</div>
```

## Comparing only the modes up to ``ℓ_\mathrm{fit}``

An [`ObjectivePlan`](@ref) built with `ℓ_fit = L` compares only the
modes with ``ℓ ≤ L`` of the two waveforms.  This is the usual way to
compare waveforms that were padded to a higher band limit to make room
for the transformation, because it leaves out the modes that were
added, and it also leaves out noise at higher ``ℓ``.  The residual is
projected onto the modes with ``ℓ ≤ L`` by the quadrature,

```math
Π^{ki}{}_{ℓm}
= \sum_j w_j\ \overline{{}_{s_k}Y'^{ℓm}_{j}}\,
  \left( μ^{ki}_{j} - c^{ki}_{j} \right),
```

where ``{}_{s}Y'^{ℓm}_{j}`` is the harmonic evaluated at pixel ``j``
of the grid itself, in the new frame, and so does not depend on
``θ``.  The objective is then the weighted sum of squares of these
modes:

```math
Φ = \frac{1}{𝒩} \sum_{i,k} v_i\, u_k
    \sum_{ℓm:\ |s_k| ≤ ℓ ≤ L} \left| Π^{ki}{}_{ℓm} \right|^2.
```

This is no longer a sum of terms that each belong to one pixel.  But
it depends on the transformed waveform only through ``Π``, which is
linear in it, so the facts above give

```math
𝒩\, \frac{∂Φ}{∂μ^{ki}_{j}}
= W_{jki}\ \overline{\,\sum_{ℓm} Π^{ki}{}_{ℓm}\ {}_{s_k}Y'^{ℓm}_{j}\,}.
```

That is the formula for ``∂ν_j/∂μ^{ki}_{j}`` again, with the residual
replaced by its projection onto the modes with ``ℓ ≤ L``, evaluated
back at the pixel.  The backward pass for each pixel starts from this
covector, and nothing else changes.  Because ``Π`` is needed before
any pixel can start, a first pass over all the pixels computes it, at
the cost of one more evaluation of the objective.

## Cost and memory

For each pixel, the objective costs one synthesis, and the gradient
costs two sums of that size: the synthesis, and the backward step
through it.  Nothing else in the gradient involves a sum over the
modes and the times together: not the interpolation and the mixing
with their eight derivatives, nor the backward steps through them, nor
the geometry with its ``P`` derivatives.  Only the last of these grows
with the number of parameters.  For each pixel it takes a number of
operations proportional to ``P\, N_\mathrm{modes}``, mostly for the
derivatives of the harmonics, where a synthesis takes
``N_\mathrm{comp}\, N_t\, N_\mathrm{modes}``.  So its share is
proportional to ``P / (N_\mathrm{comp}\, N_t)``, which is small
whenever the waveform has many more time samples than the
transformation has parameters.

Both the objective and the gradient take the pixels a block at a
time, so that the synthesis of a block — and, in the gradient, the
backward step through it — is one matrix product for each component.

The table gives the measured time for the gradient, and for
`ForwardDiff.gradient` applied to [`transform_objective`](@ref), as
multiples of the time for one evaluation of the objective.  The data
were the strain and ``ψ_4``, on 8 threads.  The column for
``ℓ_\mathrm{fit}`` uses a plan with ``ℓ_\mathrm{fit} =
ℓ_\mathrm{max}/2``, and is a multiple of that plan's own objective.

| ``ℓ_\mathrm{max}`` | ``N_t`` | ``P`` | gradient | gradient with ``ℓ_\mathrm{fit}`` | ForwardDiff |
| :---: | :---: | :---: | :---: | :---: | :---: |
| 8 | 400 | 7 | 4.3 | 5.4 | 16.9 |
| 8 | 400 | 31 | 7.5 | 7.2 | 68.5 |
| 8 | 4000 | 7 | 5.2 | 5.4 | 11.0 |
| 8 | 4000 | 31 | 5.1 | 5.4 | 41.9 |
| 16 | 400 | 7 | 3.8 | 4.8 | 10.8 |
| 16 | 400 | 31 | 5.6 | 5.3 | 46.1 |
| 16 | 4000 | 7 | 3.3 | 4.2 | 9.3 |
| 16 | 4000 | 31 | 3.6 | 4.3 | 40.7 |
| 32 | 400 | 7 | 2.8 | 3.9 | 9.1 |
| 32 | 400 | 31 | 3.9 | 4.8 | 42.6 |
| 32 | 4000 | 7 | 2.9 | 3.4 | 8.3 |
| 32 | 4000 | 31 | 2.8 | 3.3 | 43.8 |

So the gradient costs from about three to nearly eight evaluations of
the objective.  The multiple is smallest for large band limits and
long time series, where the synthesis is the bulk of the work.  At
small band limits, the matrix products make the synthesis so fast that
the interpolation and the mixing, with their eight derivatives, take
most of the gradient's time.  And for short time series the
derivatives of the geometry begin to show, as in the rows with ``N_t =
400`` and ``P = 31``.  They show more with more parameters: with ``P =
87``, for a supertranslation up to ``ℓ_α = 8``, the same tests cost
from three to seven evaluations with ``N_t = 4000``, and from seven to
thirteen with ``N_t = 400``.  ForwardDiff costs between one and a
little over two evaluations for each parameter.

Besides a few buffers for each thread, the storage is the harmonics
and the synthesized values ``σ^{kh}_{j}`` of one block of pixels,
which the backward pass overwrites with the derivatives of ``ν_j``
with respect to them.  A block holds as many pixels as fit in the memory
that the data themselves occupy, or in 128 MiB if that is larger, but
at least four for each thread.  The harmonics of each pixel are
evaluated first as plain numbers, for the synthesis, and again with
their derivatives, for the contraction, because keeping the
derivatives for a whole block would take ``P + 1`` times the memory.
A general-purpose tool for reverse mode would instead have to record
every operation of every synthesis, of which there are
``N_\mathrm{pix}\, N_\mathrm{comp}\, N_t\, N_\mathrm{modes}``.

### Against finite differences

Before [`transform_objective`](@ref) existed, the only way to get a
gradient was by finite differences through [`transform!`](@ref).  That
transformation chooses its own output times, which move with the
parameters, as does the one in the Python package `scri`.  So each
evaluation of an objective is a transformation, an interpolation onto
the times of the waveform it is compared with, and a sum of squares,
and a gradient by forward differences is ``P + 1`` of those.

The next table gives the time for a gradient by each method as a
multiple of the time for [`objective_and_gradient!`](@ref), which is
given in seconds.  The waveforms have 20,000 time samples and the
transformation has ``P = 31`` parameters, on 8 threads of an otherwise
idle machine.  The finite differences in Julia are through version 0.1
of this package, with the interpolation written as carefully as the
transformation itself.  Those in Python are not measured as gradients:
they are 32 times the measured time of one evaluation.

| ``ℓ_\mathrm{max}`` | components | gradient | ForwardDiff | forward differences | forward differences, Python `scri` |
| :---: | :--- | :---: | :---: | :---: | :---: |
| 8 | ``h`` | 0.11 s | 8.5 | 17 | 810 |
| 8 | ``h, ψ_4, ψ_3, ψ_2`` | 0.43 s | 8.6 | 10 | 2400 |
| 8 | ``h, ψ_4, …, ψ_0`` | 0.68 s | 8.3 | 9.4 | 1500 |
| 16 | ``h`` | 0.56 s | 10 | 14 | 460 |
| 16 | ``h, ψ_4, ψ_3, ψ_2`` | 1.9 s | 11 | 13 | 1800 |
| 16 | ``h, ψ_4, …, ψ_0`` | 3.0 s | 12 | 11 | 1100 |

The gradient itself costs from 3.8 to 6.5 evaluations of the objective
in these six cases.  With ``P = 7`` the factors are smaller, as a
count of evaluations suggests: from 2 to 3 for ForwardDiff, from 2 to
4 for forward differences, and from 90 to 420 for Python.  With ``P =
87`` they are from 20 to 31, from 26 to 43, and from 1600 to 7400.
Central differences cost twice as much as forward differences.  And an
interpolation written plainly, with one spline for each mode, makes
the finite differences in Julia three to five times slower again.

The Python package transforms the shear and all five Weyl components
together whenever any Weyl component is wanted, which is why its
factor is largest for the set with three of them.

A second set of measurements is closer to what a user of either
package would have written, and to how the gradient is usually used.
The data have the band limit ``ℓ_\mathrm{max}`` at which they are
transformed, but the objective compares only the modes up to ``ℓ =
8``: the gradient uses a plan with `ℓ_fit = 8`, and the finite
differences interpolate and sum only those modes.  Their interpolation
is written plainly, with one spline for each mode.  The waveforms
again have 20,000 time samples, on 8 threads of an otherwise idle
machine, and each measurement ran in a process of its own.  The
gradient's time is in seconds, and the others are multiples of it.  A
multiple marked ≈ is ``P + 1`` times the measured time of one
evaluation, for a gradient that would have taken more than ten
minutes.

| ``ℓ_\mathrm{max}`` | components | ``P`` | gradient | finite differences | finite differences, Python `scri` |
| :---: | :--- | :---: | :---: | :---: | :---: |
| 8 | ``h`` | 31 | 0.18 s | 35x | 500x |
| 8 | ``h`` | 87 | 0.18 s | 95x | 2300x |
| 8 | ``h, ψ_4, …, ψ_0`` | 31 | 1.1 s | 30x | ≈ 950x |
| 8 | ``h, ψ_4, …, ψ_0`` | 87 | 1.0 s | 84x | ≈ 3300x |
| 12 | ``h`` | 31 | 0.39 s | 21x | 400x |
| 12 | ``h`` | 87 | 0.42 s | 52x | ≈ 1600x |
| 12 | ``h, ψ_4, …, ψ_0`` | 31 | 2.2 s | 17x | ≈ 840x |
| 12 | ``h, ψ_4, …, ψ_0`` | 87 | 2.3 s | 46x | ≈ 2800x |
| 16 | ``h`` | 31 | 0.77 s | 15x | 340x |
| 16 | ``h`` | 87 | 0.84 s | 39x | ≈ 1100x |
| 16 | ``h, ψ_4, …, ψ_0`` | 31 | 4.4 s | 12x | ≈ 760x |
| 16 | ``h, ψ_4, …, ψ_0`` | 87 | 4.5 s | 34x | ≈ 2400x |

With one thread, the gradient takes 3.0 to 3.5 times as long.  The
finite differences in Julia take about as long as with eight, because
much of their time goes to the interpolation, which is written for one
thread, and because `transform!` in version 0.1 gains little from
threads.  Their multiples are then from 4 to 26.

The peak memory of each process was recorded as well.  With the strain
alone, it is 1.0 to 1.5 GB for the gradient, 1.1 to 1.9 GB for the
finite differences in Julia, and 1.1 to 2.8 GB for those in Python.
With the Weyl components, it is 2.5 to 4.4 GB for the gradient, 2.2 to
3.4 GB for the finite differences in Julia, and 11 to 42 GB for those
in Python.

All of these timings, and the scripts that produce them, are in the
directory `benchmark` of the repository: `comparison_2026-11-16.txt`
and `comparison/run.jl` for the first table, and
`headline_2026-11-17.txt` and `headline/run.jl` for the second.

## How it is checked

The tests in `test/gradient.jl` check each piece, and the whole.

- The backward step through ``ι`` is checked against the forward
  interpolation: for random complex ``X_i`` and ``Z^h``, the number
  ``X_i\, S_{j}{}^{i}{}_{h}\, Z^{h}`` is computed once by interpolating
  ``Z`` and once by applying [`spline_transpose!`](@ref
  Scri.spline_transpose!) to ``X``.  This is done on unevenly spaced
  times, from four samples up, with output times on the samples and at
  the ends.
- The backward step through ``σ`` is checked against the synthesis in
  the same way, including for components that share a spin weight.
- ``M_{j}{}^{k}{}_{k'}`` and ``q^{k}_{j}``, as [`mixing_matrix!`](@ref
  Scri.mixing_matrix!) extracts them, are checked against
  [`mix_components!`](@ref), on both null infinities and in more than
  one convention.
- The value and the gradient are checked against ForwardDiff applied
  to [`transform_objective`](@ref): at the identity and away from it,
  for components that mix, on both null infinities, with weights, with
  normalization, and with ``ℓ_\mathrm{fit}``.  The two gradients agree
  to a relative error of about ``10^{-15}``, and the same comparison
  passes in `Double64` with a tolerance of ``10^{-25}``.
- The gradient is checked against central differences, which is
  independent of ForwardDiff.

## The notation and the code

The function that does all of this is `gradient_core`, in
`src/gradient.jl`.  Its variables are named for the quantities on this
page, as far as Julia allows: the indices are dropped, and the
character `╱` stands in for the slash of a derivative, so that
``∂ν_j/∂μ^{ki}_{j}`` is `∂ν╱∂μ[k, i]`.

| On this page | In the code |
| :--- | :--- |
| ``d^{kh}{}_{ℓm}`` | `data[m, h, k]`, where `m` counts the modes |
| ``{}_{s}Y^{ℓm}_{j}`` | `Yb[m, q, s - s₀ + 1]`, for pixel `j = block[q]` of a block |
| ``σ^{kh}_{j}`` | `σb[h, q, k]` for a block, copied to `σⱼ[k, h]` for one pixel |
| ``μ^{ki}_{j}`` and ``∂μ^{ki}_{j}/∂g^n_j`` | the values and the partials of `μⱼ[k, i]` |
| ``τ_{ji}`` and ``M_{j}{}^{k}{}_{k'}(τ_{ji})`` | `τ[i]`, and `M` from [`mixing_matrix!`](@ref Scri.mixing_matrix!) |
| ``W_{jki}`` and ``c^{ki}_{j}`` | `w[j] * v[i] * u[k]` and `c[j, i, k]` |
| ``∂ν_j/∂μ^{ki}_{j}`` | `∂ν╱∂μ[k, i]` |
| ``∂ν_j/∂ι^{ki}_{j}`` | `∂ν╱∂ι[k, i]` |
| ``∂ν_j/∂σ^{kh}_{j}`` | `∂ν╱∂σ[k, h]` for one pixel, copied to `∂ν╱∂σb[h, q, k]` for a block |
| ``∂ν_j/∂\,{}_{s}Y^{ℓm}_{j}`` | `∂ν╱∂Yb[m, q, s - s₀ + 1]`, from [`synthesis_transpose_block!`](@ref Scri.synthesis_transpose_block!) |
| ``∂ν_j/∂g^n_j`` | `∂ν╱∂g[n]` |
| ``∂ν_j/∂θ^p`` | `∂ν╱∂θ[p, j]` |

The arrays `∂ν╱∂σb` and `∂ν╱∂Yb` are second names for the memory of
`σb` and `Yb`, which the backward pass overwrites.  The interpolated
values ``ι^{ki}_{j}`` are never stored, because [`transform_pixel!`](@ref
Scri.transform_pixel!) interpolates and mixes in one sweep.
