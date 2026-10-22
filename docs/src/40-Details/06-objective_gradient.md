```@meta
CurrentModule = Scri
```

# [The gradient of the objective](@id objective_gradient)

[`objective_and_gradient!`](@ref) returns the objective of
[`transform_objective`](@ref) together with its gradient with respect
to every parameter of the BMS transformation.  It costs about two
evaluations of the objective, however many parameters there are,
whereas forward-mode automatic differentiation costs about one
evaluation per parameter.  This page states the problem, sets up the
notation, shows the structure of the computation in a diagram, and
then derives each piece.

## The problem

A BMS transformation is fitted to data by minimizing the distance
between the transformed waveform and a waveform it is compared with.
[`transform_objective`](@ref) measures that distance as

```math
Φ(θ) = \frac{1}{𝒩} \sum_j w_j \sum_i v_i \sum_k u_k
       \left| {}^{k}d'^{i}_{j}(θ) - {}^{k}c^{i}_{j} \right|^2,
```

where ``θ`` holds the parameters of the transformation, and

- ``j`` labels the pixels of a grid on the sphere, which are the
  directions ``R'_j`` in the new frame, and the ``w_j`` are the weights
  of a quadrature on that grid.  The default grid, from
  [`sphere_quadrature`](@ref Scri.sphere_quadrature), is the product of
  Gauss–Legendre nodes in the cosine of the polar angle and equally
  spaced values of the azimuth.  With its weights, the sum over ``j``
  is exactly the integral over the sphere of the product of two fields
  band-limited to the grid's ``ℓ_\mathrm{max}``, such as
  ``|d' - c|^2`` when both are.
- ``i`` labels the output times ``t'_i``, and the ``v_i`` are the
  weights of the trapezoidal rule on them,

  ```math
  v_1 = \frac{t'_2 - t'_1}{2},
  \qquad
  v_i = \frac{t'_{i+1} - t'_{i-1}}{2},
  \qquad
  v_{N_{t'}} = \frac{t'_{N_{t'}} - t'_{N_{t'}-1}}{2},
  ```

  so that the sum over ``i`` approximates the integral over the output
  interval, even when the times are unevenly spaced.
- ``k`` labels the components of the data, and the ``u_k`` weight
  them.  The components have different units and typical sizes, and
  there is as yet no principled choice of these weights, so they are
  left to the caller, and are 1 by default.
- ``{}^{k}c^{i}_{j}`` is the waveform being compared with, which
  [`ObjectivePlan`](@ref) calls the target: its component ``k`` at
  pixel ``j`` and time ``t'_i``.
- ``𝒩`` is a normalization, either 1 or the norm of ``c``, which is the
  same sum with ``d'`` set to zero.
- ``{}^{k}d'^{i}_{j}(θ)`` is the transformed waveform at the same
  pixels and times.

Everything here except ``d'`` is fixed before the minimization starts:
the [`ObjectivePlan`](@ref) computes the grid, the weights, the values
of ``c``, and the normalization once, and the caller may choose a
different grid or other weights.  So ``Φ`` is a weighted and squared
``L^2`` distance between the two waveforms, over the unit sphere
``S^2`` and the output interval:

```math
Φ ≈ \frac{1}{𝒩} \int dt' \int_{S^2} dA\, \sum_k u_k
    \left| {}^{k}d' - {}^{k}c \right|^2.
```

Note that by approximating as we have above, we integrate over the
sphere directly in terms of the pixel values.  This allows us to skip
the expensive step of synthesizing the modes of ``|d' - c|^2``.

### What has to be computed

``Φ`` depends on ``θ`` only through the transformed values.  Splitting
each of them into its real and imaginary parts, the chain rule gives

```math
\frac{∂Φ}{∂θ^p}
= \sum_{i,j,k} \left[
    \frac{∂Φ}{∂\left(\operatorname{Re} {}^{k}d'^{i}_{j}\right)}
    \frac{∂\left(\operatorname{Re} {}^{k}d'^{i}_{j}\right)}{∂θ^p}
  + \frac{∂Φ}{∂\left(\operatorname{Im} {}^{k}d'^{i}_{j}\right)}
    \frac{∂\left(\operatorname{Im} {}^{k}d'^{i}_{j}\right)}{∂θ^p}
  \right].
```

The first factor of each term is immediate:

```math
\frac{∂Φ}{∂\left(\operatorname{Re} {}^{k}d'^{i}_{j}\right)}
= \frac{2}{𝒩}\, w_j\, v_i\, u_k
  \operatorname{Re}\left({}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j}\right),
\qquad
\frac{∂Φ}{∂\left(\operatorname{Im} {}^{k}d'^{i}_{j}\right)}
= \frac{2}{𝒩}\, w_j\, v_i\, u_k
  \operatorname{Im}\left({}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j}\right).
```

So we can rewrite the gradient as

```math
\frac{∂Φ}{∂θ^p}
= \frac{2}{𝒩}\, \sum_{i,j,k} \left[
  w_j\, v_i\, u_k\,
  \operatorname{Re}\left\{
  \overline{\left({}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j}\right)}
  \ \frac{∂\,{{}^{k}d'^{i}_{j}}}{∂θ^p}
  \right\} \right].
```

The factor ``\left({}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j}\right)`` is
computed directly when ``Φ`` is evaluated; this page is concerned more
with the next factor: how the transformed data at every pixel and time
responds to each parameter.  Computing that response one parameter at
a time, as forward-mode automatic differentiation does, costs about
one evaluation of the objective per parameter.  The rest of this page
shows how to obtain the whole sum for about two evaluations,
essentially regardless of the number of parameters, by carrying the
first factor backward through the computation of the transformed data
``d'``.

## Derivatives with respect to complex quantities

Most of the quantities between ``θ`` and ``Φ`` are complex, while
``θ`` and ``Φ`` are real.  As in
[ChainRules](https://juliadiff.org/ChainRulesCore.jl/stable/maths/complex.html),
each complex number is treated as the pair of its real and imaginary
parts, so that every derivative is an ordinary derivative of real
functions of real variables.  The ChainRulesCore documentation states
the convention this way:

> `ChainRules` follows the convention that `frule` applied to a function ``f(x + i y) = u(x,y) + i v(x,y)`` with perturbation ``Δx + i Δy`` returns the value and
>
> ```math
> \tfrac{∂u}{∂x} \, Δx + \tfrac{∂u}{∂y} \, Δy + i \, \Bigl( \tfrac{∂v}{∂x} \, Δx + \tfrac{∂v}{∂y} \, Δy \Bigr).
> ```
>
> Similarly, `rrule` applied to the same function returns the value and a pullback function which, when applied to the adjoint ``Δu + i Δv``, returns
>
> ```math
> Δu \, \tfrac{∂u}{∂x} + Δv \, \tfrac{∂v}{∂x} + i \, \Bigl(Δu \, \tfrac{∂u }{∂y} + Δv \, \tfrac{∂v}{∂y} \Bigr).
> ```
>
> If we interpret complex numbers as vectors in ``\mathbb{R}^2``, then `frule` (`rrule`) corresponds to multiplication with the (transposed) Jacobian of ``f(z)``, i.e. `frule` corresponds to
>
> ```math
> \begin{pmatrix}
> \tfrac{∂u}{∂x} \, Δx + \tfrac{∂u}{∂y} \, Δy
> \\
> \tfrac{∂v}{∂x} \, Δx + \tfrac{∂v}{∂y} \, Δy
> \end{pmatrix}
> =
> \begin{pmatrix}
> \tfrac{∂u}{∂x} & \tfrac{∂u}{∂y} \\
> \tfrac{∂v}{∂x} & \tfrac{∂v}{∂y} \\
> \end{pmatrix}
> \begin{pmatrix}
> Δx \\ Δy
> \end{pmatrix}
> ```
>
> and `rrule` corresponds to
>
> ```math
> \begin{pmatrix}
> \tfrac{∂u}{∂x} \, Δu + \tfrac{∂v}{∂x} \, Δv
> \\
> \tfrac{∂u}{∂y} \, Δu + \tfrac{∂v}{∂y} \, Δv
> \end{pmatrix}
> =
> \begin{pmatrix}
> \tfrac{∂u}{∂x} & \tfrac{∂u}{∂y} \\
> \tfrac{∂v}{∂x} & \tfrac{∂v}{∂y} \\
> \end{pmatrix}^\mathsf{T}
> \begin{pmatrix}
> Δu \\ Δv
> \end{pmatrix}.
> ```

So a tangent at a complex quantity ``z`` is a pair of real numbers,
and so is a cotangent, which maps tangents linearly to real numbers.
This page writes the cotangent of ``Φ`` at ``z`` as one complex
number,

```math
Ω_z = \frac{∂Φ}{∂(\operatorname{Re} z)}
    - \mathrm{i}\, \frac{∂Φ}{∂(\operatorname{Im} z)},
```

so that

```math
dΦ = \frac{∂Φ}{∂(\operatorname{Re} z)}\, d(\operatorname{Re} z)
   + \frac{∂Φ}{∂(\operatorname{Im} z)}\, d(\operatorname{Im} z)
   = \operatorname{Re}\left( Ω_z\, dz \right).
```

For complex quantities ``z^a``, the cotangent ``(Ω_z)_a`` is a
covector, and ``dΦ = \operatorname{Re}\left[(Ω_z)_a\, dz^a\right]``.
When ``z`` is real, ``Ω_z`` is the ordinary derivative.  The code uses
the same representation and the same names, so a variable such as
`Ωx` holds ``Ω_x``.

With this representation, the pullback of a complex-linear map is its
plain transpose.  Multiplication by a complex number ``A`` acts on the
real and imaginary parts of a tangent with the Jacobian

```math
\begin{pmatrix}
\operatorname{Re} A & -\operatorname{Im} A \\
\operatorname{Im} A & \operatorname{Re} A
\end{pmatrix},
```

whose transpose takes a cotangent ``(c_1, c_2)`` to

```math
\left(c_1 \operatorname{Re} A + c_2 \operatorname{Im} A,\
      c_2 \operatorname{Re} A - c_1 \operatorname{Im} A\right).
```

Packed as above, that is ``A\,(c_1 - \mathrm{i}\, c_2)``, which is
multiplication by ``A`` again.  So if ``z^a = A^a{}_b\, x^b``, with
coefficients that do not depend on the ``x^b``, then

```math
(Ω_x)_b = (Ω_z)_a\, A^a{}_b,
```

with no conjugation.  A conjugate appears only where a function is not
complex-linear, as the norm ``|z|^2`` is not.  And when complex
quantities ``z^a`` depend on a real quantity ``g``,

```math
\frac{∂Φ}{∂g} = \operatorname{Re}\left[ (Ω_z)_a\, \frac{∂z^a}{∂g} \right].
```

ChainRules packs the same cotangent the other way, as
``c_1 + \mathrm{i}\, c_2 = Ω_z^*``.  That number is the vector that
represents the covector under the Euclidean inner product
``\operatorname{Re}\left(z_1^*\, z_2\right)`` of the complex plane, and
with it the pullback of multiplication by ``A`` is multiplication by
``A^*``.  This page avoids that identification, so that every backward
step is a plain transpose.  Readers used to Wirtinger derivatives may
note that ``Ω_z = 2\, ∂Φ/∂z``.

## Notation

### Indices

Upper indices are slots of vectors, and lower indices are slots of
covectors.  An index repeated once up and once down is summed, as in
Einstein's convention, and every other sum is written out.  The pixel
``j`` and the spin weight ``s`` are labels: they select independent
pieces of the computation and are never contracted.  Fixed numbers,
such as the times and the weights, also carry their indices as labels.

| Index | Runs over | Number |
| :--- | :--- | :--- |
| ``j`` | pixels: the directions ``R'_j`` of the grid in the new frame | ``N_\mathrm{pix}`` |
| ``i`` | output times ``t'_i`` | ``N_{t'}`` |
| ``h, h', h''`` | input times ``t_h`` | ``N_t`` |
| ``k, k'`` | components of the data, of spin weight ``s_k`` | ``N_\mathrm{comp}`` |
| ``(ℓ,m)`` | modes, always written as one index in parentheses | ``N_\mathrm{modes}`` |
| ``n`` | the eight real numbers that describe the geometry of one pixel | ``8`` |
| ``p`` | parameters ``θ^p`` | ``P`` |

A component index is written in front of a waveform's symbol, as in
``{}^{k}x^{h}_{j}``, and a spin weight in front of a harmonic, as in
``{}_{s}Y^{(ℓ,m)}_{j}``.  The cotangent of ``Φ_j``, defined below, at a
quantity of pixel ``j`` is written with its indices outside
parentheses, as in ``\left(Ω_{x_j}\right)_{kh}``.

### The parameters

The parameters ``θ^p`` are

```math
θ = \left(\vec{r},\ \vec{η},\ a_1, …, a_{(ℓ_α+1)^2}\right),
```

made of a rotation vector ``\vec{r}``, a rapidity vector ``\vec{η}``,
and the real degrees of freedom ``a`` of the supertranslation's modes
up to ``ℓ_α``.  [`bms_from_parameters`](@ref) maps ``θ`` to the boost
velocity, the frame rotation, and the supertranslation:

```math
v⃗ = \tanh|\vec{η}|\, \frac{\vec{η}}{|\vec{η}|},
\qquad
R = \exp(\vec{r}/2),
\qquad
α_{(ℓ,m)} \text{ from the } a,
```

where the docstring of `bms_from_parameters` gives the layout of the
``a``.

### The geometry of each pixel

For each pixel ``j``, the transformation determines

- the rotor ``R_j`` of the corresponding direction in the original
  frame;
- the conformal factor ``κ_j``;
- the supertranslation and its derivatives there, ``α_j``, ``ðα_j``,
  and ``ð^2α_j``;
- the null-rotation parameter ``b = ðt'/2κ``, which is affine in the
  original time ``t``:

  ```math
  b_j(t) = b^{(0)}_j + t\, b^{(1)}_j
  ```

  (see [Computing ``ðt'/2κ``](@ref computing_eth_tprime_over_2kappa);
  ``ðα_j`` enters through ``b^{(0)}_j``).  The time is written ``t``,
  rather than as a retarded or an advanced time, because the same
  formulas serve on ``ℐ^+`` and on ``ℐ^-``.

From these, the geometry map ``Γ_j`` makes the two sets of numbers that
everything after it uses.  The first is the harmonics at the original
direction, for each spin weight ``s`` among the components:

```math
{}_{s}Y^{(ℓ,m)}_{j} = {}_{s}Y^{(ℓ,m)}(R_j).
```

The second is eight real numbers:

```math
\left(g^1_j, …, g^8_j\right)
= \left( κ_j^{-1},\ α_j,\
  \operatorname{Re} b^{(0)}_j,\ \operatorname{Im} b^{(0)}_j,\
  \operatorname{Re} b^{(1)}_j,\ \operatorname{Im} b^{(1)}_j,\
  \operatorname{Re} ð^2α_j,\ \operatorname{Im} ð^2α_j \right).
```

At pixel ``j``, the output time ``t'_i`` corresponds to the original
time

```math
τ_{ji} = κ_j^{-1}\, t'_i + α_j.
```

### The four steps

The objective is a sum over pixels, ``Φ = 𝒩^{-1} \sum_j Φ_j``, and each
``Φ_j`` is computed from the input by four functions in turn.  The
input is the modes ``{}^{k}d^{h}{}_{(ℓ,m)}`` of component ``k`` at
input time ``t_h``.

**The synthesis ``σ``** gives the value of component ``k`` in the
original direction ``R_j`` at each input time:

```math
{}^{k}x^{h}_{j} = {}^{k}d^{h}{}_{(ℓ,m)}\ {}_{s_k}Y^{(ℓ,m)}_{j}.
```

It is linear in the harmonics, and its Jacobian is made of the modes:
``\left(J_σ\right)^{kh}{}_{(ℓ,m)} = {}^{k}d^{h}{}_{(ℓ,m)}``, acting on the
harmonics of spin weight ``s_k``.  It costs
``N_t\, N_\mathrm{comp}\, N_\mathrm{modes}`` complex multiply-adds per
pixel, which dominates everything else except at the smallest band
limits.

**The interpolation ``ι``** evaluates the natural cubic spline through
those samples at the times ``τ_{ji}``:

```math
{}^{k}y^{i}_{j} = \left(S_j\right)^{i}{}_{h}\ {}^{k}x^{h}_{j}.
```

The coefficients ``\left(S_j\right)^{i}{}_{h}`` are real.  They depend
on the input times and on ``τ_{ji}``, and so on ``g_j``, but not on the
data.  With respect to ``x``, they are the Jacobian ``J_ι``.

**The mixing ``μ``** applies the transformation laws, which mix the
components:

```math
{}^{k}d'^{i}_{j}
= M^{k}{}_{k'}\!\left(g_j, τ_{ji}\right)\, {}^{k'}y^{i}_{j}
+ q^{k}\!\left(g_j\right).
```

The coefficients ``M^{k}{}_{k'}`` hold the conformal factors and the
mixing by the null rotation, so they depend on ``κ_j^{-1}`` and on
``b_j(τ_{ji})``.  On ``ℐ^+``, for example, the transformed ``ψ_3`` is
``κ^{-3}`` times ``ψ_3`` plus a multiple of ``ψ_4`` proportional to the
null-rotation parameter.  The offsets ``q^{k}`` are nonzero only for
the shear and the strain, and are proportional to ``κ_j^{-1}`` and to
``ð^2α_j`` or its conjugate.  [`mix_components!`](@ref) applies both,
and [`mixing_matrix!`](@ref Scri.mixing_matrix!) extracts
``M^{k}{}_{k'}`` from it.  Each output time is mixed separately, so the
Jacobian with respect to ``y`` is
``\left(J_μ\right)^{ki}{}_{k'i'} = M^{k}{}_{k'}\!\left(g_j, τ_{ji}\right) δ^{i}{}_{i'}``.

**The norm ``ν``** compares the result with ``c``:

```math
Φ_j = w_j \sum_i v_i \sum_k u_k
      \left| {}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j} \right|^2.
```

The weights ``w_j\, v_i\, u_k`` are the components of the metric that
defines this norm, and it is the only metric that the derivation of the
full objective needs.

Altogether, ``Φ_j = ν\left(μ_{g_j}\left(ι_{g_j}\left(σ\left(Y_j\right)\right)\right)\right)``,
where the subscripts mark the dependence of ``ι`` and ``μ`` on the
eight numbers.

## The structure of the computation

In each pixel, the expensive work is the synthesis, the step that sums
over every mode at every input time.  The parameters enter only
through the geometry, ``{}_{s}Y^{(ℓ,m)}_{j}`` and ``g^{n}_{j}``, which
is relatively cheap.  So the chain ``ν ∘ μ ∘ ι ∘ σ`` is differentiated
in reverse, pulling a cotangent back through the transpose of each
Jacobian, which costs one more pass of the same size.  The geometry is
differentiated forward, with ForwardDiff's dual numbers carrying
tangents through the Jacobian of ``Γ_j``.  The diagram shows the flow
for one pixel and omits its label ``j``.

```@raw html
<div class="gradient-flow-diagram" style="text-align: center; margin: 1.5em 0;">
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 1070 790" width="100%" style="max-width: 900px;" role="img"
     aria-label="The flow of the objective and its gradient through each direction of the grid">
  <style>
    .grad-flow text { fill: currentColor; text-anchor: middle; font-size: 17px; white-space: pre; }
    .grad-flow text.cap { font-size: 14.5px; }
    .grad-flow text.lab { font-size: 16px; }
    .grad-flow .mi { font-family: KaTeX_Math, "Times New Roman", serif; font-style: italic; }
    .grad-flow .mu { font-family: KaTeX_Main, "Times New Roman", serif; }
    .grad-flow .pr { font-family: "Times New Roman", "Liberation Serif", serif; }
    .grad-flow .box { stroke: currentColor; stroke-width: 1.2; fill: none; }
    .grad-flow .geom { fill: rgba(56, 132, 255, 0.17); }
    .grad-flow .rev { fill: rgba(255, 140, 0, 0.22); }
    .grad-flow .plain { fill: rgba(128, 128, 128, 0.10); }
    .grad-flow .heavy { stroke-width: 3; }
    .grad-flow .edge { stroke: currentColor; stroke-width: 1.3; fill: none; marker-end: url(#grad-flow-arrow); }
  </style>
  <defs>
    <marker id="grad-flow-arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7"
            orient="auto-start-reverse">
      <path d="M0,0 L10,5 L0,10 z" fill="currentColor"/>
    </marker>
  </defs>
  <g class="grad-flow">
    <path class="edge " d="M 52,54 L 70,54"/>
  <path class="edge " d="M 166,54 L 195,54"/>
  <path class="edge " d="M 222,96 L 222,139"/>
  <path class="edge " d="M 567,96 L 567,139"/>
  <path class="edge " d="M 222,197 L 222,273"/>
  <path class="edge " d="M 527,197 L 476,273"/>
  <path class="edge " d="M 607,197 L 658,273"/>
  <path class="edge " d="M 94,318 L 117,318"/>
  <path class="edge " d="M 327,318 L 347,318"/>
  <path class="edge " d="M 557,318 L 577,318"/>
  <path class="edge " d="M 787,318 L 805,318"/>
  <path class="edge " d="M 1019,318 L 1046,318"/>
  <text x="1052" y="323" class="lab" style="text-anchor: start"><tspan
      class="mi">Φ</tspan></text>
  <path class="edge " d="M 912,363 L 912,431"/>
  <text x="921" y="402.0" class="lab" style="text-anchor: start"><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="pr" dy="-7" dx="1.5" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mu" dy="7"> − </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">c</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan></text>
  <path class="edge rev" d="M 805,476 L 787,476"/>
  <path class="edge rev" d="M 577,476 L 557,476"/>
  <path class="edge rev" d="M 347,476 L 327,476"/>
  <path class="edge " d="M 892,521 L 862,562"/>
  <path class="edge " d="M 222,521 L 222,673"/>
  <path class="edge " d="M 800,638 L 800,673"/>
    <g>
    <rect class="box io" x="8" y="32" width="44" height="44" rx="7"/>
    <text x="30" y="59"><tspan
      class="mi">θ</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="70" y="32" width="96" height="44" rx="7"/>
    <text x="118" y="59"><tspan
      class="mi">v</tspan><tspan
      class="mu">⃗, </tspan><tspan
      class="mi">R</tspan><tspan
      class="mu">, </tspan><tspan
      class="mi">α</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="195" y="12" width="830" height="84" rx="7"/>
    <text x="610" y="38" class="cap">for each pixel of the new frame’s grid: the original direction, the conformal factor,</text>
    <text x="610" y="59" class="cap">the supertranslation and its derivatives there, and the null-rotation parameter</text>
    <text x="610" y="80"><tspan
      class="mi">R</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">,   </tspan><tspan
      class="mi">κ</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">,   </tspan><tspan
      class="mi">α</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">,   ð</tspan><tspan
      class="mi">α</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">,   ð</tspan><tspan
      class="mu" dy="-7" font-size="72%">2</tspan><tspan
      class="mi" dy="7">α</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5">,   </tspan><tspan
      class="mi">b</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-12" dx="-4.9" font-size="72%">(0)</tspan><tspan
      class="mu" dy="7"> + </tspan><tspan
      class="mi">t</tspan><tspan
      class="mu"> </tspan><tspan
      class="mi">b</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-12" dx="-4.9" font-size="72%">(1)</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="117" y="139" width="210" height="58" rx="7"/>
    <text x="222" y="162" class="cap">harmonics</text>
    <text x="222" y="184"><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan><tspan
      class="mu" dy="7">,   ∂(</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan><tspan
      class="mu" dy="7">)/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="457" y="139" width="220" height="58" rx="7"/>
    <text x="567" y="162" class="cap">eight real numbers</text>
    <text x="567" y="184"><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mu" dy="7">,   ∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mu" dy="7">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan></text>
  </g>
  <g>
    <rect class="box io" x="10" y="289" width="84" height="58" rx="7"/>
    <text x="52" y="312" class="cap">input</text>
    <text x="52" y="334"><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">h</tspan><tspan
      class="mu" dy="12" dx="-6.4" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan></text>
  </g>
  <g>
    <rect class="box plain heavy" x="117" y="273" width="210" height="90" rx="7"/>
    <text x="222" y="312" class="cap">1. synthesis, <tspan class="mi">σ</tspan></text>
    <text x="222" y="334"><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">x</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">h</tspan><tspan
      class="mu" dy="7"> = </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">h</tspan><tspan
      class="mu" dy="12" dx="-6.4" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="3" font-size="58%">k</tspan><tspan
      class="mi" dy="-8">Y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="347" y="273" width="210" height="90" rx="7"/>
    <text x="452" y="302" class="cap">2. interpolation, <tspan class="mi">ι</tspan></text>
    <text x="452" y="323"><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan><tspan
      class="mu" dy="7"> = </tspan><tspan
      class="mi">S</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">h</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">x</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">h</tspan></text>
    <text x="452" y="344"><tspan
      class="mi">τ</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5"> = </tspan><tspan
      class="mi">κ</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">−1</tspan><tspan
      class="mi" dy="7">t</tspan><tspan
      class="pr">′</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mu" dy="-5"> + </tspan><tspan
      class="mi">α</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="577" y="273" width="210" height="90" rx="7"/>
    <text x="682" y="302" class="cap">3. mixing, <tspan class="mi">μ</tspan></text>
    <text x="682" y="323"><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="pr" dy="-7" dx="1.5" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mu" dy="7"> = </tspan><tspan
      class="mi">M</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">k</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mi" dy="7">y</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan></text>
    <text x="682" y="344"><tspan
      class="mu">+ </tspan><tspan
      class="mi">q</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">k</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="805" y="273" width="214" height="90" rx="7"/>
    <text x="912" y="302" class="cap">4. norm, <tspan class="mi">ν</tspan></text>
    <text x="912" y="323"><tspan
      class="mi">Φ</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> = </tspan><tspan
      class="mi">w</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">ik</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">v</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mi" dy="-5">u</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan></text>
    <text x="912" y="344"><tspan
      class="mu">× |</tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="pr" dy="-7" dx="1.5" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mu" dy="7"> − </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">c</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan><tspan
      class="mu" dy="7">|</tspan><tspan
      class="mu" dy="-7" font-size="72%">2</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="805" y="431" width="214" height="90" rx="7"/>
    <text x="912" y="460" class="cap">back through <tspan class="mi">ν</tspan></text>
    <text x="912" y="481"><tspan
      class="mu">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">d</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">ki</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="912" y="502"><tspan
      class="mu">2</tspan><tspan
      class="mi">w</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mi" dy="-5">v</tspan><tspan
      class="mi" dy="5" font-size="72%">i</tspan><tspan
      class="mi" dy="-5">u</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="mu" dy="-5">(</tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="pr" dy="-7" dx="1.5" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mu" dy="7"> − </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">c</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan><tspan
      class="mu" dy="7">)</tspan><tspan
      class="mu" dy="-7" font-size="72%">*</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="577" y="431" width="210" height="90" rx="7"/>
    <text x="682" y="460" class="cap">back through <tspan class="mi">μ</tspan></text>
    <text x="682" y="481"><tspan
      class="mu">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">y</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="682" y="502"><tspan
      class="mu">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">d</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">ki</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">M</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">k</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">k</tspan><tspan
      class="pr" font-size="72%">′</tspan></text>
  </g>
  <g>
    <rect class="box rev" x="347" y="431" width="210" height="90" rx="7"/>
    <text x="452" y="460" class="cap">back through <tspan class="mi">ι</tspan></text>
    <text x="452" y="481"><tspan
      class="mu">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">x</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">kh</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="452" y="502"><tspan
      class="mu">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">y</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">ki</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi">S</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">i</tspan><tspan
      class="mi" dy="12" dx="-6.4" font-size="72%">h</tspan></text>
  </g>
  <g>
    <rect class="box rev heavy" x="117" y="431" width="210" height="90" rx="7"/>
    <text x="222" y="460" class="cap">back through <tspan class="mi">σ</tspan></text>
    <text x="222" y="481"><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mu" dy="-5">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">Y</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mu" dy="5" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan><tspan
      class="mu" dy="-5"> =</tspan></text>
    <text x="222" y="502"><tspan
      class="mu">Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">k</tspan><tspan
      class="mu" font-size="72%">:</tspan><tspan
      class="mi" font-size="72%">s</tspan><tspan
      class="mi" dy="3" font-size="58%">k</tspan><tspan
      class="mu" dy="-3" font-size="72%">=</tspan><tspan
      class="mi" font-size="72%">s</tspan><tspan
      class="mu" dy="-5"> (</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">x</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">kh</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">h</tspan><tspan
      class="mu" dy="12" dx="-6.4" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan></text>
  </g>
  <g>
    <rect class="box geom" x="580" y="562" width="440" height="76" rx="7"/>
    <text x="800" y="594" class="cap"><tspan class="mi">ι</tspan> and <tspan class="mi">μ</tspan> again, with derivatives in the eight numbers</text>
    <text x="800" y="616"><tspan
      class="mu">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">g</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">n</tspan><tspan
      class="mu" dy="-5"> = </tspan><tspan
      class="mu">Re</tspan><tspan
      class="mu"> [(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">d</tspan><tspan
      class="pr" font-size="72%">′</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">ki</tspan><tspan
      class="mu" dy="-5"> ∂(</tspan><tspan
      class="mi" dy="-7" font-size="72%">k</tspan><tspan
      class="mi" dy="7">d</tspan><tspan
      class="pr" dy="-7" dx="1.5" font-size="72%">′</tspan><tspan
      class="mi" font-size="72%">i</tspan><tspan
      class="mu" dy="7">)/∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mu" dy="7">]</tspan></text>
  </g>
  <g>
    <rect class="box plain" x="70" y="673" width="900" height="62" rx="7"/>
    <text x="520" y="698" class="cap">chain rule, summed over pixels</text>
    <text x="520" y="720"><tspan
      class="mu">∂</tspan><tspan
      class="mi">Φ</tspan><tspan
      class="mu">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7"> = (1/𝒩) Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">j</tspan><tspan
      class="mu" dy="-5"> [ Σ</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mu" dy="-5"> </tspan><tspan
      class="mu">Re</tspan><tspan
      class="mu"> [</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mu" dy="-5">(</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">Y</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mu" dy="5" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan><tspan
      class="mu" dy="-5"> ∂(</tspan><tspan
      class="mi" dy="5" font-size="72%">s</tspan><tspan
      class="mi" dy="-5">Y</tspan><tspan
      class="mu" dy="-7" dx="1.5" font-size="72%">(</tspan><tspan
      class="mi" font-size="72%">ℓ</tspan><tspan
      class="mu" font-size="72%">,</tspan><tspan
      class="mi" font-size="72%">m</tspan><tspan
      class="mu" font-size="72%">)</tspan><tspan
      class="mu" dy="7">)/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7">] + (</tspan><tspan
      class="mi">Ω</tspan><tspan
      class="mi" dy="5" font-size="72%">g</tspan><tspan
      class="mu" dy="-5">)</tspan><tspan
      class="mi" dy="5" font-size="72%">n</tspan><tspan
      class="mu" dy="-5"> ∂</tspan><tspan
      class="mi">g</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">n</tspan><tspan
      class="mu" dy="7">/∂</tspan><tspan
      class="mi">θ</tspan><tspan
      class="mi" dy="-7" dx="1.5" font-size="72%">p</tspan><tspan
      class="mu" dy="7"> ]</tspan></text>
  </g>

  <g class="legend">
    <rect class="box geom" x="110" y="758" width="22" height="16" rx="3"/>
    <text x="140" y="771" class="cap" style="text-anchor: start">forward mode: dual numbers, cheap</text>
    <rect class="box rev" x="420" y="758" width="22" height="16" rx="3"/>
    <text x="450" y="771" class="cap" style="text-anchor: start">reverse pass, by hand</text>
    <rect class="box plain heavy" x="650" y="758" width="22" height="16" rx="3"/>
    <text x="680" y="771" class="cap" style="text-anchor: start">a sum over every mode and time</text>
  </g>
  </g>
</svg>
</div>
```

## Where the parameters enter

Within pixel ``j``, ``θ`` enters only through ``{}_{s}Y^{(ℓ,m)}_{j}`` and
``g^{n}_{j}``, so the chain rule gives the covector that
`objective_and_gradient!` returns,

```math
\frac{∂Φ}{∂θ^p}
= \frac{1}{𝒩} \sum_j \left[
    \sum_s \operatorname{Re}\left(
      {}_{s}\left(Ω_{Y_j}\right)_{(ℓ,m)}\,
      \frac{∂\,{}_{s}Y^{(ℓ,m)}_{j}}{∂θ^p} \right)
  + \left(Ω_{g_j}\right)_{n}\, \frac{∂g^{n}_{j}}{∂θ^p}
\right],
```

where ``Ω_{Y_j}`` and ``Ω_{g_j}`` are cotangents of ``Φ_j``.  There are
three kinds of factor here, computed in three different ways:

- The Jacobian of the geometry, ``∂\,{}_{s}Y^{(ℓ,m)}_{j}/∂θ^p`` and
  ``∂g^{n}_{j}/∂θ^p``, comes from forward-mode differentiation in all
  the parameters.  It is cheap: a few operations per number and
  parameter.
- The cotangent of the harmonics, ``{}_{s}\left(Ω_{Y_j}\right)_{(ℓ,m)}``,
  comes from the reverse pass through ``ν``, ``μ``, ``ι``, and ``σ``,
  holding ``g_j`` fixed.
- The cotangent of the eight numbers, ``\left(Ω_{g_j}\right)_{n}``,
  comes from forward-mode differentiation of ``ι`` and ``μ`` alone, in
  those eight numbers.

Each ``Φ_j`` depends only on the quantities of pixel ``j``.  So each
pixel is computed and then differentiated at once, and nothing has to
be stored from one pixel to the next.

## The reverse pass

Each step pulls the cotangent of ``Φ_j`` back from the output of a
function to its input, through the transpose of the function's
Jacobian.

**Through the norm ``ν``:**

```math
\left(Ω_{d'_j}\right)_{ki}
= 2\, w_j\, v_i\, u_k
  \left( {}^{k}d'^{i}_{j} - {}^{k}c^{i}_{j} \right)^*.
```

Here a conjugate belongs, because ``ν`` is not complex-linear, and the
metric ``w_j\, v_i\, u_k`` lowers the indices of the conjugated
residual.

**Through the mixing ``μ``:**

```math
\left(Ω_{y_j}\right)_{k'i}
= \left(Ω_{d'_j}\right)_{ki}\, M^{k}{}_{k'}\!\left(g_j, τ_{ji}\right).
```

The offsets ``q^{k}`` do not depend on the data, so they contribute
nothing here.  Their dependence on the geometry is not lost: it enters
through ``Ω_{g_j}``, below.

**Through the interpolation ``ι``:**

```math
\left(Ω_{x_j}\right)_{kh}
= \left(Ω_{y_j}\right)_{ki} \left(S_j\right)^{i}{}_{h}.
```

This runs the interpolation backwards, so it needs the structure of
``S_j``.  For one pixel and one component, whose labels are dropped
here, suppose that ``t_{h_0} ≤ τ_{i} < t_{h_0+1}``.  The spline's value
at ``τ_i`` is a combination of the two samples ``x^{h_0}`` and
``x^{h_0+1}`` and of the spline's second derivatives ``\ddot{x}^{h_0}``
and ``\ddot{x}^{h_0+1}`` there:

```math
y^{i} = E^{i}{}_{h}\, x^{h} + F^{i}{}_{h}\, \ddot{x}^{h},
```

where ``E^{i}{}_{h}`` and ``F^{i}{}_{h}`` are nonzero only for
``h ∈ \{h_0, h_0+1\}``.  The second derivatives at the interior knots
solve a tridiagonal system, and vanish at the ends:

```math
Q^{h}{}_{h'}\, \ddot{x}^{h'} = B^{h}{}_{h'}\, x^{h'},
```

where ``B^{h}{}_{h'}`` takes second differences.  So

```math
S^{i}{}_{h}
= E^{i}{}_{h} + F^{i}{}_{h'} \left(Q^{-1}\right)^{h'}{}_{h''} B^{h''}{}_{h},
```

and the pullback of a cotangent ``(Ω_y)_i`` is

```math
(Ω_x)_{h} = (Ω_y)_{i}\, E^{i}{}_{h} + ζ_{h''}\, B^{h''}{}_{h},
\qquad
ζ_{h''}\, Q^{h''}{}_{h'} = (Ω_y)_{i}\, F^{i}{}_{h'}.
```

The matrix of ``Q`` is symmetric in this basis, so the system for
``ζ`` has the same matrix as the forward solve and reuses its factors.
That is an equality of matrix entries, not an identification of vectors
with covectors.  So each ``(Ω_y)_i`` is spread onto the samples and the
second derivatives on either side of ``τ_i``, the second-derivative
parts go through the tridiagonal solve, and the result is spread by the
second differences.  [`spline_transpose!`](@ref Scri.spline_transpose!)
takes ``E`` and ``F`` from the same evaluation function, and chooses
each time's interval in the same way, as the forward interpolation, so
the two agree to roundoff.

**Through the synthesis ``σ``:**

```math
{}_{s}\left(Ω_{Y_j}\right)_{(ℓ,m)}
= \sum_{k:\, s_k = s} \left(Ω_{x_j}\right)_{kh}\, {}^{k}d^{h}{}_{(ℓ,m)}.
```

The sum over ``k`` collects the components that share a spin weight,
and so share harmonics.  This is a sum of the same size as the
synthesis itself, and it is the second of the two expensive operations.

## The eight numbers of the geometry

The cotangent of ``Φ_j`` at the eight numbers of pixel ``j`` is

```math
\left(Ω_{g_j}\right)_{n}
= \operatorname{Re}\left[ \left(Ω_{d'_j}\right)_{ki}\,
  \left.\frac{∂\,{}^{k}d'^{i}_{j}}{∂g^{n}_{j}}\right|_{x} \right],
```

where the derivative holds the original time series ``x_j`` fixed.  The
eight numbers enter ``ι`` and ``μ`` in three ways: ``κ_j^{-1}`` and
``α_j`` move the interpolation times ``τ_{ji}``; ``κ_j^{-1}`` and the
null-rotation parameter set ``M^{k}{}_{k'}``, which also depends on
``τ_{ji}``; and ``κ_j^{-1}`` and ``ð^2α_j`` set the offsets ``q^{k}``.
Rather than differentiating each of these by hand, ``ι`` and ``μ`` are
run once more for the pixel, with dual numbers that carry the eight
tangent directions, which gives the columns
``∂\,{}^{k}d'^{i}_{j}/∂g^{n}_{j}`` of their forward Jacobian.  The
spline's second derivatives depend only on ``x_j``, so they stay plain
numbers.  This costs a small multiple of ``ι`` and ``μ``, with no sum
over modes, and it includes the offsets of the shear and the strain.

## Assembling the gradient

For the Jacobian of the geometry, ``θ`` is turned into dual numbers
whose partials carry the ``P`` tangent directions of the parameters.
The same functions that compute the objective's geometry then return it
with its derivatives: [`bms_from_parameters`](@ref), the boost and
rotation of the grid, [`supertranslation_on_grid`](@ref
Scri.supertranslation_on_grid), [`compute_ðt′╱2κ`](@ref
Scri.compute_ðt′╱2κ), the conformal factor, and the recursion for the
harmonics at ``R_j``.  These are contracted with
``{}_{s}\left(Ω_{Y_j}\right)_{(ℓ,m)}`` and ``\left(Ω_{g_j}\right)_{n}``
as in the formula above, and summed over the pixels, in parallel over
threads.

A quantity of the geometry can depend on the parameters in more than
one way at once.  ``ð^2α_j``, for example, changes with the
supertranslation's modes, and with the direction ``R_j`` at which it is
evaluated, which moves with the boost and the rotation, including the
spin phase of the rotor.  Both enter through the dual numbers, by the
product rule.

The harmonics depend only on the boost and the rotation, not on the
``a``, but they hold derivatives with respect to all the parameters.
That extra work is a few percent of the total.

## Comparing only the modes up to ``ℓ_\mathrm{fit}``

With ``ℓ_\mathrm{fit} = L``, the objective compares projections onto
the modes with ``ℓ ≤ L``, using the harmonics at the grid's own
directions, ``{}_{s}Y'^{(ℓ,m)}_{j} = {}_{s}Y^{(ℓ,m)}(R'_j)``:

```math
Φ = \frac{1}{𝒩} \sum_i v_i \sum_k u_k \sum_{(ℓ,m):\, |s_k| ≤ ℓ ≤ L}
    \left| {}^{k}Π^{i}{}_{(ℓ,m)} - {}^{k}Π_{c}{}^{i}{}_{(ℓ,m)} \right|^2,
\qquad
{}^{k}Π^{i}{}_{(ℓ,m)}
= \sum_j w_j \left({}_{s_k}Y'^{(ℓ,m)}_{j}\right)^* {}^{k}d'^{i}_{j},
```

where ``Π_{c}`` is the same projection of ``c``.  The projection uses
the inner product of the sphere, whose metric is the quadrature weights
``w_j``, and the conjugate of the harmonic; together they lower the
mode index.  ``Φ`` is not a sum of terms that each depend on one pixel.
But it depends on ``d'`` only through ``Π``, which is complex-linear in
it, so the rules above give the cotangent of ``𝒩 Φ`` at pixel ``j``,

```math
\left(Ω_{d'_j}\right)_{ki}
= 2\, w_j\, v_i\, u_k \sum_{(ℓ,m)}
  \left( {}^{k}Π^{i}{}_{(ℓ,m)} - {}^{k}Π_{c}{}^{i}{}_{(ℓ,m)} \right)^*
  \left({}_{s_k}Y'^{(ℓ,m)}_{j}\right)^*:
```

the conjugated residual of the projection, evaluated back at pixel
``j``.  Each pixel's reverse pass starts from this in place of the
cotangent through ``ν``.  Since it needs ``Π`` before any pixel can
start, a first pass over the pixels computes the projection, at the
cost of one evaluation of the objective.

## Cost and memory

Per pixel, the objective costs one synthesis and the gradient two sums
of that size, plus work that involves no sum over modes: the second run
of ``ι`` and ``μ`` with eight derivatives, the harmonics' derivatives,
and the geometry.  Forward-mode differentiation of the objective
instead carries a derivative for every parameter through the
synthesis.  Measured with strain and ``ψ_4`` data, ``N_t = 400``, and 8
threads, as multiples of the time for one evaluation of the objective:

| ``ℓ_\mathrm{max}`` | ``ℓ_\mathrm{fit}`` | ``P`` | gradient | ForwardDiff |
| :---: | :---: | :---: | :---: | :---: |
| 8 | — | 7 | 2.4 | 9.6 |
| 8 | — | 22 | 2.4 | 21.7 |
| 8 | 4 | 7 | 3.0 | 10.2 |
| 8 | 4 | 22 | 3.1 | 29.8 |
| 16 | — | 7 | 1.6 | 6.6 |
| 16 | — | 22 | 1.6 | 14.4 |
| 16 | 8 | 7 | 2.5 | 7.1 |
| 16 | 8 | 22 | 2.7 | 18.0 |

The gradient's multiple is smaller at the larger band limit, where the
sums over modes take a larger share of the time.  With
``ℓ_\mathrm{fit}``, the first pass adds about one evaluation.

No result of one pixel is needed by another, so each thread reuses a
few buffers the size of one pixel's time series and harmonics.  A
general reverse-mode tool would instead record every operation of every
synthesis, of which there are
``N_\mathrm{pix}\, N_t\, N_\mathrm{comp}\, N_\mathrm{modes}``.

## How it is checked

The tests in `test/gradient.jl` check each piece and the whole:

- ``X_i \left(S_j\right)^{i}{}_{h} Z^{h}``, computed as ``X_i``
  contracted with ``S_j Z`` and as ``X S_j`` from
  [`spline_transpose!`](@ref Scri.spline_transpose!) contracted with
  ``Z``, for random complex ``X_i`` and ``Z^h`` on uneven grids, from
  four knots up, with times on the knots and at the ends;
- ``M^{k}{}_{k'}`` and ``q^{k}`` against [`mix_components!`](@ref), for
  components that mix, on both null infinities and in more than one
  convention;
- the value and the gradient against ForwardDiff through
  [`transform_objective`](@ref), at the identity and away from it, for
  components that mix, on both null infinities, with weights,
  normalization, and ``ℓ_\mathrm{fit}``; they agree to about
  ``10^{-15}``, and to about ``10^{-25}`` in `Double64`;
- the gradient against central differences, independently of
  ForwardDiff.
